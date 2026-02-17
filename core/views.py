from __future__ import annotations

from functools import lru_cache
import re

import pandas as pd
from django.conf import settings
from django.http import HttpRequest, HttpResponse
from django.http import JsonResponse
from django.shortcuts import render
from django.utils.html import escape


DATA_PATH = settings.BASE_DIR / "data" / "interim_data" / "antidepressants_binding_affinities.csv"
ABOUT_MD_PATH = settings.BASE_DIR / "core" / "content" / "about.md"

# Keep only binary functional direction for the 2D axis.
DIRECTION_SIDE = {
    "Antagonist": -1.0,
    "Agonist": 1.0,
}

DIRECTION_COLORS = {
    "Agonist": "#0f766e",
    "Antagonist": "#1d4ed8",
}

ACTIVITY_COLORS = {
    "Agonist": "#0f766e",
    "Partial agonist": "#14b8a6",
    "Antagonist": "#1d4ed8",
    "Inhibitor": "#7c3aed",
    "Blocker": "#334155",
    "Ligand": "#b45309",
    "Unknown": "#6b7280",
}


@lru_cache(maxsize=1)
def _load_affinity_data() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH)
    df["activity_recoded"] = df["activity_recoded"].fillna("Unknown").astype(str)
    df["modelled_ki_nm"] = pd.to_numeric(df["modelled_ki_nm"], errors="coerce")
    df["pKi"] = pd.to_numeric(df["pKi"], errors="coerce")
    valid_ki = df["modelled_ki_nm"] > 0
    valid_pki = df["pKi"].notna()
    df = df.loc[valid_ki & valid_pki].copy()
    df = df.loc[df["activity_recoded"].isin(DIRECTION_SIDE)].copy()

    # Use source pKi from the dataset as the plotted affinity metric.
    df["pKi_modelled"] = df["pKi"].astype(float)
    return df


def _compass_xy(direction: str, pki: float, target: str) -> tuple[float, float]:
    pki_clamped = min(10.0, max(4.0, float(pki)))
    side = DIRECTION_SIDE.get(direction, 0.0)
    # pKi controls left/right placement within each side.
    distance = 30 + ((pki_clamped - 4.0) / 6.0) * 390.0
    x = 460 + side * distance
    # Straight-line 2D axis: no vertical encoding.
    y = 162
    return round(x, 2), round(y, 2)


def _build_view_payload(df: pd.DataFrame, selected_drug: str, selected_targets: list[str]) -> dict:
    drug_df = df.loc[df["drug"] == selected_drug].copy()
    available_targets = sorted(drug_df["target"].unique().tolist())
    selected_targets = [target for target in selected_targets if target in available_targets]
    if selected_targets:
        drug_df = drug_df.loc[drug_df["target"].isin(selected_targets)].copy()
    else:
        drug_df = drug_df.iloc[0:0].copy()

    points = []
    for row in drug_df.itertuples(index=False):
        x, y = _compass_xy(row.activity_recoded, row.pKi_modelled, row.target)
        original_activity = str(row.activity).strip() if pd.notna(row.activity) else "Unknown"
        points.append(
            {
                "x": x,
                "y": y,
                "target": row.target,
                "direction": row.activity_recoded,
                "activity": original_activity,
                "pki_modelled": round(float(row.pKi_modelled), 2),
                "color": ACTIVITY_COLORS.get(original_activity, ACTIVITY_COLORS["Unknown"]),
            }
        )

    strongest = (
        drug_df.sort_values("pKi_modelled", ascending=False)
        .head(10)[["target", "activity", "activity_recoded", "pKi_modelled"]]
        .to_dict("records")
    )

    for item in strongest:
        original_activity = str(item["activity"]).strip() if pd.notna(item["activity"]) else "Unknown"
        item["activity"] = original_activity
        item["activity_color"] = ACTIVITY_COLORS.get(original_activity, ACTIVITY_COLORS["Unknown"])
        item["pKi_modelled"] = round(float(item["pKi_modelled"]), 2)

    activity_legend = []
    for activity in sorted(drug_df["activity"].dropna().astype(str).str.strip().unique().tolist()):
        activity_legend.append(
            {"label": activity, "color": ACTIVITY_COLORS.get(activity, ACTIVITY_COLORS["Unknown"])}
        )

    return {
        "selected_drug": selected_drug,
        "available_targets": available_targets,
        "selected_targets": selected_targets,
        "points": points,
        "strongest_targets": strongest,
        "activity_legend": activity_legend,
    }


def home(request: HttpRequest) -> HttpResponse:
    df = _load_affinity_data()
    drugs = sorted(df["drug"].unique().tolist())
    drug_targets_map = {
        str(drug): sorted(group["target"].astype(str).unique().tolist())
        for drug, group in df.groupby("drug")
    }
    selected_drug = request.GET.get("drug", "")
    if selected_drug not in drugs:
        selected_drug = ""
    payload = _build_view_payload(df, selected_drug, request.GET.getlist("target"))

    context = {
        "drugs": drugs,
        "direction_colors": DIRECTION_COLORS,
        "drug_targets_map": drug_targets_map,
        "current_page": "home",
        **payload,
    }
    return render(request, "core/compass.html", context)


def _render_markdown_inline(text: str) -> str:
    content = escape(text)
    content = re.sub(r"`([^`]+)`", r"<code>\1</code>", content)
    content = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", content)
    content = re.sub(
        r"\[([^\]]+)\]\((https?://[^)]+)\)",
        r'<a href="\2" target="_blank" rel="noopener noreferrer">\1</a>',
        content,
    )
    return content


def _render_markdown_basic(markdown_text: str) -> str:
    lines = markdown_text.splitlines()
    html_parts: list[str] = []
    paragraph_lines: list[str] = []
    blockquote_lines: list[str] = []
    code_lines: list[str] = []
    in_ul = False
    in_ol = False
    in_code = False
    code_lang = ""

    def flush_paragraph() -> None:
        nonlocal paragraph_lines
        if paragraph_lines:
            paragraph = " ".join(part.strip() for part in paragraph_lines if part.strip())
            if paragraph:
                html_parts.append(f"<p>{_render_markdown_inline(paragraph)}</p>")
            paragraph_lines = []

    def flush_blockquote() -> None:
        nonlocal blockquote_lines
        if blockquote_lines:
            quote = " ".join(part.strip() for part in blockquote_lines if part.strip())
            if quote:
                html_parts.append(f"<blockquote><p>{_render_markdown_inline(quote)}</p></blockquote>")
            blockquote_lines = []

    def flush_code() -> None:
        nonlocal code_lines, code_lang
        if code_lines:
            language_class = f' class="language-{escape(code_lang)}"' if code_lang else ""
            code_text = escape("\n".join(code_lines))
            html_parts.append(f"<pre><code{language_class}>{code_text}</code></pre>")
            code_lines = []
            code_lang = ""

    def close_lists() -> None:
        nonlocal in_ul, in_ol
        if in_ul:
            html_parts.append("</ul>")
            in_ul = False
        if in_ol:
            html_parts.append("</ol>")
            in_ol = False

    for raw_line in lines:
        line = raw_line.rstrip()
        stripped = line.strip()

        if in_code:
            if stripped.startswith("```"):
                flush_code()
                in_code = False
            else:
                code_lines.append(raw_line)
            continue

        fence_match = re.match(r"^```(\w+)?\s*$", stripped)
        if fence_match:
            flush_paragraph()
            flush_blockquote()
            close_lists()
            in_code = True
            code_lang = fence_match.group(1) or ""
            continue

        if not stripped:
            flush_paragraph()
            flush_blockquote()
            close_lists()
            continue

        heading_match = re.match(r"^(#{1,3})\s+(.*)$", stripped)
        if heading_match:
            flush_paragraph()
            flush_blockquote()
            close_lists()
            level = len(heading_match.group(1))
            text = _render_markdown_inline(heading_match.group(2).strip())
            html_parts.append(f"<h{level}>{text}</h{level}>")
            continue

        if re.match(r"^(-{3,}|\*{3,}|_{3,})\s*$", stripped):
            flush_paragraph()
            flush_blockquote()
            close_lists()
            html_parts.append("<hr>")
            continue

        quote_match = re.match(r"^>\s?(.*)$", stripped)
        if quote_match:
            flush_paragraph()
            close_lists()
            blockquote_lines.append(quote_match.group(1))
            continue
        flush_blockquote()

        ordered_match = re.match(r"^\d+\.\s+(.*)$", stripped)
        if ordered_match:
            flush_paragraph()
            if in_ul:
                html_parts.append("</ul>")
                in_ul = False
            if not in_ol:
                html_parts.append("<ol>")
                in_ol = True
            html_parts.append(f"<li>{_render_markdown_inline(ordered_match.group(1).strip())}</li>")
            continue

        unordered_match = re.match(r"^[-*]\s+(.*)$", stripped)
        if unordered_match:
            flush_paragraph()
            if in_ol:
                html_parts.append("</ol>")
                in_ol = False
            if not in_ul:
                html_parts.append("<ul>")
                in_ul = True
            html_parts.append(f"<li>{_render_markdown_inline(unordered_match.group(1).strip())}</li>")
            continue

        paragraph_lines.append(stripped)

    if in_code:
        flush_code()
    flush_paragraph()
    flush_blockquote()
    close_lists()
    return "\n".join(html_parts)


def about(request: HttpRequest) -> HttpResponse:
    try:
        markdown_text = ABOUT_MD_PATH.read_text(encoding="utf-8")
    except FileNotFoundError:
        markdown_text = "# About\n\nAbout content file not found."
    about_html = _render_markdown_basic(markdown_text)
    return render(request, "core/about.html", {"current_page": "about", "about_html": about_html})


def axis_data(request: HttpRequest) -> JsonResponse:
    df = _load_affinity_data()
    drugs = sorted(df["drug"].unique().tolist())
    selected_drug = request.GET.get("drug", "")
    if selected_drug not in drugs:
        selected_drug = ""

    payload = _build_view_payload(df, selected_drug, request.GET.getlist("target"))
    return JsonResponse(payload)
