from __future__ import annotations

from functools import lru_cache
import re

import pandas as pd
from django.conf import settings
from django.http import FileResponse
from django.http import HttpRequest, HttpResponse
from django.http import JsonResponse
from django.shortcuts import render
from django.urls import reverse
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


@lru_cache(maxsize=4)
def _load_affinity_data_cached(dataset_mtime_ns: int) -> pd.DataFrame:
    # `dataset_mtime_ns` is only used as a cache key to invalidate on file changes.
    _ = dataset_mtime_ns
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


def _load_affinity_data() -> pd.DataFrame:
    return _load_affinity_data_cached(DATA_PATH.stat().st_mtime_ns)


def _compass_xy(direction: str, pki: float) -> tuple[float, float]:
    pki_clamped = min(10.0, max(4.0, float(pki)))
    side = DIRECTION_SIDE.get(direction, 0.0)
    # pKi controls left/right placement within each side.
    distance = 30 + ((pki_clamped - 4.0) / 6.0) * 390.0
    x = 460 + side * distance
    # Straight-line 2D axis: no vertical encoding.
    y = 162
    return round(x, 2), round(y, 2)


def _build_view_payload(df: pd.DataFrame, selected_drug: str, selected_receptor: str) -> dict:
    all_receptors = sorted(df["target"].unique().tolist())
    all_drugs = sorted(df["drug"].unique().tolist())
    if selected_receptor not in all_receptors:
        selected_receptor = ""

    receptor_df = df.loc[df["target"] == selected_receptor].copy() if selected_receptor else df.iloc[0:0].copy()
    available_drugs = sorted(receptor_df["drug"].unique().tolist()) if selected_receptor else []
    if selected_drug not in available_drugs:
        selected_drug = ""

    points = []
    rows = []
    activity_values: set[str] = set()
    ranking_active = bool(selected_receptor and selected_drug)

    base_pki = 0.0
    if ranking_active:
        base_row = receptor_df.loc[receptor_df["drug"] == selected_drug].head(1)
        if base_row.empty:
            ranking_active = False
            selected_drug = ""
        else:
            base_pki = float(base_row.iloc[0]["pKi_modelled"])

    if not receptor_df.empty:
        for row in receptor_df.itertuples(index=False):
            original_activity = str(row.activity).strip() if pd.notna(row.activity) else "Unknown"
            activity_values.add(original_activity)
            x, y = _compass_xy(row.activity_recoded, row.pKi_modelled)
            is_reference = bool(selected_drug) and row.drug == selected_drug
            points.append(
                {
                    "x": x,
                    "y": y,
                    "drug": row.drug,
                    "receptor": row.target,
                    "direction": row.activity_recoded,
                    "activity": original_activity,
                    "pki_modelled": round(float(row.pKi_modelled), 2),
                    "is_reference": is_reference,
                    "color": "#0f172a" if is_reference else ACTIVITY_COLORS.get(original_activity, ACTIVITY_COLORS["Unknown"]),
                }
            )

            if ranking_active and not is_reference:
                other_pki = float(row.pKi_modelled)
                delta_pki = abs(other_pki - base_pki)
                rows.append(
                    {
                        "drug": row.drug,
                        "activity": original_activity,
                        "activity_color": ACTIVITY_COLORS.get(original_activity, ACTIVITY_COLORS["Unknown"]),
                        "pKi_modelled": round(other_pki, 2),
                        "delta_pki": round(delta_pki, 2),
                    }
                )
            elif not ranking_active:
                rows.append(
                    {
                        "drug": row.drug,
                        "activity": original_activity,
                        "activity_color": ACTIVITY_COLORS.get(original_activity, ACTIVITY_COLORS["Unknown"]),
                        "pKi_modelled": round(float(row.pKi_modelled), 2),
                        "delta_pki": None,
                    }
                )

    if ranking_active:
        rows.sort(key=lambda item: item["delta_pki"], reverse=True)
    else:
        rows.sort(key=lambda item: item["pKi_modelled"], reverse=True)
    points.sort(key=lambda point: point["is_reference"])

    activity_legend = [
        {"label": activity, "color": ACTIVITY_COLORS.get(activity, ACTIVITY_COLORS["Unknown"])}
        for activity in sorted(activity_values)
    ]
    if points and ranking_active:
        activity_legend.insert(0, {"label": "Selected drug", "color": "#0f172a"})

    return {
        "selected_drug": selected_drug,
        "available_drugs": available_drugs,
        "all_drugs": all_drugs,
        "available_receptors": all_receptors,
        "selected_receptor": selected_receptor,
        "ranking_active": ranking_active,
        "points": points,
        "comparison_rows": rows,
        "activity_legend": activity_legend,
    }


def home(request: HttpRequest) -> HttpResponse:
    df = _load_affinity_data()
    receptor_drugs_map = {
        str(receptor): sorted(group["drug"].astype(str).unique().tolist())
        for receptor, group in df.groupby("target")
    }
    selected_drug = request.GET.get("drug", "")
    selected_receptor = request.GET.get("receptor") or request.GET.get("target", "")
    payload = _build_view_payload(df, selected_drug, selected_receptor)

    context = {
        "drugs": payload["available_drugs"],
        "direction_colors": DIRECTION_COLORS,
        "receptor_drugs_map": receptor_drugs_map,
        "current_page": "home",
        **payload,
    }
    return render(request, "core/compass.html", context)


def dataset_download(request: HttpRequest) -> FileResponse:
    response = FileResponse(DATA_PATH.open("rb"), content_type="text/csv")
    response["Content-Disposition"] = f'attachment; filename="{DATA_PATH.name}"'
    return response


def dataset(request: HttpRequest) -> HttpResponse:
    context = {
        "current_page": "dataset",
        "dataset_download_url": request.build_absolute_uri(reverse("dataset_download")),
        "dataset_filename": DATA_PATH.name,
    }
    return render(request, "core/dataset.html", context)


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

    selected_receptor = request.GET.get("receptor") or request.GET.get("target", "")
    payload = _build_view_payload(df, selected_drug, selected_receptor)
    return JsonResponse(payload)
