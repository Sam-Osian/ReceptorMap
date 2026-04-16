from __future__ import annotations

import csv
import json
import math
import shutil
import tempfile
import uuid
from functools import lru_cache
import re
from pathlib import Path

import pandas as pd
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.http import FileResponse
from django.http import HttpRequest, HttpResponse
from django.http import JsonResponse
from django.shortcuts import render
from django.urls import reverse
from django.utils.html import escape
from django.utils import timezone
from django.views.decorators.http import require_POST


DATA_PATH = settings.BASE_DIR / "data" / "antidepressants_binding_affinities.csv"
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

ACTIVITY_TO_DIRECTION = {
    "agonist": "Agonist",
    "full agonist": "Agonist",
    "partial agonist": "Agonist",
    "antagonist": "Antagonist",
    "inverse agonist": "Antagonist",
    "inhibitor": "Antagonist",
    "inhibition": "Antagonist",
    "blocker": "Antagonist",
}


def _derive_activity_recoded(activity: object) -> str:
    return ACTIVITY_TO_DIRECTION.get(str(activity or "").strip().lower(), "")


def _normalise_direction(value: object) -> str | None:
    text = str(value or "").strip().lower()
    if text == "agonist":
        return "Agonist"
    if text == "antagonist":
        return "Antagonist"
    return None


def _derive_direction(activity_recoded: object, activity: object) -> str | None:
    recoded = _normalise_direction(activity_recoded)
    if recoded:
        return recoded
    return _derive_activity_recoded(activity) or None


@lru_cache(maxsize=4)
def _load_affinity_data_cached(dataset_mtime_ns: int) -> pd.DataFrame:
    # `dataset_mtime_ns` is only used as a cache key to invalidate on file changes.
    _ = dataset_mtime_ns
    df = pd.read_csv(DATA_PATH)
    df["activity"] = df["activity"].fillna("").astype(str).str.strip()
    df["activity_recoded"] = df["activity_recoded"].fillna("").astype(str).str.strip()
    df["direction"] = [
        _derive_direction(recoded, activity)
        for recoded, activity in zip(df["activity_recoded"], df["activity"])
    ]
    df["modelled_ki_nm"] = pd.to_numeric(df["modelled_ki_nm"], errors="coerce")
    df["pKi"] = pd.to_numeric(df["pKi"], errors="coerce")
    valid_ki = df["modelled_ki_nm"] > 0
    valid_pki = df["pKi"].notna()
    df = df.loc[valid_ki & valid_pki].copy()
    df = df.loc[df["direction"].isin(DIRECTION_SIDE)].copy()

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
            original_activity = str(row.activity).strip() if pd.notna(row.activity) else ""
            if not original_activity:
                original_activity = "Unknown"
            activity_values.add(original_activity)
            x, y = _compass_xy(row.direction, row.pKi_modelled)
            is_reference = bool(selected_drug) and row.drug == selected_drug
            points.append(
                {
                    "x": x,
                    "y": y,
                    "drug": row.drug,
                    "receptor": row.target,
                    "direction": row.direction,
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


# ─────────────────────────────────────────────────────────────────────────────
# CSV Editor views
# ─────────────────────────────────────────────────────────────────────────────

FULL_CSV_COLUMNS = [
    "drug", "drug_class", "target", "activity", "activity_recoded",
    "ki_lower_nm", "ki_upper_nm", "modelled_ki_nm", "ki_is_range",
    "log10_ki_nm", "pKi", "source_url", "evidence_note", "last_verified_at",
]

EDITABLE_FIELDS = frozenset({
    "drug_class", "activity",
    "ki_lower_nm", "ki_upper_nm",
    "source_url", "evidence_note",
})

NUMERIC_FIELDS = frozenset({"ki_lower_nm", "ki_upper_nm", "modelled_ki_nm"})

ACTIVITY_OPTIONS = [
    "", "Agonist", "Antagonist", "Full agonist", "Inhibitor",
    "Partial agonist", "Blocker", "Ligand",
]


def _can_edit_csv(user) -> bool:
    return user.is_superuser or user.has_perm("core.can_edit_csv")


def _load_full_rows() -> list[dict]:
    """Load all rows from the CSV with every column preserved as strings."""
    rows: list[dict] = []
    with DATA_PATH.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            row = {col: (raw.get(col) or "").strip() for col in FULL_CSV_COLUMNS}
            rows.append(row)
    return rows


def _write_full_rows(rows: list[dict]) -> None:
    """Atomically write rows back to the CSV (temp-file + rename)."""
    tmp_fd, tmp_name = tempfile.mkstemp(
        suffix=".tmp", dir=DATA_PATH.parent, prefix="_editor_"
    )
    tmp_path = Path(tmp_name)
    try:
        with open(tmp_fd, mode="w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=FULL_CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        tmp_path.replace(DATA_PATH)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _recompute_derived(row: dict) -> None:
    """Recompute log10_ki_nm and pKi from the current modelled_ki_nm string."""
    val_str = row.get("modelled_ki_nm", "")
    try:
        ki = float(val_str)
        if ki > 0:
            log10 = math.log10(ki)
            row["log10_ki_nm"] = f"{log10:.15g}"
            row["pKi"] = f"{9.0 - log10:.15g}"
            return
    except (ValueError, TypeError):
        pass
    # Can't compute — leave existing values intact


def _get_drug_class_options() -> list[str]:
    classes: set[str] = set()
    try:
        for row in _load_full_rows():
            dc = row.get("drug_class", "")
            if dc:
                classes.add(dc)
    except Exception:
        pass
    return [""] + sorted(classes)


def _invalidate_caches() -> None:
    _load_affinity_data_cached.cache_clear()
    try:
        from core.admin import _load_dataset_rows_cached
        _load_dataset_rows_cached.cache_clear()
    except Exception:
        pass


@login_required(login_url="/editor/login/")
def editor_view(request: HttpRequest) -> HttpResponse:
    if not _can_edit_csv(request.user):
        return render(request, "core/editor_access_denied.html", status=403)

    rows = _load_full_rows()
    drug_class_options = _get_drug_class_options()
    unique_drugs = sorted({r["drug"] for r in rows})

    dropdown_options = {
        "drug_class": drug_class_options,
        "activity": ACTIVITY_OPTIONS,
    }

    context = {
        "current_page": "editor",
        "rows": rows,
        "total_rows": len(rows),
        "drug_count": len(unique_drugs),
        "dropdown_options": dropdown_options,
    }
    return render(request, "core/editor.html", context)


@login_required(login_url="/editor/login/")
@require_POST
def editor_save_api(request: HttpRequest) -> JsonResponse:
    if not _can_edit_csv(request.user):
        return JsonResponse({"error": "Permission denied"}, status=403)

    try:
        payload = json.loads(request.body)
    except (json.JSONDecodeError, ValueError):
        return JsonResponse({"error": "Invalid JSON body"}, status=400)

    change_list: list[dict] = payload.get("changes", [])
    if not change_list:
        return JsonResponse({"error": "No changes provided"}, status=400)

    # Per-field validation
    drug_class_opts = set(_get_drug_class_options())
    dropdown_validators: dict[str, set] = {
        "activity": set(ACTIVITY_OPTIONS),
        "drug_class": drug_class_opts,
    }
    errors: list[str] = []

    for i, ch in enumerate(change_list, 1):
        field = str(ch.get("field", ""))
        new_val = str(ch.get("newVal", "") or "")

        if field not in EDITABLE_FIELDS:
            errors.append(f"Change {i}: field '{field}' is not editable.")
            continue

        if field in NUMERIC_FIELDS and new_val:
            try:
                float(new_val)
            except ValueError:
                errors.append(f"Change {i}: invalid number for '{field}': {new_val!r}")

    if errors:
        return JsonResponse({"error": "Validation errors", "details": errors}, status=400)

    # Load current CSV
    all_rows = _load_full_rows()
    row_index: dict[tuple[str, str], dict] = {
        (r["drug"], r["target"]): r for r in all_rows
    }

    # Group changes by (drug, target)
    row_changes: dict[tuple[str, str], dict[str, str]] = {}
    for ch in change_list:
        drug = str(ch.get("drug", ""))
        target = str(ch.get("target", ""))
        field = str(ch.get("field", ""))
        new_val = str(ch.get("newVal", "") or "")
        key = (drug, target)
        row_changes.setdefault(key, {})[field] = new_val

    missing = [f"{d}/{t}" for d, t in row_changes if (d, t) not in row_index]
    if missing:
        return JsonResponse({"error": "Rows not found", "details": missing}, status=400)

    session_id = uuid.uuid4().hex
    audit_entries: list[dict] = []
    updated_rows: list[dict] = []
    verified_date = timezone.localdate().isoformat()

    for (drug, target), field_map in row_changes.items():
        row = row_index[(drug, target)]
        ki_touched = False
        row_changed = False

        for field, new_val in field_map.items():
            old_val = row.get(field, "")
            if old_val == new_val:
                continue
            row[field] = new_val
            row_changed = True
            audit_entries.append({
                "drug": drug, "target": target,
                "field": field, "old": old_val, "new": new_val,
            })
            if field in NUMERIC_FIELDS:
                ki_touched = True

        if ki_touched:
            # ki_is_range is derived: True only when both bounds are provided.
            lo_str = row.get("ki_lower_nm", "")
            hi_str = row.get("ki_upper_nm", "")
            is_range = bool(lo_str and hi_str)
            derived_range_val = "True" if is_range else "False"
            old_range_val = row.get("ki_is_range", "")
            if old_range_val != derived_range_val:
                row["ki_is_range"] = derived_range_val
                audit_entries.append({
                    "drug": drug, "target": target,
                    "field": "ki_is_range", "old": old_range_val, "new": derived_range_val,
                })

            if is_range and lo_str and hi_str:
                try:
                    lo, hi = float(lo_str), float(hi_str)
                    if lo > 0 and hi > 0:
                        geom = math.sqrt(lo * hi)
                        row["modelled_ki_nm"] = f"{geom:.15g}"
                except (ValueError, TypeError):
                    pass
            _recompute_derived(row)

        if "activity" in field_map:
            derived_activity_recoded = _derive_activity_recoded(row.get("activity", ""))
            old_activity_recoded = row.get("activity_recoded", "")
            if old_activity_recoded != derived_activity_recoded:
                row["activity_recoded"] = derived_activity_recoded
                row_changed = True
                audit_entries.append({
                    "drug": drug, "target": target,
                    "field": "activity_recoded", "old": old_activity_recoded, "new": derived_activity_recoded,
                })

        if row_changed:
            old_verified = row.get("last_verified_at", "")
            if old_verified != verified_date:
                row["last_verified_at"] = verified_date
                audit_entries.append({
                    "drug": drug, "target": target,
                    "field": "last_verified_at", "old": old_verified, "new": verified_date,
                })
            updated_rows.append({
                "drug": drug,
                "target": target,
                "activity_recoded": row.get("activity_recoded", ""),
                "last_verified_at": row.get("last_verified_at", ""),
            })

    try:
        _write_full_rows(all_rows)
    except Exception as exc:
        return JsonResponse({"error": f"Failed to write CSV: {exc}"}, status=500)

    _invalidate_caches()

    # Persist audit log
    from .models import CsvEditRecord  # local import avoids circular
    CsvEditRecord.objects.bulk_create([
        CsvEditRecord(
            session_id=session_id,
            drug=e["drug"],
            target=e["target"],
            field_name=e["field"],
            old_value=e["old"],
            new_value=e["new"],
            edited_by=request.user,
        )
        for e in audit_entries
    ])

    return JsonResponse({
        "success": True,
        "applied": len(audit_entries),
        "session": session_id,
        "updated_rows": updated_rows,
    })


@login_required(login_url="/editor/login/")
def editor_audit_api(request: HttpRequest) -> JsonResponse:
    if not _can_edit_csv(request.user):
        return JsonResponse({"error": "Permission denied"}, status=403)

    from .models import CsvEditRecord
    qs = CsvEditRecord.objects.select_related("edited_by").order_by("-edited_at")[:300]
    records = [
        {
            "id": r.id,
            "session": r.session_id[:8],
            "drug": r.drug,
            "target": r.target,
            "field": r.field_name,
            "old": r.old_value,
            "new": r.new_value,
            "by": r.edited_by.username if r.edited_by else "—",
            "at": r.edited_at.strftime("%Y-%m-%d %H:%M:%S UTC"),
        }
        for r in qs
    ]
    return JsonResponse({"records": records})
