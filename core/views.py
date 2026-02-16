from __future__ import annotations

from functools import lru_cache

import pandas as pd
from django.conf import settings
from django.http import HttpRequest, HttpResponse
from django.http import JsonResponse
from django.shortcuts import render


DATA_PATH = settings.BASE_DIR / "data" / "interim_data" / "antidepressants_binding_affinities.csv"

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
        **payload,
    }
    return render(request, "core/compass.html", context)


def axis_data(request: HttpRequest) -> JsonResponse:
    df = _load_affinity_data()
    drugs = sorted(df["drug"].unique().tolist())
    selected_drug = request.GET.get("drug", "")
    if selected_drug not in drugs:
        selected_drug = ""

    payload = _build_view_payload(df, selected_drug, request.GET.getlist("target"))
    return JsonResponse(payload)
