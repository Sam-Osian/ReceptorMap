from __future__ import annotations

import math
from functools import lru_cache

import pandas as pd
from django.conf import settings
from django.http import HttpRequest, HttpResponse
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


@lru_cache(maxsize=1)
def _load_affinity_data() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH)
    df["activity_recoded"] = df["activity_recoded"].fillna("Unknown").astype(str)
    df["modelled_ki_nm"] = pd.to_numeric(df["modelled_ki_nm"], errors="coerce")
    valid_ki = df["modelled_ki_nm"] > 0
    df = df.loc[valid_ki].copy()
    df = df.loc[df["activity_recoded"].isin(DIRECTION_SIDE)].copy()

    # Model pKi from modelled Ki in nM: pKi = 9 - log10(Ki[nM]).
    df["pKi_modelled"] = 9 - df["modelled_ki_nm"].astype(float).map(math.log10)
    return df


def _compass_xy(direction: str, pki: float, target: str) -> tuple[float, float]:
    pki_clamped = min(10.0, max(4.0, float(pki)))
    side = DIRECTION_SIDE.get(direction, 0.0)
    # pKi controls left/right placement within each side.
    distance = 20 + ((pki_clamped - 4.0) / 6.0) * 200.0
    x = 240 + side * distance
    # Straight-line 2D axis: no vertical encoding.
    y = 240
    return round(x, 2), round(y, 2)


def home(request: HttpRequest) -> HttpResponse:
    df = _load_affinity_data()
    drugs = sorted(df["drug"].unique().tolist())
    selected_drug = request.GET.get("drug") or (drugs[0] if drugs else "")

    drug_df = df.loc[df["drug"] == selected_drug].copy()
    points = []
    for row in drug_df.itertuples(index=False):
        x, y = _compass_xy(row.activity_recoded, row.pKi_modelled, row.target)
        points.append(
            {
                "x": x,
                "y": y,
                "target": row.target,
                "direction": row.activity_recoded,
                "activity": row.activity,
                "pki_modelled": round(float(row.pKi_modelled), 2),
                "color": DIRECTION_COLORS.get(row.activity_recoded, "#374151"),
            }
        )

    strongest = (
        drug_df.sort_values("pKi_modelled", ascending=False)
        .head(10)[["target", "activity_recoded", "pKi_modelled"]]
        .to_dict("records")
    )

    for item in strongest:
        item["pKi_modelled"] = round(float(item["pKi_modelled"]), 2)

    drug_summaries = (
        df.groupby(["drug", "drug_class"], as_index=False)
        .agg(mean_pki=("pKi_modelled", "mean"), targets=("target", "nunique"))
        .sort_values("mean_pki", ascending=False)
    )
    summary_rows = drug_summaries.to_dict("records")
    for row in summary_rows:
        row["mean_pki"] = round(float(row["mean_pki"]), 2)

    context = {
        "drugs": drugs,
        "selected_drug": selected_drug,
        "points": points,
        "strongest_targets": strongest,
        "summary_rows": summary_rows,
        "direction_colors": DIRECTION_COLORS,
    }
    return render(request, "core/compass.html", context)
