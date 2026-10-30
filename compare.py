"""Two result workbooks of the same images -> how much they agree.

  python compare.py <reference.xlsx> <other.xlsx>

Both are anatomy_summary_*.xlsx (from Make results). Typical uses: auto pass vs one person's review
(how good is the auto pass?), or reviewer vs reviewer on the same images (how much do careful people differ?).
Images are matched on image_file; pixel columns are compared, so a missing or different um/px can't fake a
difference. Writes compare_<other>_vs_<reference>.xlsx next to <other>.

Per trait: n images, how many changed, mean of each, mean difference (other - reference) with 95% limits of
agreement (Bland & Altman 1986, Lancet 327:307), mean absolute % difference, and Lin's concordance correlation
coefficient (Lin 1989, Biometrics 45:255; 1 = identical, it drops with both scatter and bias).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

TRAITS = {                                        # per_image column -> plain name
    "n_files": "cortical cell files (median of quarters)",
    "n_cortex_cells": "cortical cells",
    "mean_cortex_cell_area_px": "mean cortical cell area (px)",
    "cortex_width_px": "cortex width (px)",
    "root_area_px": "root area (px)",
    "stele_area_px": "stele area (px)",
    "cortex_area_px": "cortex area (px)",
    "mx_n": "metaxylem vessels",
    "mx_mean_area_px": "mean metaxylem area (px)",
    "n_air_spaces": "air spaces",
    "air_space_total_area_px": "air space area (px)",
}


def ccc(x, y):
    """Lin's concordance correlation coefficient."""
    if len(x) < 2:
        return np.nan
    sxy = np.cov(x, y, ddof=1)[0, 1]
    d = np.var(x, ddof=1) + np.var(y, ddof=1) + (x.mean() - y.mean()) ** 2
    return 2 * sxy / d if d > 0 else 1.0


def compare(ref_file, other_file):
    ref, oth = (pd.read_excel(f, "per_image") for f in (ref_file, other_file))
    m = ref.merge(oth, on="image_file", suffixes=("_ref", "_other"))
    if m.empty:
        raise SystemExit("No image_file in common: are these results of the same images?")
    rows, wide = [], m[["image_file", "status_ref", "status_other"]].copy()
    for col, name in TRAITS.items():
        if f"{col}_ref" not in m:
            continue
        x, y = m[f"{col}_ref"].astype(float), m[f"{col}_other"].astype(float)
        ok = x.notna() & y.notna()
        x, y = x[ok].to_numpy(), y[ok].to_numpy()
        d = y - x
        pct = np.abs(d) / np.where(x == 0, np.nan, np.abs(x)) * 100
        rows.append({"trait": name, "column": col, "n_images": len(x), "n_changed": int((d != 0).sum()),
                     "mean_reference": x.mean(), "mean_other": y.mean(), "mean_difference": d.mean(),
                     "loa_lower": d.mean() - 1.96 * d.std(ddof=1), "loa_upper": d.mean() + 1.96 * d.std(ddof=1),
                     "mean_abs_pct_difference": np.nanmean(pct) if len(pct) else np.nan, "lin_ccc": ccc(x, y)})
        wide[f"{col}_ref"], wide[f"{col}_other"] = m[f"{col}_ref"], m[f"{col}_other"]
        wide[f"{col}_diff"] = m[f"{col}_other"] - m[f"{col}_ref"]
    summary = pd.DataFrame(rows)
    only = sorted(set(ref.image_file) ^ set(oth.image_file))
    readme = pd.DataFrame({"what": [
        f"reference: {Path(ref_file).name}", f"other: {Path(other_file).name}",
        f"{len(m)} images matched; status in other: "
        + ", ".join(f"{k} {v}" for k, v in m.status_other.value_counts().items()),
        f"in only one of the two (left out): {', '.join(only) if only else 'none'}",
        "summary: one row per trait. difference = other - reference. loa = 95% limits of agreement "
        "(Bland & Altman 1986). mean_abs_pct_difference = mean |other - reference| / reference x 100. "
        "lin_ccc = Lin's concordance (Lin 1989): 1 identical, lower with scatter or bias",
        "per_image: reference, other and difference for every trait, one row per image",
        "With ~10 images these describe this set; they are not a test."]})
    out = Path(other_file).with_name(f"compare_{Path(other_file).stem}_vs_{Path(ref_file).stem}.xlsx")
    with pd.ExcelWriter(out) as w:
        readme.to_excel(w, sheet_name="README", index=False)
        summary.round(3).to_excel(w, sheet_name="summary", index=False)
        wide.to_excel(w, sheet_name="per_image", index=False)
    return out, summary, m


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out, summary, m = compare(sys.argv[1], sys.argv[2])
    print(summary[["trait", "n_images", "n_changed", "mean_abs_pct_difference", "lin_ccc"]]
          .round(3).to_string(index=False))
    print(f"\n{len(m)} images -> {out}")
