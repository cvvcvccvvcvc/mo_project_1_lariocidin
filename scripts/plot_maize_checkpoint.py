#!/usr/bin/env python3
"""Small review figure; no classification thresholds are fitted to this plot."""
import matplotlib
matplotlib.use("Agg")
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from profile_search import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--language", choices=("ru", "en"), default="ru")
    args = parser.parse_args()
    english = args.language == "en"
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    provenance = pd.read_csv(ROOT / "results/tables/maize_provenance_audit.tsv", sep="\t")
    candidates = pd.read_csv(ROOT / "results/tables/maize_curated_candidates.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.2]})
    names = ["primary_metadata_supports_metagenomic_origin", "culture_conflicts_with_MAG_label", "retrieved_origin_not_resolved"]
    counts = provenance.provenance_status.value_counts().reindex(names)
    labels = (["Metagenomic origin\nsupported", "Cell cultures despite\nMAG label", "Origin unresolved"]
              if english else ["Метагеномное происхождение\nподдерживается метаданными", "Клеточные культуры:\nпротиворечие пометке MAG", "Происхождение\nне разрешено"])
    colors = ["#327990", "#B35644", "#A8A6A0"]
    axes[0].barh(range(3), counts, color=colors, height=.58)
    axes[0].set_yticks(range(3), labels)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 310)
    axes[0].set_xlabel("Records in fixed sample" if english else "Записей в фиксированной выборке")
    axes[0].set_title("A. 399 MAG-labelled records" if english else "А. 399 записей, помеченных каталогом как MAG", loc="left", fontsize=11)
    for i, count in enumerate(counts):
        axes[0].text(count + 5, i, str(count), va="center")
    statuses = ["AsnB_favored_target_unconfirmed", "other_lasso_affinity_target_unconfirmed", "context_supported_canonical_culture_not_MAG"]
    labels = (["AsnB-favoured (37)", "Other-lasso affinity (10)", "Supported target BGC (3)\nall cell cultures"]
              if english else ["Ближе к AsnB (37)", "Ближе к другим\nlasso-системам (10)", "Целевой BGC-контекст (3)\nвсе — клеточные культуры"])
    for y, status in enumerate(statuses):
        values = candidates.loc[candidates.curation_status.eq(status), "sequence_score"].sort_values().to_numpy()
        offset = (np.arange(len(values)) % 7 - 3) * .055
        axes[1].scatter(values, y + offset, color=["#90938F", "#C4923C", "#B35644"][y], s=32, alpha=.85)
    axes[1].set_xscale("log")
    axes[1].set_xticks([80, 100, 250, 500, 1000], ["80", "100", "250", "500", "1000"])
    axes[1].set_yticks(range(3), labels)
    axes[1].invert_yaxis()
    axes[1].set_ylim(2.55, -.5)
    axes[1].axvline(80, color="#555555", linestyle=":", linewidth=1)
    axes[1].set_xlim(74, 1150)
    axes[1].set_xlabel("HMM bitscore (log scale)" if english else "HMM bitscore (логарифмическая шкала)")
    axes[1].set_title("B. All 50 hits: primary curation" if english else "Б. Все 50 кандидатов: первичная курация", loc="left", fontsize=11)
    footnote = ("The 37 + 10 hits are unconfirmed targets, not proven exclusions.\n"
                "No supported MAG hit in this sample does not prove natural absence." if english else
                "37 + 10: целевое семейство не подтверждено, но строгое исключение не доказано.\nНоль подтверждённых MAG-находок ≠ доказанное отсутствие семейства.")
    fig.text(.02, .025, footnote, fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, .12, 1, 1), w_pad=3)
    out = ROOT / "results/figures"
    out.mkdir(exist_ok=True)
    suffix = "_en" if english else ""
    fig.savefig(out / f"s3_search_checkpoint{suffix}.png", dpi=170)
    fig.savefig(out / f"s3_search_checkpoint{suffix}.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
