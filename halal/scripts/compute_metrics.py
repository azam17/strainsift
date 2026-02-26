#!/usr/bin/env python3
"""Compute paper metrics from benchmark TSV files.

Reads TSV files from benchmark_results/ and computes:
- MAE per experiment type
- R-squared (true vs estimated)
- Bland-Altman statistics (bias, limits of agreement)
- Sensitivity/specificity/F1 for species detection
- Ablation: MAE improvement with vs without correction
- Multi-dataset real data metrics (Denay + OPSON X)
"""

import os
import sys
import csv
from collections import defaultdict
import math

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "benchmark_results")


def load_tsv(filename):
    """Load a TSV file, return list of dicts."""
    path = os.path.join(RESULTS_DIR, filename)
    if not os.path.exists(path):
        print(f"  [SKIP] {filename} not found")
        return []
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def compute_mae(rows, nontrivial_only=True):
    """Compute Mean Absolute Error from rows with abs_error column.

    When nontrivial_only=True (default), excludes rows where both true_w and
    est_w are effectively zero, which would otherwise dilute the MAE with
    trivially-correct non-detections.
    """
    if nontrivial_only:
        errors = [float(r["abs_error"]) for r in rows
                  if "abs_error" in r
                  and (float(r.get("true_w", 0)) > 0
                       or float(r.get("est_w", 0)) > 0.001)]
    else:
        errors = [float(r["abs_error"]) for r in rows if "abs_error" in r]
    if not errors:
        return float("nan")
    return sum(errors) / len(errors)


def compute_r_squared(rows, nontrivial_only=True):
    """Compute R-squared between true_w and est_w."""
    pairs = [(float(r["true_w"]), float(r["est_w"]))
             for r in rows if "true_w" in r and "est_w" in r
             and (not nontrivial_only
                  or float(r.get("true_w", 0)) > 0
                  or float(r.get("est_w", 0)) > 0.001)]
    if len(pairs) < 2:
        return float("nan")
    mean_true = sum(t for t, _ in pairs) / len(pairs)
    ss_tot = sum((t - mean_true) ** 2 for t, _ in pairs)
    ss_res = sum((t - e) ** 2 for t, e in pairs)
    if ss_tot < 1e-15:
        return 1.0 if ss_res < 1e-15 else 0.0
    return 1.0 - ss_res / ss_tot


def bland_altman(rows, nontrivial_only=True):
    """Compute Bland-Altman statistics: bias, upper/lower limits of agreement."""
    diffs = [float(r["est_w"]) - float(r["true_w"])
             for r in rows if "true_w" in r and "est_w" in r
             and (not nontrivial_only
                  or float(r.get("true_w", 0)) > 0
                  or float(r.get("est_w", 0)) > 0.001)]
    if len(diffs) < 2:
        return float("nan"), float("nan"), float("nan")
    mean_diff = sum(diffs) / len(diffs)
    sd_diff = math.sqrt(sum((d - mean_diff) ** 2 for d in diffs) / (len(diffs) - 1))
    return mean_diff, mean_diff - 1.96 * sd_diff, mean_diff + 1.96 * sd_diff


def detection_metrics(rows, threshold=0.001):
    """Compute sensitivity, specificity, F1 for species detection.

    A species is 'truly present' if true_w > 0, and 'detected' if est_w > threshold.
    """
    tp = fp = fn = tn = 0
    for r in rows:
        true_w = float(r.get("true_w", 0))
        est_w = float(r.get("est_w", 0))
        actually_present = true_w > 0
        detected = est_w > threshold

        if actually_present and detected:
            tp += 1
        elif not actually_present and detected:
            fp += 1
        elif actually_present and not detected:
            fn += 1
        else:
            tn += 1

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    specificity = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    f1 = (2 * precision * sensitivity / (precision + sensitivity)
          if (precision + sensitivity) > 0 else float("nan"))
    return {
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "sensitivity": sensitivity, "specificity": specificity,
        "precision": precision, "F1": f1,
    }


def advanced_comparison():
    """Compare standard vs advanced inference methods side-by-side."""
    pairs = [
        ("Beef + Pork", "binary_beef_pork.tsv", "binary_beef_pork_advanced.tsv"),
        ("Beef + Horse", "binary_beef_horse.tsv", "binary_beef_horse_advanced.tsv"),
        ("Chicken + Pork", "binary_chicken_pork.tsv", "binary_chicken_pork_advanced.tsv"),
    ]

    any_found = False
    std_all = []
    adv_all = []

    for name, std_file, adv_file in pairs:
        std_rows = load_tsv(std_file)
        adv_rows = load_tsv(adv_file)
        if not std_rows or not adv_rows:
            continue

        any_found = True
        std_all.extend(std_rows)
        adv_all.extend(adv_rows)

        std_mae = compute_mae(std_rows)
        adv_mae = compute_mae(adv_rows)
        std_r2 = compute_r_squared(std_rows)
        adv_r2 = compute_r_squared(adv_rows)
        std_bias, std_lo, std_hi = bland_altman(std_rows)
        adv_bias, adv_lo, adv_hi = bland_altman(adv_rows)
        std_det = detection_metrics(std_rows)
        adv_det = detection_metrics(adv_rows)

        print(f"\n--- {name}: Standard vs Advanced ---")
        print(f"  {'Metric':<25} {'Standard':>12} {'Advanced':>12} {'Delta':>12}")
        print(f"  {'-'*25} {'-'*12} {'-'*12} {'-'*12}")
        print(f"  {'MAE':<25} {std_mae:>12.4f} {adv_mae:>12.4f} {adv_mae - std_mae:>+12.4f}")
        print(f"  {'R-squared':<25} {std_r2:>12.4f} {adv_r2:>12.4f} {adv_r2 - std_r2:>+12.4f}")
        print(f"  {'Bland-Altman bias':<25} {std_bias:>12.4f} {adv_bias:>12.4f} {adv_bias - std_bias:>+12.4f}")
        loa_width_std = std_hi - std_lo if not math.isnan(std_hi) else float("nan")
        loa_width_adv = adv_hi - adv_lo if not math.isnan(adv_hi) else float("nan")
        print(f"  {'LoA width':<25} {loa_width_std:>12.4f} {loa_width_adv:>12.4f} {loa_width_adv - loa_width_std:>+12.4f}")
        print(f"  {'Sensitivity':<25} {std_det['sensitivity']:>12.4f} {adv_det['sensitivity']:>12.4f}")
        print(f"  {'Specificity':<25} {std_det['specificity']:>12.4f} {adv_det['specificity']:>12.4f}")
        print(f"  {'F1':<25} {std_det['F1']:>12.4f} {adv_det['F1']:>12.4f}")

    if not any_found:
        return

    # Aggregate comparison
    if std_all and adv_all:
        print(f"\n{'=' * 60}")
        print("  Aggregate: Standard vs Advanced Inference")
        print(f"{'=' * 60}")

        std_mae = compute_mae(std_all)
        adv_mae = compute_mae(adv_all)
        std_r2 = compute_r_squared(std_all)
        adv_r2 = compute_r_squared(adv_all)
        std_det = detection_metrics(std_all)
        adv_det = detection_metrics(adv_all)

        print(f"  {'Metric':<25} {'Standard':>12} {'Advanced':>12} {'Winner':>12}")
        print(f"  {'-'*25} {'-'*12} {'-'*12} {'-'*12}")
        mae_winner = "Advanced" if adv_mae < std_mae else "Standard" if std_mae < adv_mae else "Tie"
        r2_winner = "Advanced" if adv_r2 > std_r2 else "Standard" if std_r2 > adv_r2 else "Tie"
        f1_winner = "Advanced" if adv_det["F1"] > std_det["F1"] else "Standard" if std_det["F1"] > adv_det["F1"] else "Tie"
        print(f"  {'MAE':<25} {std_mae:>12.4f} {adv_mae:>12.4f} {mae_winner:>12}")
        print(f"  {'R-squared':<25} {std_r2:>12.4f} {adv_r2:>12.4f} {r2_winner:>12}")
        print(f"  {'F1':<25} {std_det['F1']:>12.4f} {adv_det['F1']:>12.4f} {f1_winner:>12}")


def main():
    print("=" * 60)
    print("  SpeciesID Paper Metrics")
    print("=" * 60)

    if not os.path.isdir(RESULTS_DIR):
        print(f"Error: {RESULTS_DIR} not found. Run run_paper_benchmark.sh first.")
        sys.exit(1)

    # --- Binary mixtures ---
    for name, filename in [
        ("Beef + Pork", "binary_beef_pork.tsv"),
        ("Beef + Horse", "binary_beef_horse.tsv"),
        ("Chicken + Pork", "binary_chicken_pork.tsv"),
    ]:
        rows = load_tsv(filename)
        if not rows:
            continue
        mae = compute_mae(rows)
        r2 = compute_r_squared(rows)
        bias, loa_lo, loa_hi = bland_altman(rows)
        det = detection_metrics(rows)

        print(f"\n--- {name} ---")
        print(f"  MAE:          {mae:.4f} ({mae*100:.2f} pp)")
        print(f"  R-squared:    {r2:.4f}")
        print(f"  Bland-Altman: bias={bias:.4f}, LoA=[{loa_lo:.4f}, {loa_hi:.4f}]")
        print(f"  Detection:    TP={det['TP']}, FP={det['FP']}, "
              f"FN={det['FN']}, TN={det['TN']}")
        print(f"  Sensitivity:  {det['sensitivity']:.4f}")
        print(f"  Specificity:  {det['specificity']:.4f}")
        print(f"  F1:           {det['F1']:.4f}")

        # Per-ratio MAE
        by_ratio = defaultdict(list)
        for r in rows:
            by_ratio[r["experiment"]].append(float(r["abs_error"]))
        print("  Per-ratio MAE:")
        for exp, errs in sorted(by_ratio.items()):
            print(f"    {exp}: {sum(errs)/len(errs):.4f}")

    # --- Ternary ---
    rows = load_tsv("ternary.tsv")
    if rows:
        print(f"\n--- Ternary Mixtures ---")
        mae = compute_mae(rows)
        r2 = compute_r_squared(rows)
        print(f"  MAE:       {mae:.4f} ({mae*100:.2f} pp)")
        print(f"  R-squared: {r2:.4f}")

    # --- Trace detection ---
    rows = load_tsv("trace_detection.tsv")
    if rows:
        print(f"\n--- Trace Detection ---")
        by_depth = defaultdict(list)
        for r in rows:
            exp = r["experiment"]
            # Extract trace fraction and read depth
            by_depth[exp].append(r)

        for exp_name in sorted(by_depth.keys()):
            exp_rows = by_depth[exp_name]
            det = detection_metrics(exp_rows, threshold=0.001)
            mae = compute_mae(exp_rows)
            print(f"  {exp_name}: MAE={mae:.4f}, "
                  f"Sensitivity={det['sensitivity']:.2f}, "
                  f"TP={det['TP']}, FN={det['FN']}")

    # --- Ablation: degradation ---
    rows = load_tsv("ablation_degradation.tsv")
    if rows:
        print(f"\n--- Ablation: Degradation Correction ---")
        with_deg = [r for r in rows if r.get("mode", "").strip() == "with_degradation"]
        without_deg = [r for r in rows if r.get("mode", "").strip() == "no_degradation"]
        if with_deg and without_deg:
            mae_with = compute_mae(with_deg)
            mae_without = compute_mae(without_deg)
            improvement = ((mae_without - mae_with) / mae_without * 100
                           if mae_without > 0 else 0)
            print(f"  MAE without degradation: {mae_without:.4f}")
            print(f"  MAE with degradation:    {mae_with:.4f}")
            print(f"  Improvement:             {improvement:.1f}%")

    # --- Performance ---
    rows = load_tsv("performance.tsv")
    if rows:
        print(f"\n--- Computational Performance ---")
        by_depth = defaultdict(list)
        for r in rows:
            rpm = int(r["n_reads_per_marker"])
            by_depth[rpm].append(float(r["wall_seconds"]))
        print(f"  {'Reads/marker':>15}  {'Mean time (s)':>15}  {'Std (s)':>10}")
        for rpm in sorted(by_depth.keys()):
            times = by_depth[rpm]
            mean_t = sum(times) / len(times)
            if len(times) > 1:
                std_t = math.sqrt(sum((t - mean_t)**2 for t in times) / (len(times) - 1))
            else:
                std_t = 0
            total_reads = rpm * 3  # 3 markers
            print(f"  {total_reads:>15}  {mean_t:>15.3f}  {std_t:>10.3f}")

    # --- Advanced vs Standard comparison ---
    advanced_comparison()

    # --- Summary ---
    print(f"\n{'=' * 60}")
    print("  Summary for Paper")
    print(f"{'=' * 60}")

    all_binary = []
    for fn in ["binary_beef_pork.tsv", "binary_beef_horse.tsv", "binary_chicken_pork.tsv"]:
        all_binary.extend(load_tsv(fn))
    if all_binary:
        overall_mae = compute_mae(all_binary)
        overall_r2 = compute_r_squared(all_binary)
        overall_det = detection_metrics(all_binary)
        print(f"  Overall binary MAE:    {overall_mae:.4f} ({overall_mae*100:.2f} pp)")
        print(f"  Overall binary R^2:    {overall_r2:.4f}")
        print(f"  Overall detection F1:  {overall_det['F1']:.4f}")
        print(f"  Detection: {overall_det['TP']} TP, {overall_det['FP']} FP, "
              f"{overall_det['FN']} FN, {overall_det['TN']} TN")

    # --- Real data: multi-dataset ---
    real_data_metrics_multi()


def _dataset_metrics(rows, details, dataset_label):
    """Compute and print metrics for a single real-data dataset."""
    if not rows:
        return {}

    # Overview
    n_samples = len(details) if details else 0
    total_reads = sum(int(d.get("total_reads", 0)) for d in details) if details else 0
    classified_reads = sum(int(d.get("classified_reads", 0)) for d in details) if details else 0
    classified_pct = classified_reads / total_reads * 100 if total_reads > 0 else 0

    print(f"\n  Samples analyzed:     {n_samples}")
    print(f"  Total reads:          {total_reads:,}")
    print(f"  Classified reads:     {classified_reads:,} ({classified_pct:.1f}%)")

    # Per-category breakdown
    categories = defaultdict(list)
    for r in rows:
        categories[r["category"]].append(r)

    cat_labels = {
        "spike": "Spike (known single species, in DB)",
        "spike_outofdb": "Spike (species NOT in DB)",
        "spike_unknown": "Spike (numbered, species TBD)",
        "lgc": "LGC Certified Reference Materials",
        "equiden": "Equiden (equine mixtures)",
        "gemisch": "Gemisch (multi-species mixtures)",
        "proficiency": "Proficiency Test Samples (DLA/LVU)",
        "exoten": "Exotic Species Samples",
        "lippold": "Lippold (boiled sausage, multi-species)",
        "real_product": "Real Market Products",
        "mock": "Mock Mixtures (known composition)",
    }

    cat_summary = {}
    for cat in ["spike", "spike_outofdb", "spike_unknown", "lgc", "equiden",
                 "gemisch", "proficiency", "exoten", "lippold", "real_product", "mock"]:
        cat_rows = categories.get(cat, [])
        if not cat_rows:
            continue

        print(f"\n--- {cat_labels.get(cat, cat)} ---")

        # Group by sample
        by_sample = defaultdict(list)
        for r in cat_rows:
            by_sample[r["accession"]].append(r)

        n_cat_samples = len(by_sample)
        cat_class_rates = []

        for acc in sorted(by_sample.keys()):
            sample_rows = by_sample[acc]
            sample_name = sample_rows[0]["sample_name"]
            expected = sample_rows[0]["expected"]
            species_list = [(r["species"], float(r["weight_pct"])) for r in sample_rows]
            species_list.sort(key=lambda x: -x[1])

            top_species = ", ".join(f"{sp} ({w:.1f}%)" for sp, w in species_list[:3])
            n_total = int(sample_rows[0].get("total_reads", 0))
            n_class = int(sample_rows[0].get("classified_reads", 0))
            class_pct = n_class / n_total * 100 if n_total > 0 else 0
            cat_class_rates.append(class_pct)
            print(f"  {acc} ({sample_name}): {n_class}/{n_total} classified ({class_pct:.0f}%)")
            print(f"    Expected: {expected}")
            print(f"    Detected: {top_species}")

        if cat_class_rates:
            min_rate = min(cat_class_rates)
            max_rate = max(cat_class_rates)
            mean_rate = sum(cat_class_rates) / len(cat_class_rates)
            cat_summary[cat] = {
                "n": n_cat_samples,
                "class_min": min_rate,
                "class_max": max_rate,
                "class_mean": mean_rate,
            }

    # Detection metrics for spike samples
    spike_rows = [r for r in rows if r["category"] == "spike"]
    tp = fp = fn = 0
    if spike_rows:
        print(f"\n--- Spike Detection Metrics ---")

        by_sample = defaultdict(list)
        for r in spike_rows:
            by_sample[r["accession"]].append(r)

        for acc, sample_rows in by_sample.items():
            expected_str = sample_rows[0]["expected"]
            expected_set = set(expected_str.split(",")) if expected_str not in ("unknown", "mixture") else set()
            detected_set = set(r["species"] for r in sample_rows)

            for sp in expected_set:
                if sp in detected_set:
                    tp += 1
                else:
                    fn += 1
            for sp in detected_set:
                if sp not in expected_set:
                    fp += 1

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
        f1 = (2 * precision * sensitivity / (precision + sensitivity)
              if (precision + sensitivity) > 0 else float("nan"))

        print(f"  True Positives:  {tp}")
        print(f"  False Positives: {fp}")
        print(f"  False Negatives: {fn}")
        print(f"  Sensitivity:     {sensitivity:.3f}")
        print(f"  Precision:       {precision:.3f}")
        print(f"  F1 Score:        {f1:.3f}")

    # Out-of-database samples
    ood_rows = [r for r in rows if r["category"] == "spike_outofdb"]
    if ood_rows:
        print(f"\n--- Out-of-Database Species ---")
        by_sample = defaultdict(list)
        for r in ood_rows:
            by_sample[r["accession"]].append(r)
        for acc, sample_rows in by_sample.items():
            sample_name = sample_rows[0]["sample_name"]
            expected = sample_rows[0]["expected"]
            n_class = int(sample_rows[0].get("classified_reads", 0))
            n_total = int(sample_rows[0].get("total_reads", 0))
            class_pct = n_class / n_total * 100 if n_total > 0 else 0
            species_list = [(r["species"], float(r["weight_pct"])) for r in sample_rows]
            species_list.sort(key=lambda x: -x[1])
            top = ", ".join(f"{sp} ({w:.1f}%)" for sp, w in species_list[:3])
            print(f"  {acc} ({sample_name}): true={expected}")
            print(f"    Classified: {n_class}/{n_total} ({class_pct:.0f}%), detected as: {top}")

    return {
        "n_samples": n_samples,
        "total_reads": total_reads,
        "classified_reads": classified_reads,
        "classified_pct": classified_pct,
        "spike_tp": tp,
        "spike_fp": fp,
        "spike_fn": fn,
        "cat_summary": cat_summary,
    }


def real_data_metrics_multi():
    """Compute detection metrics from multiple real data benchmark datasets."""
    datasets = [
        ("denay", "Denay et al. (2023) - PRJEB57117"),
        ("opsonx", "OPSON X (Kappel et al. 2023) - PRJNA926813"),
    ]

    # Also try legacy filenames for backward compatibility
    legacy_rows = load_tsv("real_data.tsv")
    legacy_details = load_tsv("real_data_details.tsv")

    all_stats = []
    any_found = False

    for tag, label in datasets:
        rows = load_tsv(f"real_data_{tag}.tsv")
        details = load_tsv(f"real_data_{tag}_details.tsv")
        if not rows:
            continue

        any_found = True
        print(f"\n{'=' * 60}")
        print(f"  Real Data: {label}")
        print(f"{'=' * 60}")

        stats = _dataset_metrics(rows, details, label)
        stats["label"] = label
        stats["tag"] = tag
        all_stats.append(stats)

    # Fall back to legacy single-dataset format
    if not any_found and legacy_rows:
        print(f"\n{'=' * 60}")
        print("  Real Data: Denay et al. (2023)")
        print(f"{'=' * 60}")
        stats = _dataset_metrics(legacy_rows, legacy_details, "Denay et al. (2023)")
        stats["label"] = "Denay et al. (2023)"
        stats["tag"] = "denay"
        all_stats.append(stats)

    # Aggregate summary across all datasets
    if len(all_stats) > 1:
        print(f"\n{'=' * 60}")
        print("  Multi-Study Aggregate Summary")
        print(f"{'=' * 60}")

        total_samples = sum(s["n_samples"] for s in all_stats)
        total_reads = sum(s["total_reads"] for s in all_stats)
        total_classified = sum(s["classified_reads"] for s in all_stats)
        total_pct = total_classified / total_reads * 100 if total_reads > 0 else 0

        total_tp = sum(s["spike_tp"] for s in all_stats)
        total_fp = sum(s["spike_fp"] for s in all_stats)
        total_fn = sum(s["spike_fn"] for s in all_stats)

        print(f"\n  {'Dataset':<40} {'Samples':>8} {'Reads':>12} {'Classified':>12} {'Rate':>8}")
        print(f"  {'-'*40} {'-'*8} {'-'*12} {'-'*12} {'-'*8}")
        for s in all_stats:
            print(f"  {s['label']:<40} {s['n_samples']:>8} {s['total_reads']:>12,} "
                  f"{s['classified_reads']:>12,} {s['classified_pct']:>7.1f}%")
        print(f"  {'TOTAL':<40} {total_samples:>8} {total_reads:>12,} "
              f"{total_classified:>12,} {total_pct:>7.1f}%")

        if (total_tp + total_fn) > 0:
            agg_sens = total_tp / (total_tp + total_fn)
            agg_prec = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else float("nan")
            agg_f1 = (2 * agg_prec * agg_sens / (agg_prec + agg_sens)
                      if (agg_prec + agg_sens) > 0 else float("nan"))
            print(f"\n  Aggregate spike detection:")
            print(f"    TP={total_tp}, FP={total_fp}, FN={total_fn}")
            print(f"    Sensitivity: {agg_sens:.3f}")
            print(f"    Precision:   {agg_prec:.3f}")
            print(f"    F1:          {agg_f1:.3f}")

        # Per-category summary table across datasets
        print(f"\n  Per-category summary (for paper Table 6):")
        print(f"  {'Dataset':<20} {'Category':<30} {'N':>4} {'Class. rate':>12}")
        print(f"  {'-'*20} {'-'*30} {'-'*4} {'-'*12}")
        for s in all_stats:
            tag = s["tag"]
            for cat, info in sorted(s.get("cat_summary", {}).items()):
                rate_str = f"{info['class_min']:.0f}--{info['class_max']:.0f}%"
                if info["class_min"] == info["class_max"]:
                    rate_str = f"{info['class_mean']:.0f}%"
                print(f"  {tag:<20} {cat:<30} {info['n']:>4} {rate_str:>12}")

    elif len(all_stats) == 1:
        # Single dataset — print category summary
        s = all_stats[0]
        if s.get("cat_summary"):
            print(f"\n  Per-category summary:")
            print(f"  {'Category':<30} {'N':>4} {'Class. rate':>12}")
            print(f"  {'-'*30} {'-'*4} {'-'*12}")
            for cat, info in sorted(s["cat_summary"].items()):
                rate_str = f"{info['class_min']:.0f}--{info['class_max']:.0f}%"
                if info["class_min"] == info["class_max"]:
                    rate_str = f"{info['class_mean']:.0f}%"
                print(f"  {cat:<30} {info['n']:>4} {rate_str:>12}")


if __name__ == "__main__":
    main()
