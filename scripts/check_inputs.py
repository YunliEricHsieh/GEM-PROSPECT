"""Read-only, standard-library checks for documented workflow inputs.

These checks inspect file presence, table schemas, and supplied artifact joins.
They do not load a scientific model, execute optimization, or validate an
installation. Paths are resolved against this script's checkout, not cwd.
"""

import argparse
import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent
RATIO_COLUMNS = [
    "auto_mixo", "auto_hetero", "mixo_hetero", "auto_auto_CO2",
    "mixo_mixo_CO2", "mixo_hypo10_mixo", "mixo_hypo25_mixo", "mixo_hypo75_mixo",
]
MODEL_PATHS = [
    f"Data/pciCre1355/NDLadpraw_{mode}_Rep1.xml"
    for mode in ("Autotrophic", "Mixotrophic", "Heterotrophic")
]
CARBON_SOURCES = [
    "Glucose", "Mannitol", "Glucosamine", "Glycerol", "Maltose", "Gluconate",
    "Xylose", "Sorbitol", "Ribose", "Succinate", "Galactose", "Lactate",
    "Alanine", "Pyruvate", "Oxoglutarate", "Acetate",
]


def requirements(workflow):
    """Return (relative path, required headers, delimiter) input contracts."""
    if workflow == "preprocessing":
        return [
            ("Data/Mutant_phenotypes_table.xlsx", [], None),
            ("Data/mart_Cre_Uniprot.txt", ["Gene Name", "UniProt ID"], "\t"),
        ]
    if workflow == "ecoli":
        return [
            ("Data/Ecoli/iML1515.mat", [], None),
            ("Data/Ecoli/iML1515.xml", [], None),
            ("Data/Ecoli/Ecoli_gene_essentiality.csv", ["Gene_ID", "Essentiality_0_1"], ";"),
        ] + [(f"Results/FVA/Ecoli/{source}.csv", ["RxnID", "minFlux", "maxFlux"], ",")
             for source in CARBON_SOURCES]
    items = [(p, [], None) for p in MODEL_PATHS]
    items.append(("Data/Mutant_phenotypes_table_filtered_final.csv",
                  ["GeneID", "UniProtID"] + RATIO_COLUMNS, ","))
    items += [
        (f"Data/Reactions/list_of_{kind}_rxns.csv", ["RxnIndex", "RxnID"], ",")
        for kind in ("irreversible", "reversible")
    ]
    if workflow == "quick-start":
        return items
    items.append(("Data/GO_table_filtered.txt", ["UniProtID", "GO_Name", "GO_ID", "GO_Aspect"], "\t"))
    for mode in ("auto", "mixo", "hetero"):
        items.append((f"Results/FVA/{mode}_FVA_10p.csv", ["RxnID", "minFlux", "maxFlux"], ","))
        for suffix in ("", "_re"):
            items.append((f"Results/flux_sampling/{mode}_sampling{suffix}.csv",
                          ["RxnIndex", "RxnID", "meanFlux", "stdFlux"], ","))
    for name in ("without_GPR", "with_proteins_in_both", "with_at_least_one_protein_in_both"):
        headers = ["RxnIndex", "RxnID"] + ([] if name == "without_GPR" else ["Enzymes"])
        items.append((f"Data/Reactions/list_of_rxns_{name}.csv", headers, ","))
    for suffix in ("", "_Re"):
        items.append((f"Results/screens/Max_flux_screen_8{suffix}.csv", ["EnzymeID"], ","))
    return items


def read_rows(relative_path, delimiter=","):
    """Read a supplied text table without importing scientific packages."""
    with (ROOT / relative_path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def inspect(workflow):
    problems = []
    tables = {}
    contracts = requirements(workflow)
    for relative_path, headers, delimiter in contracts:
        path = ROOT / relative_path
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"Missing or empty: {relative_path}")
            continue
        if delimiter is not None:
            try:
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    actual_headers = next(csv.reader(handle, delimiter=delimiter), [])
                missing = sorted(set(headers) - set(actual_headers))
                if missing:
                    problems.append(f"Missing headers in {relative_path}: {', '.join(missing)}")
                else:
                    tables[relative_path] = actual_headers
            except (OSError, UnicodeError, csv.Error) as error:
                problems.append(f"Cannot read {relative_path}: {error}")

    phenotype_path = "Data/Mutant_phenotypes_table_filtered_final.csv"
    if phenotype_path in tables:
        rows = read_rows(phenotype_path)
        example_rows = [r for r in rows if r["UniProtID"] == "A0A2K3E4Q7"]
        if workflow == "quick-start" and len(example_rows) != 1:
            problems.append("The quick-start phenotype ID must have exactly one row.")
        if workflow == "quick-start" and len(example_rows) == 1:
            try:
                import math
                if not all(math.isfinite(float(example_rows[0][c])) for c in RATIO_COLUMNS):
                    problems.append("Quick-start phenotype ratios are not all finite.")
            except ValueError:
                problems.append("Quick-start phenotype ratios are not numeric.")
        print(f"Processed phenotype rows: {len(rows)}")

    # Compare column IDs with maps and sampling references, without thresholding
    # or altering results. This is useful when inspecting a fresh clone.
    for kind, suffix, sampling_suffix in (("irreversible", "", ""), ("reversible", "_Re", "_re")):
        map_path = f"Data/Reactions/list_of_{kind}_rxns.csv"
        result_path = f"Results/screens/Max_flux_screen_8{suffix}.csv"
        if map_path not in tables:
            continue
        mapping = {r["RxnIndex"]: r["RxnID"] for r in read_rows(map_path)}
        if workflow == "quick-start":
            expected = "3SPYRSPh" if kind == "irreversible" else "10FTHFtx_REV"
            if expected not in mapping.values():
                problems.append(f"Example reaction {expected} missing from {map_path}")
        if result_path not in tables:
            continue
        columns = set(tables[result_path]) - {"EnzymeID"}
        unmapped = sorted(columns - set(mapping))
        if unmapped:
            problems.append(f"Unmapped reaction columns in {result_path}: {unmapped[:5]}")
        for mode in ("auto", "mixo", "hetero"):
            sample_path = f"Results/flux_sampling/{mode}_sampling{sampling_suffix}.csv"
            if sample_path not in tables:
                continue
            reference = {r["RxnIndex"]: r["RxnID"] for r in read_rows(sample_path)}
            mismatched = sorted(c for c in columns if c not in reference or reference[c] != mapping.get(c))
            if mismatched:
                problems.append(f"Reference/index mismatch in {sample_path}: {mismatched[:5]}")

    for problem in problems:
        print(f"ERROR: {problem}", file=sys.stderr)
    if problems:
        print("See docs/input_formats.md and docs/reproducing_manuscript.md.", file=sys.stderr)
        return 1
    print(f"PASS: {workflow}: {len(contracts)} required input files checked.")
    print("This does not verify solver execution, scientific validity, or dependency versions.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow", choices=(
        "quick-start", "chlamydomonas", "ecoli", "preprocessing"),
        default="quick-start")
    return inspect(parser.parse_args().workflow)


if __name__ == "__main__":
    sys.exit(main())
