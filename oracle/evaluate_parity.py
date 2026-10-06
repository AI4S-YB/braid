#!/usr/bin/env python3
"""Differential CLI audit; writes inputs, commands, outputs and JSON evidence.

Build plena first. Pass a local copy of upstream flair_combine.py and its
bed_to_gtf.py from the same revision. Their hashes are recorded, not trusted
as proof of an upstream revision. No production implementation is changed.

Example:
  python3 oracle/evaluate_parity.py --flair-source /tmp/flair_combine.py \
    --gtf-source /tmp/bed_to_gtf.py --out /tmp/plena-parity-audit

Unused upstream imports are stubbed and fail if called; combine and GTF
conversion execute the supplied source without rewriting it. For successful
runs, compare the complete output file set and bytes. For failed runs,
record exit codes, stderr and partial files without claiming output parity.
"""

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import runpy
import subprocess
import sys
import types


def oracle_main(source, gtf_source, args):
    def unused(*args, **kwargs):
        raise RuntimeError("Audit stub called: install the real upstream dependency")

    for name in ("pipettor", "pysam", "flair.flair_transcriptome"):
        module = types.ModuleType(name)
        module.__getattr__ = lambda name: unused
        sys.modules[name] = module
    flair = types.ModuleType("flair")
    flair.FlairInputDataError = type("FlairInputDataError", (Exception,), {})
    sys.modules["flair"] = flair
    spec = importlib.util.spec_from_file_location("flair.bed_to_gtf", gtf_source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    sys.argv = [str(source), *args]
    runpy.run_path(str(source), run_name="__main__")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bed_row(name, start, end, intron=None, strand="+", chrom="chr1"):
    if intron is None:
        sizes, rels = [end - start], [0]
    else:
        sizes, rels = [intron[0] - start, end - intron[1]], [0, intron[1] - start]
    return "\t".join(map(str, [chrom, start, end, name, 40, strand, start, end, 0,
                                len(sizes), ",".join(map(str, sizes)) + ",",
                                ",".join(map(str, rels)) + ","])) + "\n"


def audit(args):
    repo = Path(__file__).resolve().parents[1]
    binary = (repo / "target/debug/plena").resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = Path(args.flair_source).resolve()
    gtf_source = Path(args.gtf_source).resolve()
    oracle = [sys.executable, str(Path(__file__).resolve()), "--oracle",
              str(source), str(gtf_source)]
    results = []

    def execute(cmd, directory):
        directory.mkdir()
        p = subprocess.run(cmd, capture_output=True, env={**os.environ, "PYTHONHASHSEED": "0"})
        (directory / "stdout.txt").write_bytes(p.stdout)
        (directory / "stderr.txt").write_bytes(p.stderr)
        (directory / "command.json").write_text(json.dumps(cmd, indent=2))
        return p.returncode

    def compare(name, rows, flags=(), golden=None):
        case = out / name
        case.mkdir(exist_ok=True)
        manifest = case / "manifest.tsv"
        manifest.write_text("\n".join(rows) + "\n")
        codes = {}
        for label, cmd in (("upstream", oracle), ("plena", [str(binary), "flair", "combine"])):
            directory = case / label
            codes[label] = execute([*cmd, "-m", str(manifest), "-o", str(directory / "combined"), *flags], directory)
        files = {label: {p.name: p.read_bytes() for p in (case / label).glob("combined*")}
                 for label in codes}
        differences = [name for name in sorted(files["upstream"].keys() | files["plena"].keys())
                       if files["upstream"].get(name) != files["plena"].get(name)]
        result = {"case": name, "exit_codes": codes, "different_files": differences,
                  "status": "match" if codes == {"upstream": 0, "plena": 0} and not differences
                  else "both_failed" if all(codes.values()) else "mismatch"}
        if golden:
            expected = {p.name: p.read_bytes() for p in golden.glob("combined*")}
            result["oracle_matches_frozen_golden"] = files["upstream"] == expected
        results.append(result)
        return case

    fixtures = repo / "tests/parity/flair_combine"
    official = fixtures / "official/in"
    for name, n in (("bed-only", 1), ("bed-fa", 2), ("bed-fa-map", 3)):
        paths = [official / f"collapse.{suffix}" for suffix in
                 ("isoforms.bed", "isoforms.fa", "isoform.read.map.txt")]
        compare("official-" + name, ["\t".join(["A1", "isoforms", *map(str, paths[:n])])],
                golden=fixtures / "official" / name)
    synthetic = fixtures / "synthetic/in"
    for name, samples, flags in (
        ("usageandlongest", [("S1", "s1"), ("S2", "s2")], ["-c"]),
        ("usageonly", [("S1", "s1"), ("S2", "s2")], ["-f", "usageonly"]),
        ("filternone", [("S1", "s1"), ("S2", "s2")], ["-f", "none"]),
        ("dropped", [("D", "drop")], []),
        ("dropped_none", [("D", "drop")], ["-f", "none"]),
        ("numeric", [("N", "num")], ["-f", "3"]),
        ("se_off", [("E", "se")], []),
        ("se_on", [("E", "se")], ["-s"]),
        ("fusion", [("F", "fus")], ["-c"]),
    ):
        rows = ["\t".join([label, "fusionisoform" if stem == "fus" else "isoforms",
                           *[str(synthetic / (stem + ext)) for ext in (".bed", ".fa", ".map")]])
                for label, stem in samples]
        compare("synthetic-" + name, rows, flags, fixtures / "synthetic" / name)

    case = out / "consecutive-filtered-chains"
    case.mkdir()
    (case / "in.bed").write_text("".join(bed_row(f"t{i}_ENSG{i}", i*1000+100, i*1000+400,
                                                    (i*1000+200, i*1000+300)) for i in range(3)))
    (case / "in.map").write_text("t0_ENSG0\ta\nt1_ENSG1\tb\nt2_ENSG2\tc\n" +
                                  "\n".join(f"bulk_ENSG{i}\t" + ",".join(["r"]*99) for i in (1, 2)) + "\n")
    compare(case.name, [f"S\tisoforms\t{case}/in.bed\t\t{case}/in.map"])

    case = out / "mixed-fasta-missing-unused-file"
    case.mkdir()
    (case / "in.bed").write_text(bed_row("t_ENSG1", 100, 400, (200, 300)))
    compare(case.name, [f"A\tisoforms\t{case}/in.bed\t{case}/absent.fa", f"B\tisoforms\t{case}/in.bed"])

    case = out / "fasta-blank-line"
    case.mkdir()
    (case / "in.bed").write_text(bed_row("t_ENSG1", 100, 400, (200, 300)))
    (case / "in.fa").write_text(">t_ENSG1\nACGT\n\nACGT\n")
    compare(case.name, [f"A\tisoforms\t{case}/in.bed\t{case}/in.fa"])

    rng = random.Random(20261004)
    for seed in range(args.random_cases):
        case = out / f"random-{seed:03d}"
        case.mkdir()
        rows = []
        flags = ["-w", str(rng.choice([0, 10, 50, 200, 201])), "-p", str(rng.choice([0, 10, 25, 50])),
                 "-f", rng.choice(["usageandlongest", "usageonly", "none", "3", "10"])]
        if rng.choice([True, False]):
            flags.append("-s")
        if rng.choice([True, False]):
            flags.append("-c")
        for sample in range(3):
            bed, fa, maps = [], [], []
            for i in range(20):
                chain = rng.randrange(5)
                start = chain*2000 + rng.choice([100, 110, 111, 150, 300])
                end = chain*2000 + rng.choice([900, 1000, 1100])
                intron = None if chain == 4 else (chain*2000+400, chain*2000+700)
                name = f"t{sample}-{i}_ENSG{chain}"
                strand = "+" if chain % 2 else "-"
                bed.append(bed_row(name, start, end, intron, strand))
                fa.append(f">{name}\n" + "A"*(end-start-(300 if intron else 0)) + "\n")
                maps.append(name + "\t" + ",".join(f"r{j}" for j in range(rng.choice([1, 2, 10, 50]))) + "\n")
            stem = case / f"s{sample}"
            for ext, text in (("bed", bed), ("fa", fa), ("map", maps)):
                stem.with_suffix("." + ext).write_text("".join(text))
            rows.append(f"S{sample}\tisoforms\t{stem}.bed\t{stem}.fa\t{stem}.map")
        compare(case.name, rows, flags)

    golden_failures = [r["case"] for r in results if r.get("oracle_matches_frozen_golden") is False]
    summary = {"source_sha256": {str(source): sha(source), str(gtf_source): sha(gtf_source)},
               "binary_sha256": sha(binary), "random_seed": 20261004,
               "oracle_golden_failures": golden_failures,
               "totals": {s: sum(r["status"] == s for r in results)
                          for s in ("match", "mismatch", "both_failed")}, "cases": results}
    (out / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["totals"]))
    for result in results:
        if result["status"] == "mismatch":
            print(json.dumps(result))
    print("Evidence:", out / "results.json")
    return int(bool(golden_failures) or any(r["status"] == "mismatch" for r in results))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--oracle":
        oracle_main(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--flair-source", required=True)
        parser.add_argument("--gtf-source", required=True)
        parser.add_argument("--out", required=True)
        parser.add_argument("--random-cases", type=int, default=100)
        sys.exit(audit(parser.parse_args()))
