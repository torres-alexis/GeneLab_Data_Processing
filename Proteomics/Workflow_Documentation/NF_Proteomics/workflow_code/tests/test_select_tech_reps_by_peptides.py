"""Offline synthetic regressions; run with Python's standard-library unittest."""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
SPEC = importlib.util.spec_from_file_location("pick_best", BIN / "select_tech_reps_by_peptides.py")
pick = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pick)
METRIC = "distinct_identified_peptide_sequences"


def lfq(name, source="mouse1", **extra):
    return {"Sample Name": name, "Source Name": source, "Has Tech Reps": "TRUE",
            "data_file": name + ".mzML", "Factor Value[condition]": "flight", **extra}


def tmt(run, **extra):
    return {"run": run, "plex": "plex1", "TechRepMixture": "1", "fraction": "1",
            "Has Tech Reps": "TRUE", "data_file": run + ".mzML", **extra}


def select(rows, counts, mode="lfq"):
    cols = list(dict.fromkeys(k for row in rows for k in row))
    return pick.select_rows(rows, cols, counts, mode, METRIC)


class SelectionTests(unittest.TestCase):
    def test_best_not_first_and_order_preserved(self):
        rows = [lfq("a"), lfq("independent", "mouse2"), lfq("b")]
        kept, log = select(rows, {"a": 10, "b": 30})
        self.assertEqual(kept, [rows[1], rows[2]])
        self.assertEqual([r["row_number"] for r in log], [1, 2, 3])
        self.assertEqual(log[0]["kept_row_number"], 3)
        self.assertEqual(log[2]["quant_sample"], "b")
        self.assertEqual(log[0]["reason"], "most_peptides")

    def test_tie_uses_sheet_order(self):
        rows = [lfq("b"), lfq("a")]
        kept, log = select(rows, {"a": 3, "b": 3})
        self.assertEqual(kept, [rows[0]])
        self.assertEqual(log[0]["reason"], "tie_keep_first")

    def test_zero_is_not_missing(self):
        kept, log = select([lfq("a"), lfq("b")], {"a": 0, "b": 0})
        self.assertEqual(kept[0]["Sample Name"], "a")
        self.assertEqual(log[0]["reason"], "all_zero_keep_first")

    def test_missing_candidate_fails(self):
        for counts in ({"a": 10}, {}, {"a": 0}):
            with self.subTest(counts=counts), self.assertRaisesRegex(ValueError, "unambiguous"):
                select([lfq("a"), lfq("b")], counts)

    def test_nontechnical_rows_never_collapsed(self):
        for value in ("", "FALSE", "NO"):
            rows = [lfq("same", **{"Has Tech Reps": value})] * 2
            self.assertEqual(select(rows, {})[0], rows)
        rows = [{"run": "same"}, {"run": "same"}]
        self.assertEqual(select(rows, {})[0], rows)

    def test_invalid_boolean_fails(self):
        with self.assertRaisesRegex(ValueError, "invalid Has Tech Reps"):
            select([lfq("a", **{"Has Tech Reps": "perhaps"})], {})

    def test_flagged_and_unflagged_never_group(self):
        rows = [lfq("a"), lfq("b"), lfq("c", **{"Has Tech Reps": "FALSE"})]
        self.assertEqual(select(rows, {"a": 1, "b": 2})[0], rows[1:])

    def test_conditions_sources_and_fractions_separate(self):
        rows = [lfq("a", fraction="1"), lfq("b", fraction="1"), lfq("c", fraction="2"),
                lfq("d", **{"Factor Value[condition]": "ground"}), lfq("e", "mouse2")]
        self.assertEqual(select(rows, {"a": 2, "b": 3})[0], rows[1:])

    def test_lfq_fraction_guard_reused(self):
        rows = [lfq("same", data_file="a.mzML"), lfq("same", data_file="b.mzML")]
        with self.assertRaisesRegex(SystemExit, "requires a non-empty fraction"):
            select(rows, {"a": 1, "b": 2})

    def test_missing_source_rejected(self):
        with self.assertRaisesRegex(SystemExit, "requires Source Name"):
            select([lfq("a", source=""), lfq("b", source="")], {"a": 1, "b": 2})

    def test_tmt_groups(self):
        rows = [tmt("a"), tmt("b"), tmt("c", fraction="2"), tmt("d", plex="plex2"),
                tmt("e", TechRepMixture="2")]
        self.assertEqual(select(rows, {"a": 1, "b": 2}, "tmt")[0], rows[1:])

    def test_missing_tmt_keys_rejected(self):
        for key in ("plex", "TechRepMixture", "fraction"):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "TMT technical runs require"):
                select([tmt("a", **{key: ""})], {}, "tmt")

    def test_shared_aggregated_alias_rejected(self):
        rows = [lfq("pooled", data_file="a.mzML", fraction="1"),
                lfq("pooled", data_file="b.mzML", fraction="1")]
        with self.assertRaisesRegex(ValueError, "unambiguous"):
            select(rows, {"pooled": 20})
        self.assertEqual(select(rows, {"pooled": 20, "a": 1, "b": 2})[0], [rows[1]])

    def test_multiple_distinct_matching_aliases_rejected(self):
        rows = [lfq("a", run="other"), lfq("b")]
        with self.assertRaisesRegex(ValueError, "unambiguous"):
            select(rows, {"a": 1, "other": 1, "b": 2})

    def test_original_sample_not_used(self):
        rows = [lfq("a", **{"Original Sample Name": "old"}), lfq("b")]
        with self.assertRaisesRegex(ValueError, "unambiguous"):
            select(rows, {"old": 1, "b": 2})

    def test_filename_case_and_dots(self):
        rows = [lfq("a", data_file="/some/path/first.run.MzML"),
                lfq("b", data_file="/some/path/second.run.RAW")]
        self.assertEqual(select(rows, {"first.run": 1, "second.run": 2})[0], [rows[1]])


class TableTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def table(self, name, fields, rows, delimiter="\t"):
        path = self.root / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter=delimiter)
            writer.writerow(fields)
            writer.writerows(rows)
        return path

    def combined(self, rows, fields=None):
        return self.table("combined.tsv", fields or ["Peptide Sequence", "a Spectral Count", "b Spectral Count"], rows)

    def psm(self, name, rows):
        return self.table(name, ["Spectrum", "Peptide", "Is Decoy", "Is Contaminant"], rows)

    def test_counts_distinct_sequences_not_psms_or_duplicate_rows(self):
        path = self.combined([["AA", 100, 0], ["AA", 3, 0], ["BB", 0, 1], ["CC", 0, 2]])
        self.assertEqual(pick.per_sample_peptide_counts(path), {"a": 1, "b": 2})

    def test_protein_table_rejected(self):
        path = self.combined([["protein", 20]], ["Protein", "a Spectral Count"])
        with self.assertRaisesRegex(ValueError, "protein tables"):
            pick.per_sample_peptide_counts(path)

    def test_spectral_counts_preferred_to_mbr_intensities(self):
        path = self.combined([["AA", 1, 0, 0, 100]], ["Peptide Sequence", "a Spectral Count", "b Spectral Count", "a Intensity", "b Intensity"])
        self.assertEqual(pick.per_sample_peptide_counts(path, True), {"a": 1, "b": 0})

    def test_intensity_requires_opt_in_and_ignores_maxlfq(self):
        path = self.combined([["AA", 1, 2, 3]], ["Peptide Sequence", "a Intensity", "a MaxLFQ Intensity", "Combined Intensity"])
        with self.assertRaisesRegex(ValueError, "allow-intensity-fallback"):
            pick.per_sample_peptide_counts(path)
        self.assertEqual(pick.per_sample_peptide_counts(path, True), {"a": 1})

    def test_sample_names_containing_metric_words_preserved(self):
        path = self.combined([["AA", 1, 1, 1]], ["Peptide Sequence", "UniqueMouse Spectral Count", "TotalMouse Spectral Count", "CombinedMouse Spectral Count"])
        self.assertEqual(len(pick.per_sample_peptide_counts(path)), 3)

    def test_bad_counts_fail(self):
        for value in ("oops", "inf", "-1"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "invalid measurement"):
                pick.per_sample_peptide_counts(self.combined([["AA", value, 0]]))

    def test_duplicate_headers_fail(self):
        path = self.combined([["AA", 1, 2]], ["Peptide Sequence", "a Spectral Count", "a Spectral Count"])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            pick.per_sample_peptide_counts(path)

    def test_blank_peptide_fails(self):
        with self.assertRaisesRegex(ValueError, "missing Peptide Sequence"):
            pick.per_sample_peptide_counts(self.combined([["", 1, 2]]))

    def test_psm_deduplicates_files_scans_charges_and_modifications(self):
        path = self.psm("psm.tsv", [["a.run.001.001.2", "AA", "false", "false"],
                                     ["a.run.002.002.3", "AA", "false", "false"],
                                     ["a.run.003.003.2", "BB", "false", "false"]])
        self.assertEqual(pick.per_run_psm_counts([path, path]), {"a.run": 2})

    def test_psm_decoys_contaminants_and_zero(self):
        path = self.psm("psm.tsv", [["a.1.1.2", "AA", "true", "false"],
                                     ["b.2.2.2", "BB", "false", "true"],
                                     ["c.3.3.2", "CC", "false", "false"]])
        self.assertEqual(pick.per_run_psm_counts([path]), {"a": 0, "b": 0, "c": 1})

    def test_psm_invalid_spectrum_or_flags_fail(self):
        for row in (["broken", "AA", "false", "false"], ["a.1.1.2", "AA", "maybe", "false"]):
            with self.subTest(row=row), self.assertRaises(ValueError):
                pick.per_run_psm_counts([self.psm("psm.tsv", [row])])

    def test_psm_required_columns(self):
        path = self.table("psm.tsv", ["Spectrum", "Peptide"], [["a.1.1.2", "AA"]])
        with self.assertRaisesRegex(ValueError, "filtered PSM column"):
            pick.per_run_psm_counts([path])

    def cli(self, extra, quant=None, rows=None):
        rows = rows or [lfq("a"), lfq("b")]
        self.sheet = self.table("sheet.csv", rows[0].keys(), [r.values() for r in rows], ",")
        quant = quant or self.combined([["AA", 1, 0], ["BB", 0, 1], ["CC", 0, 1]])
        return subprocess.run([sys.executable, str(BIN / "select_tech_reps_by_peptides.py"),
                               "--mode", "lfq", "--input", str(self.sheet), "--quant", str(quant),
                               *extra], text=True, capture_output=True)

    def test_combined_explicit_junk_excluded(self):
        path = self.combined([["AA", "contam_sp|X|JUNK", 2, 0],
                              ["BB", "rev_sp|Y|DECOY", 1, 0],
                              ["CC", "sp|Z|TARGET", 0, 1]],
                             ["Peptide Sequence", "Protein", "a Spectral Count", "b Spectral Count"])
        self.assertEqual(pick.per_sample_peptide_counts(path), {"a": 0, "b": 1})

    def test_malformed_quantitative_rows_fail(self):
        path = self.combined([["AA", 1, 2, "extra"]])
        with self.assertRaisesRegex(ValueError, "wrong number of fields"):
            pick.per_sample_peptide_counts(path)
        path = self.psm("psm.tsv", [["a.1.1.2", "AA", "false", "false", "extra"]])
        with self.assertRaisesRegex(ValueError, "wrong number of fields"):
            pick.per_run_psm_counts([path])
        path = self.psm("psm.tsv", [["a.1.1.2", "AA", "false"]])
        with self.assertRaisesRegex(ValueError, "wrong number of fields"):
            pick.per_run_psm_counts([path])

    def test_psm_duplicate_unrelated_headers_supported(self):
        path = self.table("psm.tsv", ["Spectrum", "Peptide", "Is Decoy", "Is Contaminant", "Intensity NA", "Intensity NA"],
                          [["a.1.1.2", "AA", "false", "false", 1, 2]])
        self.assertEqual(pick.per_run_psm_counts([path]), {"a": 1})

    def test_cli_psm_and_provenance(self):
        rows = [tmt("a"), tmt("b")]
        sheet = self.table("sheet.csv", rows[0].keys(), [r.values() for r in rows], ",")
        psm = self.psm("psm.tsv", [["a.1.1.2", "AA", "false", "false"],
                                   ["b.2.2.2", "BB", "false", "false"],
                                   ["b.3.3.2", "CC", "false", "false"]])
        out, log = self.root / "out.csv", self.root / "log.tsv"
        result = subprocess.run([sys.executable, str(BIN / "select_tech_reps_by_peptides.py"),
                                 "--mode", "tmt", "--input", str(sheet), "--psm", str(psm),
                                 "--output", str(out), "--log", str(log)], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with out.open() as handle:
            self.assertEqual([r["run"] for r in csv.DictReader(handle)], ["b"])
        with log.open() as handle:
            entry = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(entry["evidence_type"], "filtered_psm")
        self.assertEqual(entry["evidence_files"], str(psm.resolve()))

    def test_cli_writes_auditable_output(self):
        out, log = self.root / "selected.csv", self.root / "logs" / "selected.tsv"
        result = self.cli(["--output", str(out), "--log", str(log)])
        self.assertEqual(result.returncode, 0, result.stderr)
        with out.open() as handle:
            self.assertEqual([r["Sample Name"] for r in csv.DictReader(handle)], ["b"])
        with log.open() as handle:
            audit = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(audit[1]["peptide_count"], "2")
        self.assertEqual(audit[1]["selection_metric"], METRIC)

    def test_validation_failure_does_not_touch_outputs(self):
        out, log = self.root / "out.csv", self.root / "log.tsv"
        out.write_text("existing output")
        log.write_text("existing log")
        path = self.combined([["AA", 1]], ["Peptide Sequence", "a Spectral Count"])
        result = self.cli(["--output", str(out), "--log", str(log)], quant=path)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(out.read_text(), "existing output")
        self.assertEqual(log.read_text(), "existing log")

    def test_input_overwrite_rejected(self):
        result = self.cli(["--output", str(self.root / "sheet.csv"), "--log", str(self.root / "log.tsv")])
        self.assertEqual(result.returncode, 2)
        self.assertIn("distinct", result.stderr)
        self.assertTrue(self.sheet.read_text().startswith("Sample Name"))

    def test_pipeline_default_unchanged(self):
        config = (BIN.parent / "nextflow.config").read_text()
        self.assertRegex(config, r"tech_rep\s*=\s*'first'")


if __name__ == "__main__":
    unittest.main()
