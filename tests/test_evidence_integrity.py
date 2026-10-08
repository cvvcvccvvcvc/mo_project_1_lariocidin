"""Regression cases for the scientific reporting and provenance bugs found in review."""
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from assess_free_bgc_trees import validate_taxa
from audit_maize_provenance import validate_accessions
from summarize_hgt_trees import apply_holm, au_constrained_row, split_clade, support_values


class EvidenceIntegrity(unittest.TestCase):
    def test_significant_result_changes_verdict_and_holm_controls_family(self):
        rows = [dict(role=role, constrained_au_p=p) for role, p in
                zip(("B1", "B2", "C", "D"), (.001, .02, None, .4))]
        apply_holm(rows)
        self.assertEqual(rows[0]["h0_decision"], "constraint_rejected")
        self.assertAlmostEqual(rows[0]["holm_adjusted_p"], .004)
        self.assertEqual(rows[1]["h0_decision"], "constraint_not_rejected")
        self.assertAlmostEqual(rows[1]["holm_adjusted_p"], .06)
        self.assertIsNone(rows[2]["holm_adjusted_p"])
        self.assertEqual(rows[2]["h0_decision"], "free_tree_compatible")

    def test_missing_gene_cannot_silently_reduce_test_family(self):
        with self.assertRaises(ValueError):
            apply_holm([dict(role="B1", constrained_au_p=.001)])

    def test_au_parser_uses_named_table(self):
        report = "2 0 0 0 0 0 0 0 0 0 0\nUSER TREES\n----------\n"
        report += "Tree logL deltaL bp-RELL p-KH p-SH c-ELW p-AU\n-----\n"
        report += "1 -100 0 .9 + .9 + 1 + .9 + .99 +\n"
        report += "2 -105 5 .1 - .1 - .1 - .1 - .003 -\n\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "au.iqtree"
            path.write_text(report)
            self.assertEqual(au_constrained_row(path), (-105., 5., .003))
            path.write_text(report.replace("p-AU", "not-AU"))
            with self.assertRaises(ValueError):
                au_constrained_row(path)

    def test_missing_or_duplicate_tree_taxa_rejected(self):
        for observed in (("A", "B"), ("A", "B", "B", "C")):
            with self.assertRaises(ValueError):
                validate_taxa(observed, {"A", "B", "C"}, "example")

    def test_split_is_unrooted_and_absence_has_no_support(self):
        tree = Phylo.read(StringIO("(A,(B,(C,D)90/99));"), "newick")
        self.assertIsNotNone(split_clade(tree, {"A", "B"}))
        self.assertIsNone(split_clade(tree, {"A", "C"}))
        self.assertEqual(support_values(None), (None, None))

    def test_biosample_identity_is_checked_not_just_xml_format(self):
        def samples(*ids):
            return [ET.Element("BioSample", accession=identifier) for identifier in ids]
        validate_accessions(samples("SAMN2", "SAMN1"), ["SAMN1", "SAMN2"])
        for received in (samples("SAMN3"), samples("SAMN1"), samples("SAMN1", "SAMN1", "SAMN2")):
            with self.assertRaises(ValueError):
                validate_accessions(received, ["SAMN1", "SAMN2"])


if __name__ == "__main__":
    unittest.main()
