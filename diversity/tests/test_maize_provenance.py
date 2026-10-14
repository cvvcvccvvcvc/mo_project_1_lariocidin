"""Keep BioSample identity checks for the retained Maize baseline."""
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/steps"))
from audit_maize_provenance import validate_accessions


class MaizeProvenance(unittest.TestCase):
    def test_biosample_identity_is_checked_not_just_xml_format(self):
        def samples(*ids):
            return [ET.Element("BioSample", accession=identifier) for identifier in ids]
        validate_accessions(samples("SAMN2", "SAMN1"), ["SAMN1", "SAMN2"])
        for received in (samples("SAMN3"), samples("SAMN1"), samples("SAMN1", "SAMN1", "SAMN2")):
            with self.assertRaises(ValueError):
                validate_accessions(received, ["SAMN1", "SAMN2"])
