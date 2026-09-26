#!/usr/bin/env python3
"""Checks that incomplete or corrupted historical evidence fails closed."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
import analyze


class EvidenceRejectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="historical-check-")
        self.root = Path(self.temp.name)
        self.original_root = analyze.ROOT
        self.bundle = self.root / "bundle"
        shutil.copytree(self.original_root, self.bundle)
        analyze.ROOT = self.bundle

    def tearDown(self):
        analyze.ROOT = self.original_root
        self.temp.cleanup()

    def test_corrupted_compressed_record(self):
        path = next((self.bundle / "raw/july13").glob("*.time.gz"))
        path.write_bytes(path.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "compressed input hash/size mismatch"):
            analyze.unpack(self.root / "restore")

    def test_duplicate_record(self):
        path = self.bundle / "inputs.json"
        doc = json.loads(path.read_text())
        doc["files"].append(doc["files"][0])
        path.write_text(json.dumps(doc))
        with self.assertRaisesRegex(ValueError, "duplicate bundle member"):
            analyze.unpack(self.root / "restore")

    def test_nonzero_exit_status(self):
        restored = self.root / "restore"
        analyze.unpack(restored)
        path = restored / "july23/clean_j16.time"
        path.write_text(path.read_text().replace("Exit status: 0", "Exit status: 1"))
        with self.assertRaises(ValueError):
            analyze.jul23(restored / "july23", analyze.module("analyze_dystrophin_bench"))

    def test_missing_repeat(self):
        restored = self.root / "restore"
        analyze.unpack(restored)
        (restored / "july13/rep3.time").unlink()
        dyst = analyze.module("analyze_dystrophin_bench")
        with self.assertRaises(ValueError):
            dyst.collect_observations(restored / "july13", dyst.GROUPS)


if __name__ == "__main__":
    unittest.main()
