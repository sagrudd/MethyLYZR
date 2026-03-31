import importlib.util
import unittest
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "bam2feather.py"
SPEC = importlib.util.spec_from_file_location("bam2feather", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class Bam2FeatherHelpersTest(unittest.TestCase):
    def test_build_sites_index_creates_forward_and_reverse_lookups(self):
        sites = pd.DataFrame(
            [
                ("chr1", 10, 11, "cg0001"),
                ("chr1", 20, 21, "cg0002"),
                ("chr2", 30, 31, "cg0003"),
            ],
            columns=["chromosome", "start", "end", "epic_id"],
        )

        index = MODULE.build_sites_index(sites)

        self.assertEqual(index["chr1"]["forward"].loc[10, "epic_id"], "cg0001")
        self.assertEqual(index["chr1"]["reverse"].loc[21, "epic_id"], "cg0002")
        self.assertEqual(index["chr2"]["forward"].loc[30, "epic_id"], "cg0003")

    def test_iter_unique_methylation_alignments_filters_non_unique_and_empty(self):
        alignments = {
            "read1": [("chr1", [], {("C", 0, "m"): [(1, 255)]}, False, "st", "rg", "rn", "read1", 10, 100, 60)],
            "read2": [
                ("chr1", [], {("C", 0, "m"): [(2, 255)]}, False, "st", "rg", "rn", "read2", 10, 100, 60),
                ("chr1", [], {("C", 0, "m"): [(3, 255)]}, False, "st", "rg", "rn", "read2", 10, 100, 60),
            ],
            "read3": [("chr1", [], {}, False, "st", "rg", "rn", "read3", 10, 100, 60)],
        }

        result = list(MODULE.iter_unique_methylation_alignments(alignments))

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][7], "read1")

    def test_chunk_records_splits_into_fixed_sizes(self):
        chunks = list(MODULE.chunk_records(list(range(7)), 3))
        self.assertEqual(chunks, [[0, 1, 2], [3, 4, 5], [6]])


if __name__ == "__main__":
    unittest.main()
