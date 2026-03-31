import importlib.util
import multiprocessing
import tempfile
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

        self.assertEqual(index["chr1"]["forward"][10], "cg0001")
        self.assertEqual(index["chr1"]["reverse"][21], "cg0002")
        self.assertEqual(index["chr2"]["forward"][30], "cg0003")

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

    def test_chunk_unique_methylation_alignments_streams_unique_reads(self):
        alignments = {
            "read1": [("chr1", [], {("C", 0, "m"): [(1, 255)]}, False, "st", "rg", "rn", "read1", 10, 100, 60)],
            "read2": [("chr1", [], {}, False, "st", "rg", "rn", "read2", 10, 100, 60)],
            "read3": [("chr1", [], {("C", 0, "m"): [(2, 255)]}, False, "st", "rg", "rn", "read3", 10, 100, 60)],
        }

        chunks = list(MODULE.chunk_unique_methylation_alignments(alignments, 1))

        self.assertEqual([[chunk[0][7]] for chunk in chunks], [["read1"], ["read3"]])

    def test_process_methylation_data_maps_modified_positions_without_pandas_joins(self):
        sites_index = {
            "chr1": {
                "forward": {101: "cg101"},
                "reverse": {205: "cg205"},
            }
        }
        read_row = (
            "chr1",
            [(1, 101), (2, 102), (3, 103)],
            {("C", 0, "m"): [(1, 255), (2, 10)]},
            False,
            "2026-03-31T19:30:00Z",
            "runA",
            "rn",
            "read1",
            12,
            1234,
            60,
        )

        rows = MODULE.process_methylation_data(read_row, sites_index)

        self.assertEqual(
            rows,
            [["cg101", 1.0, 1, 1, "read1", "2026-03-31T19:30:00Z", "runA", 12, 1234, 60]],
        )

    def test_write_chunk_file_writes_feather_with_worker_label(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [
                ["cg101", 1.0, 1, 1, "read1", "2026-03-31T19:30:00Z", "runA", 12, 1234, 60],
                ["cg102", 0.5, 1, 0, "read2", "2026-03-31T19:31:00Z", "runA", 13, 1235, 55],
            ]
            chunk_path = MODULE.write_chunk_file(tmpdir, "worker-00", 0, rows)
            self.assertTrue(Path(chunk_path).exists())
            self.assertEqual(Path(chunk_path).name, "worker-00-chunk-00000.feather")
            merged = pd.read_feather(chunk_path)
            self.assertEqual(list(merged["epic_id"]), ["cg101", "cg102"])

    def test_methylation_reader_flushes_pending_rows_when_sentinel_arrives(self):
        methylation_queue = multiprocessing.Queue()
        finished = multiprocessing.Value("i", 0, lock=True)
        analysed = multiprocessing.Value("i", 0, lock=True)
        sites_index = {"chr1": {"forward": {101: "cg101"}, "reverse": {}}}
        read_row = (
            "chr1",
            [(1, 101)],
            {("C", 0, "m"): [(1, 255)]},
            False,
            "2026-03-31T19:30:00Z",
            "runA",
            "rn",
            "read1",
            12,
            1234,
            60,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            methylation_queue.put([read_row])
            methylation_queue.put(None)

            MODULE.methylation_reader(
                methylation_queue,
                sites_index,
                finished,
                analysed,
                tmpdir,
                100,
                "worker-00",
            )

            chunk_paths = sorted(Path(tmpdir).glob("*.feather"))
            self.assertEqual([path.name for path in chunk_paths], ["worker-00-chunk-00000.feather"])
            written = pd.read_feather(chunk_paths[0])
            self.assertEqual(list(written["epic_id"]), ["cg101"])
            self.assertEqual(finished.value, 1)
            self.assertEqual(analysed.value, 1)


if __name__ == "__main__":
    unittest.main()
