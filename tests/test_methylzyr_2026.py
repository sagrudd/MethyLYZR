import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "MethyLYZR.py"
SPEC = importlib.util.spec_from_file_location("methylzyr", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class MethyLYZR2026Test(unittest.TestCase):
    def test_get_adaptive_params_interpolates_alpha_and_basecount(self):
        alpha, basecount = MODULE.get_adaptive_params(4000)

        self.assertEqual(alpha, 7.5)
        self.assertAlmostEqual(basecount, 988.2)

    def test_calibration_preserves_top_class_and_normalizes_posteriors(self):
        posteriors = pd.Series([0.9, 0.09, 0.01], index=["class_a", "class_b", "class_c"])

        calibrated, temperature = MODULE.calibrate_posteriors(posteriors, n_cpg=500)

        self.assertGreater(temperature, MODULE.CALIB_T_BASE_DEFAULT)
        self.assertEqual(calibrated.idxmax(), "class_a")
        self.assertAlmostEqual(calibrated.sum(), 1.0)
        self.assertLess(calibrated.iloc[0], posteriors.iloc[0])

    def test_predict_from_fingerprint_defaults_match_original_parameters(self):
        newX = pd.Series([1, 0], index=["cg1", "cg2"])
        feature_ids = pd.Series(["cg1", "cg2"])
        centroids = pd.DataFrame(
            {
                "class_a": [0.8, 0.2],
                "class_b": [0.3, 0.7],
            },
            index=["cg1", "cg2"],
        )
        weights = pd.DataFrame(
            {
                "class_a": [0.2, 0.3],
                "class_b": [0.4, 0.1],
            },
            index=["cg1", "cg2"],
        )
        prior = pd.Series([0.5, 0.5], index=["class_a", "class_b"])
        noise = np.array([0.05, 0.05])
        read_weights = np.array([1.0, 1.0])

        default_result = MODULE.predict_from_fingerprint(
            newX, feature_ids, centroids, weights, noise, prior, read_weights
        )
        explicit_original_result = MODULE.predict_from_fingerprint(
            newX, feature_ids, centroids, weights, noise, prior, read_weights,
            relief_alpha=-1.0, basecount=MODULE.BASECOUNT
        )

        np.testing.assert_allclose(
            default_result["posterior"].to_numpy(),
            explicit_original_result["posterior"].to_numpy(),
        )


if __name__ == "__main__":
    unittest.main()
