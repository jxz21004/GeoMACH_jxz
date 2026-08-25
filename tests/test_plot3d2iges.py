import tempfile
import unittest
from pathlib import Path

import numpy as np

from GeoMACH.utilities import P3dSurfClass

try:
    from GeoMACH.utilities.plot3d2iges import fit_plot3d_to_iges
except ImportError:
    fit_plot3d_to_iges = None


def _write_surface_plot3d(path, surface):
    ni, nj = surface.shape[:2]
    vals = ["1", str(ni), str(nj), "1"]
    for idim in range(3):
        vals.extend(f"{x:.16e}" for x in surface[..., idim].reshape(-1, order="F"))
    path.write_text("\n".join(vals) + "\n")


def _make_surface(ni=5, nj=4):
    surface = np.empty((ni, nj, 3), dtype=float, order="F")
    for j in range(nj):
        for i in range(ni):
            u = i / (ni - 1)
            v = j / (nj - 1)
            surface[i, j] = [u, v, 0.1 * u * (1.0 - u) + 0.05 * v]
    return surface


@unittest.skipIf(fit_plot3d_to_iges is None, "GeoMACH compiled extensions are not built")
class TestPlot3D2IGES(unittest.TestCase):
    def test_fit_and_export_default_meters(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            p3d = tmp / "surface.xyz"
            iges = tmp / "surface.igs"
            surface = _make_surface()
            _write_surface_plot3d(p3d, surface)

            original = P3dSurfClass(p3d).SurfGen()[0]
            bse = fit_plot3d_to_iges(p3d, iges)
            fitted = bse.vec["pt_str"](0)
            error = np.max(np.linalg.norm(fitted - original, axis=2))
            self.assertLess(error, 1e-10)

            lines = iges.read_text().splitlines()
            self.assertTrue(lines)
            self.assertTrue(all(len(line) == 80 for line in lines))
            self.assertTrue(any(line.startswith("     128") for line in lines))
            self.assertEqual(lines[-1][72], "T")
            self.assertIn("4HSLOT,1.,6,1HM,", lines[3][:72])

    def test_output_units_feet(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            p3d = tmp / "surface.xyz"
            iges = tmp / "surface_ft.igs"
            _write_surface_plot3d(p3d, _make_surface())
            fit_plot3d_to_iges(p3d, iges, units="ft")
            lines = iges.read_text().splitlines()
            self.assertIn("4HSLOT,1.,4,2HFT,", lines[3][:72])
            self.assertTrue(all(len(line) == 80 for line in lines))

    def test_input_feet_to_output_meters_scales_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            p3d = tmp / "surface_ft.xyz"
            iges = tmp / "surface_m.igs"
            surface_ft = _make_surface()
            _write_surface_plot3d(p3d, surface_ft)

            bse = fit_plot3d_to_iges(p3d, iges, input_units="ft", units="m")
            fitted_m = bse.vec["pt_str"](0)
            expected_m = surface_ft * 0.3048
            error = np.max(np.linalg.norm(fitted_m - expected_m, axis=2))
            self.assertLess(error, 1e-10)
            self.assertIn("4HSLOT,1.,6,1HM,", iges.read_text().splitlines()[3][:72])


if __name__ == "__main__":
    unittest.main()
