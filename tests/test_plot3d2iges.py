import tempfile
import unittest
from pathlib import Path

import numpy as np

from GeoMACH.utilities import P3dSurfClass

try:
    from GeoMACH.utilities.plot3d2iges import fit_plot3d_to_iges
except ImportError:
    fit_plot3d_to_iges = None


@unittest.skipIf(fit_plot3d_to_iges is None, "GeoMACH compiled extensions are not built")
class TestPlot3D2IGES(unittest.TestCase):
    def test_fit_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            p3d = tmp / "surface.xyz"
            iges = tmp / "surface.igs"

            ni, nj = 5, 4
            surface = np.empty((ni, nj, 3), dtype=float, order="F")
            for j in range(nj):
                for i in range(ni):
                    u = i / (ni - 1)
                    v = j / (nj - 1)
                    surface[i, j] = [u, v, 0.1 * u * (1.0 - u) + 0.05 * v]

            vals = ["1", str(ni), str(nj), "1"]
            for idim in range(3):
                vals.extend(f"{x:.16e}" for x in surface[..., idim].reshape(-1, order="F"))
            p3d.write_text("\n".join(vals) + "\n")

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


if __name__ == "__main__":
    unittest.main()
