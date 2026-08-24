import tempfile
import unittest
from pathlib import Path

import numpy as np

from GeoMACH.utilities import P3dSurfClass, readplot3d


class TestPlot3DIO(unittest.TestCase):
    def _write_surface(self, path):
        blocks = []
        for iblk, (ni, nj) in enumerate([(3, 2), (2, 4)]):
            block = np.empty((ni, nj, 1, 3), dtype=float, order="F")
            for j in range(nj):
                for i in range(ni):
                    block[i, j, 0] = [10 * iblk + i, j, 100 + 10 * iblk + i + 0.1 * j]
            blocks.append(block)

        values = ["2"]
        for block in blocks:
            values.extend([str(block.shape[0]), str(block.shape[1]), "1"])
        for block in blocks:
            for idim in range(3):
                values.extend(f"{value:.12E}" for value in block[..., idim].reshape(-1, order="F"))

        with path.open("w") as stream:
            for i in range(0, len(values), 7):
                stream.write(" ".join(values[i : i + 7]) + "\n")
        return blocks

    def test_multiblock_arbitrary_line_wrapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "surface.xyz"
            reference = self._write_surface(path)
            mesh = readplot3d().getmesh(path)
            self.assertEqual(mesh.nblk, 2)
            self.assertEqual(mesh.ni, [3, 2])
            self.assertEqual(mesh.nj, [2, 4])
            self.assertEqual(mesh.nk, [1, 1])
            for actual, expected in zip(mesh.blocks, reference):
                np.testing.assert_allclose(actual, expected)

            surf = P3dSurfClass(path)
            surfaces = surf.SurfGen()
            for actual, expected in zip(surfaces, reference):
                np.testing.assert_allclose(actual, expected[:, :, 0, :])


if __name__ == "__main__":
    unittest.main()
