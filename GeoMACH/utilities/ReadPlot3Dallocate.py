"""Formatted multi-block Plot3D reader used by the GeoMACH surface utilities.

This module keeps the historical ``readplot3d`` API while removing the fixed
200x200x1x1000 temporary allocation from the later utility code.  Coordinates
are read block-by-block, with i varying fastest as required by Plot3D.
"""

from __future__ import annotations

from itertools import islice
from pathlib import Path

import numpy as np


def _tokens(filename):
    """Yield whitespace-delimited Plot3D tokens without loading the file twice."""
    with Path(filename).open("r") as stream:
        for line in stream:
            # Accept Fortran D exponents as well as normal E exponents.
            for token in line.replace("D", "E").replace("d", "e").split():
                yield token


def _read_values(token_iter, count, dtype=float):
    values = np.fromiter((dtype(v) for v in islice(token_iter, count)), dtype=dtype, count=count)
    if values.size != count:
        raise ValueError(f"Unexpected end of Plot3D file: expected {count} values, got {values.size}.")
    return values


class readplot3d:
    """Read a formatted multi-block Plot3D grid.

    The legacy attributes ``nblk``, ``npts``, ``ni``, ``nj``, ``nk``, ``pts``
    and ``ptsindex`` are preserved.  ``blocks`` additionally exposes one
    ``(ni, nj, nk, 3)`` NumPy array per block.
    """

    def __init__(self):
        self.nblk = 0
        self.npts = 0
        self.ni = []
        self.nj = []
        self.nk = []
        self.pts = None
        self.ptsindex = None
        self.blocks = []

    def getmesh(self, meshfilename):
        tokens = _tokens(meshfilename)
        try:
            self.nblk = int(next(tokens))
        except StopIteration as err:
            raise ValueError(f"Empty Plot3D file: {meshfilename}") from err

        if self.nblk <= 0:
            raise ValueError(f"Invalid Plot3D block count: {self.nblk}")

        dims = _read_values(tokens, 3 * self.nblk, dtype=int).reshape(self.nblk, 3)
        if np.any(dims <= 0):
            raise ValueError(f"Plot3D block dimensions must be positive; got {dims.tolist()}")

        self.ni = dims[:, 0].astype(int).tolist()
        self.nj = dims[:, 1].astype(int).tolist()
        self.nk = dims[:, 2].astype(int).tolist()
        self.npts = int(np.prod(dims, axis=1).sum())

        self.blocks = []
        pts = np.empty((self.npts, 3), dtype=float)
        ptsindex = np.empty((self.npts, 4), dtype=int)
        offset = 0

        for iblk, (ni, nj, nk) in enumerate(dims):
            ni, nj, nk = int(ni), int(nj), int(nk)
            nblock = ni * nj * nk
            block = np.empty((ni, nj, nk, 3), dtype=float, order="F")
            for idim in range(3):
                coord = _read_values(tokens, nblock, dtype=float)
                block[..., idim] = coord.reshape((ni, nj, nk), order="F")
            self.blocks.append(block)

            flat = block.reshape((nblock, 3), order="F")
            pts[offset : offset + nblock] = flat
            local = np.indices((ni, nj, nk), dtype=int).reshape(3, -1, order="F").T + 1
            ptsindex[offset : offset + nblock, :3] = local
            ptsindex[offset : offset + nblock, 3] = iblk + 1
            offset += nblock

        # Extra non-whitespace tokens almost always indicate a dimension mismatch.
        try:
            extra = next(tokens)
        except StopIteration:
            extra = None
        if extra is not None:
            raise ValueError(f"Plot3D file contains extra data after the expected coordinates (first extra token: {extra!r}).")

        self.pts = pts
        self.ptsindex = ptsindex
        return self
