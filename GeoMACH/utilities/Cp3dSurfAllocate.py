"""Surface-patch helpers for formatted Plot3D grids.

The public class names are retained from the later GeoMACH utility code so
existing scripts can keep using them.
"""

from __future__ import annotations

import numpy as np

from .ReadPlot3Dallocate import readplot3d


class P3dSurfClass:
    """Represent a multi-block Plot3D surface grid."""

    def __init__(self, meshfilename):
        mesh = readplot3d().getmesh(meshfilename)
        self.meshfilename = meshfilename
        self.nsp = mesh.nblk
        self.npts = mesh.npts
        self.ni = list(mesh.ni)
        self.nj = list(mesh.nj)
        self.nk = list(mesh.nk)
        self.pts = mesh.pts
        self.ptsindex = mesh.ptsindex
        self.blocks = mesh.blocks

        # These must be instance attributes; the historical class attributes
        # made different P3dSurfClass instances share patches accidentally.
        self.surf3d = []
        self.surf2d = []
        self.surf = []

    def _require_surface_blocks(self):
        bad = [(i + 1, nk) for i, nk in enumerate(self.nk) if nk != 1]
        if bad:
            raise ValueError(
                "IGES/surface conversion requires Plot3D surface blocks with nk=1; "
                f"non-surface blocks: {bad}"
            )

    def Surf3dGen(self):
        self._require_surface_blocks()
        self.surf3d = [SurfPatch3D(block, iblk) for iblk, block in enumerate(self.blocks)]
        return self.surf3d

    def Surf2dGen(self):
        self._require_surface_blocks()
        self.surf2d = [SurfPatch2D(block[:, :, 0, :], iblk) for iblk, block in enumerate(self.blocks)]
        return self.surf2d

    def SurfGen(self):
        self._require_surface_blocks()
        self.surf = [np.array(block[:, :, 0, :], copy=True, order="F") for block in self.blocks]
        return self.surf

    def surf_update(self):
        """Synchronize ``pts`` from the current arrays in ``surf``."""
        if not self.surf:
            raise RuntimeError("SurfGen() must be called before surf_update().")
        offset = 0
        for iblk, surface in enumerate(self.surf):
            nblock = self.ni[iblk] * self.nj[iblk]
            self.pts[offset : offset + nblock] = surface.reshape((nblock, 3), order="F")
            offset += nblock


class SurfPatch3D:
    def __init__(self, block, sind):
        block = np.asarray(block)
        if block.ndim != 4 or block.shape[-1] != 3:
            raise ValueError("SurfPatch3D expects an (ni, nj, nk, 3) coordinate array.")
        if block.shape[2] != 1:
            raise ValueError("SurfPatch3D requires nk=1.")
        self.ind = sind + 1
        self.ni, self.nj, self.nk = block.shape[:3]
        self.pts = np.array(block, copy=True, order="F")
        self.dvs = np.zeros_like(self.pts)

    def dvs_update(self, dvs):
        dvs = np.asarray(dvs)
        if dvs.shape != self.pts.shape:
            raise ValueError(f"dvs shape {dvs.shape} does not match patch shape {self.pts.shape}.")
        self.dvs = dvs

    def surf_update(self):
        self.pts = self.pts + self.dvs


class SurfPatch2D:
    def __init__(self, surface, sind):
        surface = np.asarray(surface)
        if surface.ndim != 3 or surface.shape[-1] != 3:
            raise ValueError("SurfPatch2D expects an (ni, nj, 3) coordinate array.")
        self.ind = sind + 1
        self.ni, self.nj = surface.shape[:2]
        self.pts = np.array(surface, copy=True, order="F")
        self.dvs = np.zeros_like(self.pts)

    def dvs_update(self, dvs):
        dvs = np.asarray(dvs)
        if dvs.shape != self.pts.shape:
            raise ValueError(f"dvs shape {dvs.shape} does not match patch shape {self.pts.shape}.")
        self.dvs = dvs

    def surf_update(self):
        self.pts = self.pts + self.dvs
