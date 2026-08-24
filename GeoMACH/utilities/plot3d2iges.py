"""Fit a GeoMACH BSE model to a surface Plot3D grid and export IGES."""

from __future__ import annotations

import argparse

import numpy as np
from scipy.sparse.linalg import lsqr

from GeoMACH.BSE.BSEmodel import BSEmodel
from .Cp3dSurfAllocate import P3dSurfClass


def fit_plot3d_to_iges(plot3d_filename, iges_filename, export_debug=False, lsqr_kwargs=None):
    """Convert a formatted surface Plot3D file to an IGES B-spline model.

    This preserves the later GeoMACH conversion algorithm: the structured
    Plot3D points initialize ``pt_str``; BSE's assembled sparse mappings define
    the unique surface-point/control-point relation; and LSQR solves the three
    coordinate systems for the B-spline control points.
    """
    surf_p3d = P3dSurfClass(plot3d_filename)
    surfaces = surf_p3d.SurfGen()
    bsurf = BSEmodel(surfaces)

    for isurf, surface in enumerate(surfaces):
        nu, nv = surface.shape[:2]
        bsurf.set_bspline_option("num_pt", isurf, "u", nu)
        bsurf.set_bspline_option("num_pt", isurf, "v", nv)
        bsurf.set_bspline_option("num_cp", isurf, "u", max(4, nu))
        bsurf.set_bspline_option("num_cp", isurf, "v", max(4, nv))

    bsurf.assemble()

    for isurf, surface in enumerate(surfaces):
        bsurf.vec["pt_str"](isurf)[:, :, :] = surface

    bsurf.apply_jacobian("pt", "d(pt)/d(pt_str)", "pt_str")

    A = bsurf.jac["d(pt)/d(pt_str)"].dot(
        bsurf.jac["d(pt_str)/d(cp_str)"].dot(bsurf.jac["d(cp_str)/d(cp)"])
    )
    rhs = bsurf.vec["pt"].array
    cp = bsurf.vec["cp"].array

    kwargs = {"atol": 1.0e-12, "btol": 1.0e-12} if lsqr_kwargs is None else dict(lsqr_kwargs)
    for idim in range(3):
        cp[:, idim] = lsqr(A, rhs[:, idim], **kwargs)[0]

    bsurf.apply_jacobian("cp_str", "d(cp_str)/d(cp)", "cp")
    bsurf.apply_jacobian("pt_str", "d(pt_str)/d(cp_str)", "cp_str")

    if export_debug:
        bsurf.vec["pt_str"].export_tec_str()
        bsurf.vec["cp"].export_tec_scatter()

    bsurf.vec["cp_str"].export_IGES(filename=iges_filename)
    return bsurf


def _build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plot3d",
        required=True,
        help="Input formatted surface Plot3D file. Every block must have nk=1.",
    )
    parser.add_argument("--iges", required=True, help="Output IGES/IGS surface file.")
    parser.add_argument("--debug-output", action="store_true", help="Also export fitted Tecplot/control-point files.")
    return parser


def main(argv=None):
    args = _build_parser().parse_args(argv)
    fit_plot3d_to_iges(args.plot3d, args.iges, export_debug=args.debug_output)


if __name__ == "__main__":
    main()
