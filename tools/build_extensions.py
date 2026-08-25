#!/usr/bin/env python3
"""Build GeoMACH's historical F2PY extensions without numpy.distutils.

The build layout follows the current MDO Lab pattern used by packages with
compiled helpers: compile the extension modules in-place first, then install
with setuptools/pip.  This script intentionally preserves the original source
partitioning into BSElib, PGMlib, PSMlib, QUADlib, CDTlib, and BLSlib.

Linux/gfortran is the primary supported path, matching the historical GeoMACH
and current MACH ecosystem deployment environment.  ``CC``, ``FC``, ``CFLAGS``
and ``FFLAGS`` can be overridden from the environment.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys
import sysconfig
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

EXTENSIONS = {
    "BSElib": (
        ROOT / "GeoMACH/BSE",
        [
            "src/BSE/compute_topology.f90",
            "src/BSE/compute_indices.f90",
            "src/BSE/compute_in_jacobian.f90",
            "src/BSE/compute_df_jacobian.f90",
            "src/BSE/compute_cp_jacobian.f90",
            "src/BSE/compute_pt_jacobian.f90",
            "src/BSE/compute_bs_jacobian.f90",
            "src/BSE/compute_sc_jacobian.f90",
            "src/BSE/compute_projection.f90",
            "src/BSE/bspline_knot.f90",
            "src/BSE/bspline_param.f90",
            "src/BSE/bspline_basis.f90",
        ],
    ),
    "PGMlib": (
        ROOT / "GeoMACH/PGM",
        [
            "src/PGM/parameter/compute_bspline.f90",
            "src/PGM/parameter/bspline_knot.f90",
            "src/PGM/parameter/bspline_param.f90",
            "src/PGM/parameter/bspline_basis.f90",
            "src/PGM/primitive/computeAngles.f90",
            "src/PGM/primitive/computeRotations.f90",
            "src/PGM/primitive/computeRtnMtx.f90",
            "src/PGM/primitive/computeSections.f90",
            "src/PGM/primitive/computeShape.f90",
            "src/PGM/interpolant/computeCone.f90",
            "src/PGM/interpolant/computeJunction.f90",
            "src/PGM/interpolant/computeTip.f90",
            "src/PGM/interpolant/interpolant.f90",
        ],
    ),
    "PSMlib": (
        ROOT / "GeoMACH/PSM",
        [
            "src/PSM/GFEM/computeProjtnInputs.f90",
            "src/PSM/GFEM/computePreviewSurfaces.f90",
            "src/PSM/GFEM/computeEdgeLengths.f90",
            "src/PSM/GFEM/computeFaceDimensions.f90",
            "src/PSM/GFEM/importMembers.f90",
            "src/PSM/GFEM/computePreviewMembers.f90",
            "src/PSM/GFEM/computeMemberTopology.f90",
            "src/PSM/GFEM/computeAdjoiningEdges.f90",
            "src/PSM/GFEM/computeFaceEdges.f90",
            "src/PSM/GFEM/computeGroupIntersections.f90",
            "src/PSM/GFEM/computeGroupSplits.f90",
            "src/PSM/GFEM/computeIntersectionVerts.f90",
            "src/PSM/GFEM/computeSurfaces.f90",
            "src/PSM/GFEM/computeSurfaceProjections.f90",
            "src/PSM/GFEM/computeMemberEdges.f90",
            "src/PSM/GFEM/computeMemberNodes.f90",
            "src/PSM/GFEM/computeMembers.f90",
            "src/PSM/GFEM/removeDuplicateNodes.f90",
            "src/PSM/GFEM/removeRightQuads.f90",
            "src/PSM/GFEM/identifySymmNodes.f90",
            "src/PSM/GFEM/misc.f90",
        ],
    ),
    "QUADlib": (
        ROOT / "GeoMACH/PSM",
        [
            "src/PSM/QUAD/importEdges.f90",
            "src/PSM/QUAD/reorderCollinear.f90",
            "src/PSM/QUAD/addIntersectionPts.f90",
            "src/PSM/QUAD/addEdgePts.f90",
            "src/PSM/QUAD/addInteriorPts.f90",
            "src/PSM/QUAD/splitEdges.f90",
            "src/PSM/QUAD/removeDegenerateEdges.f90",
            "src/PSM/QUAD/removeDuplicateEdges.f90",
            "src/PSM/QUAD/removeDuplicateQuads.f90",
            "src/PSM/QUAD/removeDuplicateTriangles.f90",
            "src/PSM/QUAD/removeDuplicateVerts.f90",
            "src/PSM/QUAD/computeAdjMap.f90",
            "src/PSM/QUAD/computeTriangles.f90",
            "src/PSM/QUAD/computeQuads.f90",
            "src/PSM/QUAD/computeConstraints.f90",
            "src/PSM/QUAD/computeQuadDominant.f90",
            "src/PSM/QUAD/splitTrisNQuads.f90",
            "src/PSM/QUAD/computeQuad2Edge.f90",
            "src/PSM/QUAD/removeInvalidQuads.f90",
        ],
    ),
    "CDTlib": (
        ROOT / "GeoMACH/PSM",
        [
            "src/PSM/CDT/addNode.f90",
            "src/PSM/CDT/computeCDT.f90",
            "src/PSM/CDT/constraints.f90",
            "src/PSM/CDT/delaunay.f90",
            "src/PSM/CDT/delete.f90",
            "src/PSM/CDT/misc.f90",
            "src/PSM/CDT/nearest.f90",
            "src/PSM/CDT/output.f90",
            "src/PSM/CDT/postProcess.f90",
        ],
    ),
    "BLSlib": (ROOT / "GeoMACH/PSM", ["src/PSM/BLS/assembleMtx.f90"]),
}


def _run(cmd, cwd=None, verbose=False):
    if verbose:
        print("+", shlex.join(str(v) for v in cmd))
    subprocess.run([str(v) for v in cmd], cwd=cwd, check=True)


def _compiler_command(env_name, default):
    return shlex.split(os.environ.get(env_name, default))


def build_extension(name, verbose=False):
    dest, rel_sources = EXTENSIONS[name]
    sources = [ROOT / src for src in rel_sources]
    cc = _compiler_command("CC", "cc")
    fc = _compiler_command("FC", "gfortran")
    cflags = shlex.split(os.environ.get("CFLAGS", "-O2 -fPIC"))
    fflags = shlex.split(os.environ.get("FFLAGS", "-O2 -fPIC"))

    py_include = Path(sysconfig.get_paths()["include"])
    np_include = Path(np.get_include())
    f2py_src = Path(np.__file__).resolve().parent / "f2py" / "src"
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"

    with tempfile.TemporaryDirectory(prefix=f"geomach-{name}-") as tmp_name:
        tmp = Path(tmp_name)
        _run(
            [sys.executable, "-m", "numpy.f2py", "--lower", *sources, "-m", name],
            cwd=tmp,
            verbose=verbose,
        )

        objects = []
        module_obj = tmp / f"{name}module.o"
        _run(
            [*cc, *cflags, "-c", tmp / f"{name}module.c", f"-I{py_include}", f"-I{np_include}", f"-I{f2py_src}", "-o", module_obj],
            verbose=verbose,
        )
        objects.append(module_obj)

        fortran_obj = tmp / "fortranobject.o"
        _run(
            [*cc, *cflags, "-c", f2py_src / "fortranobject.c", f"-I{py_include}", f"-I{np_include}", f"-I{f2py_src}", "-o", fortran_obj],
            verbose=verbose,
        )
        objects.append(fortran_obj)

        for index, source in enumerate(sources):
            obj = tmp / f"src_{index:03d}.o"
            _run([*fc, *fflags, "-c", source, "-o", obj], verbose=verbose)
            objects.append(obj)

        for wrapper in (tmp / f"{name}-f2pywrappers.f", tmp / f"{name}-f2pywrappers2.f90"):
            if wrapper.exists() and wrapper.stat().st_size:
                obj = tmp / f"{wrapper.stem}.o"
                _run([*fc, *fflags, "-c", wrapper, "-o", obj], verbose=verbose)
                objects.append(obj)

        output = dest / f"{name}{ext_suffix}"
        output.parent.mkdir(parents=True, exist_ok=True)
        _run([*fc, "-shared", "-o", output, *objects], verbose=verbose)
        print(f"Built {output.relative_to(ROOT)}")


def clean():
    for name, (dest, _) in EXTENSIONS.items():
        for path in dest.glob(f"{name}*.so"):
            path.unlink()
            print(f"Removed {path.relative_to(ROOT)}")


def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "extensions",
        nargs="*",
        metavar="EXTENSION",
        help=(
            "Extensions to build (default: all). "
            f"Available: {', '.join(sorted(EXTENSIONS))}"
        ),
    )
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args(argv)

    invalid = set(args.extensions) - set(EXTENSIONS)
    if invalid:
        parser.error(
            f"invalid extension(s): {', '.join(sorted(invalid))}. "
            f"Choose from: {', '.join(sorted(EXTENSIONS))}"
        )

    if args.clean:
        clean()
        return

    targets = args.extensions or list(EXTENSIONS)

    for name in targets:
        build_extension(name, verbose=args.verbose)


if __name__ == "__main__":
    main()
