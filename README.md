# GeoMACH Python 3 migration

This tree is a compatibility-focused Python 3 migration of the historical
GeoMACH source.  The BSE, PGM and PSM algorithms and all six original Fortran
extension groups are retained.

## Build

The old source used `numpy.distutils`, which is no longer available in modern
NumPy.  Following the current MDO Lab style used by compiled utility packages,
compile the in-place extensions first and then install the Python package:

```bash
make
pip install -e .
```

The build script uses `CC` and `FC` from the environment when supplied and
otherwise defaults to `cc` and `gfortran`.

MPI coupling through `MACHconfiguration` additionally requires:

```bash
pip install -e '.[mpi]'
```

## Surface Plot3D to IGES

A later GeoMACH utility path has been integrated into the package.  For a
formatted multi-block surface Plot3D file (each block must have `nk=1`):

```bash
geomach_plot3d2iges --plot3d surface.xyz --iges surface.igs
```

or

```bash
python -m GeoMACH.utilities.plot3d2iges --plot3d surface.xyz --iges surface.igs
```

The conversion constructs a BSE model, solves the sparse B-spline fitting
problem with LSQR, and writes IGES entity type 128 B-spline surfaces.

IGES output units default to meters:

```bash
geomach_plot3d2iges --plot3d surface.xyz --iges surface.igs
```

To declare another output unit without rescaling the numeric coordinates:

```bash
geomach_plot3d2iges --plot3d surface.xyz --iges surface_ft.igs --units ft
```

If the Plot3D numeric coordinates are in a different physical unit, specify
both units so the coordinates are rescaled before the BSE fit. For example,
feet to meters:

```bash
geomach_plot3d2iges --plot3d surface_ft.xyz --iges surface_m.igs \
    --input-units ft --units m
```

Supported IGES units are `in`, `mm`, `ft`, `mi`, `m`, `km`, `mil`, `um`,
`cm`, and `uin`. When `--input-units` is omitted, Plot3D coordinates are
assumed to already use the selected output unit.

For `surface.cgns`, first use MDO Lab `cgnsUtilities` to export a formatted
surface Plot3D file, then run the command above.  Keeping CGNS parsing in
`cgnsUtilities` avoids duplicating the CGNS implementation inside GeoMACH.

## Migration policy

The migration deliberately avoids replacing GeoMACH with pyGeo or reducing it
to an IGES converter.  BSE topology/derivatives, PGM aircraft components,
PSM structural modeling, QUAD/CDT/BLS meshing helpers, examples, and historical
export functions remain in the source tree.
