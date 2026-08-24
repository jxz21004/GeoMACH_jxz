# Python 3 migration status

## Integrated later files

- Later `BSEmodel.py` and `BSEvec.py` were used as the Python 3 baseline.
- Later Plot3D/IGES utility code was integrated under `GeoMACH.utilities`.

## Python 3 compatibility fixes

- `xrange` -> `range`.
- Python 2 `print` statements -> `print()`.
- Python 2 implicit relative imports -> explicit package-relative imports.
- Python 2 dictionary-view indexing -> `list(...values())` where sequence
  semantics are required.
- Mixed tabs/spaces normalized without changing block structure.
- Integer identity comparisons (`is 0`) changed to value comparisons.
- `MACHconfiguration` is now an optional MPI import so base geometry use does
  not require mpi4py.

## Build modernization

The six historical F2PY extension modules are preserved:

1. `GeoMACH.BSE.BSElib`
2. `GeoMACH.PGM.PGMlib`
3. `GeoMACH.PSM.PSMlib`
4. `GeoMACH.PSM.QUADlib`
5. `GeoMACH.PSM.CDTlib`
6. `GeoMACH.PSM.BLSlib`

`tools/build_extensions.py` replaces the removed `numpy.distutils` build path
and has been exercised with Python 3.13, NumPy 2.x and gfortran.

## Plot3D -> IGES changes

The numerical conversion is intentionally the same later-GeoMACH route:
structured surface points -> BSE sparse maps -> LSQR control-point solve ->
IGES entity 128.

The Plot3D reader was made robust without reducing capability:

- no fixed `200 x 200 x 1 x 1000` allocations;
- arbitrary whitespace/line wrapping;
- Fortran `D` exponents accepted;
- block-specific dimensions and dynamic allocation;
- explicit rejection of volume blocks (`nk != 1`) for surface conversion;
- historical reader and patch class names retained.

## Validation completed in this migration snapshot

- All Python sources compile under Python 3.
- All six original Fortran extension groups compile and import on the test
  Linux environment.
- Historical `wing.py`, `conventional.py`, `supersonic.py`,
  `trussbraced_wingstrut.py`, `trussbraced_full.py`, and
  `RC_aircraft_tutorial.py` execute under Python 3 in the test environment.
- A two-block formatted Plot3D test with different block dimensions converts
  through BSE/LSQR to IGES; reconstruction error was about `6e-14`.

## Still to validate against real user data

- `cgnsUtilities` conversion of the user's actual `surface.cgns` to surface
  Plot3D.
- CAD import of the resulting real-aircraft `.igs` in the user's downstream
  CAD/viewer workflow.
- Full numerical/derivative regression against an executable Python 2 GeoMACH
  environment, if one is available.
