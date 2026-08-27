from pathlib import Path
import tempfile

import numpy as np
from .plot3d2iges import fit_plot3d_to_iges


def cgns2iges(
    cgns_file,
    iges_file,
    grid_type="auto",
    units="m",
    input_units=None,
    debug=False,
):
    """
    Convert a structured CGNS surface/volume grid to IGES.
    Notice that if CGNS is volume file, we only extract wall surfaces.
    
    Parameters
    ----------
    cgns_file : str
        Input structured CGNS file.

    iges_file : str
        Output IGES file.

    grid_type : {"auto", "surface", "volume"}
        "surface": CGNS itself contains surface zones.
        "volume" : extract wall surfaces from the volume grid.
        "auto"   : determine from block dimensions.

    units : str
        IGES output units.

    input_units : str or None
        CGNS coordinate units. If None, assumed equal to `units`.

    debug : bool
        Print GeoMACH fitting information.
    """
    try:
        from cgnsutilities.cgnsutilities import readGrid
    except ImportError as err:
        raise ImportError(
            "cgns2iges requires cgnsUtilities. "
            "Install cgnsUtilities to use CGNS-to-IGES conversion."
        ) from err
    grid = readGrid(str(cgns_file))

    if grid_type == "auto":
        is_surface = all(
            np.any(np.asarray(block.dims) == 1)
            for block in grid.blocks
        )
        grid_type = "surface" if is_surface else "volume"

    if grid_type not in ("surface", "volume"):
        raise ValueError(
            "grid_type must be 'auto', 'surface', or 'volume'"
        )

    with tempfile.TemporaryDirectory() as tmp:
        plot3d_file = Path(tmp) / "surface.xyz"

        if grid_type == "surface":
            # CGNS surface zones -> Plot3D
            grid.writePlot3d(str(plot3d_file))

        else:
            # Volume CGNS -> wall surface Plot3D
            grid.extractSurface(str(plot3d_file))

        return fit_plot3d_to_iges(
            plot3d_filename=str(plot3d_file),
            iges_filename=str(iges_file),
            units=units,
            input_units=input_units,
            export_debug=debug,
        )

def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert a structured CGNS surface/volume grid to IGES."
    )

    parser.add_argument(
        "--cgns",
        help="Input CGNS file.",
    )
    parser.add_argument(
        "--iges",
        help="Output IGES file.",
    )
    parser.add_argument(
        "--grid-type",
        choices=("auto", "surface", "volume"),
        default="auto",
        help="CGNS grid type (default: auto).",
    )
    parser.add_argument(
        "--units",
        default="m",
        help="IGES output units (default: m).",
    )
    parser.add_argument(
        "--input-units",
        default=None,
        help="Input CGNS coordinate units (default: same as output units).",
    )
    parser.add_argument(
        "--debug-output",
        action="store_true",
        help="Print GeoMACH fitting diagnostics.",
    )

    args = parser.parse_args(argv)

    cgns2iges(
        cgns_file=args.cgns,
        iges_file=args.iges,
        grid_type=args.grid_type,
        units=args.units,
        input_units=args.input_units,
        debug=args.debug_output,
    )


if __name__ == "__main__":
    main()