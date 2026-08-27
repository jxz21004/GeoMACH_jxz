from setuptools import find_packages, setup

setup(
    name="GeoMACH",
    version="0.1.0",
    description=(
        "Geometry-centric MDAO of Aircraft Configurations with High fidelity: "
        "parametric geometry and structural modeling with analytic derivatives."
    ),
    url="https://github.com/hwangjt/GeoMACH",
    author="John Hwang and GeoMACH contributors",
    license="LGPL",
    packages=find_packages("."),
    include_package_data=True,
    package_data={
        "GeoMACH.PGM": ["airfoils/*.dat", "*.so"],
        "GeoMACH.BSE": ["*.so"],
        "GeoMACH.PSM": ["*.so"],
    },
    install_requires=["numpy>=1.25", "scipy>=1.11"],
    extras_require={"mpi": ["mpi4py>=3.1"]},
    python_requires=">=3.10",
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
    ],
    entry_points={
        "console_scripts": [
            "geomach_plot3d2iges=GeoMACH.utilities.plot3d2iges:main",
            "geomach_cgns2iges = GeoMACH.utilities.cgns2iges:main",
        ]
    },
)
