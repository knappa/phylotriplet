from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
import platform
import subprocess


def find_eigen3():
    """Return the include directory containing Eigen/Dense."""
    # 1. Explicit override via environment variable
    env_dir = os.environ.get("EIGEN3_INCLUDE_DIR")
    if env_dir:
        if os.path.exists(os.path.join(env_dir, "Eigen", "Dense")):
            return env_dir
        raise RuntimeError(
            f"EIGEN3_INCLUDE_DIR is set to '{env_dir}' but Eigen/Dense was not found there"
        )

    # 2. pkg-config (works on Linux and macOS with Homebrew/MacPorts)
    try:
        result = subprocess.run(
            ["pkg-config", "--cflags-only-I", "eigen3"],
            capture_output=True, text=True, check=True,
        )
        for flag in result.stdout.split():
            if flag.startswith("-I"):
                path = flag[2:]
                if os.path.exists(os.path.join(path, "Eigen", "Dense")):
                    return path
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    # 3. Common installation prefixes
    candidates = [
        "/usr/include/eigen3",           # Debian/Ubuntu
        "/usr/local/include/eigen3",     # manual install
        "/opt/homebrew/include/eigen3",  # Homebrew on Apple Silicon
        "/opt/local/include/eigen3",     # MacPorts
        r"C:\Program Files\Eigen3\include\eigen3",  # Windows (vcpkg / manual)
    ]
    for path in candidates:
        if os.path.exists(os.path.join(path, "Eigen", "Dense")):
            return path

    raise RuntimeError(
        "Eigen3 not found. Install it (e.g. 'apt install libeigen3-dev', "
        "'brew install eigen') or set the EIGEN3_INCLUDE_DIR environment variable."
    )


def openmp_flags():
    """Return (compile_args, link_args, extra_include_dirs) for OpenMP.

    On macOS, Apple Clang requires libomp installed via Homebrew and uses
    -Xpreprocessor -fopenmp / -lomp rather than -fopenmp / -lgomp.
    """
    if platform.system() == "Darwin":
        # Locate libomp via brew; fall back to common prefixes.
        omp_prefix = None
        try:
            result = subprocess.run(
                ["brew", "--prefix", "libomp"],
                capture_output=True, text=True, check=True,
            )
            omp_prefix = result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            for path in ["/opt/homebrew", "/usr/local"]:
                if os.path.exists(os.path.join(path, "opt", "libomp", "include")):
                    omp_prefix = os.path.join(path, "opt", "libomp")
                    break
        if omp_prefix is None:
            raise RuntimeError(
                "libomp not found. Install it with: brew install libomp"
            )
        inc = os.path.join(omp_prefix, "include")
        lib = os.path.join(omp_prefix, "lib")
        compile_args = ["-Xpreprocessor", "-fopenmp", "-I" + inc]
        link_args = ["-L" + lib, "-lomp"]
        return compile_args, link_args
    else:
        return ["-fopenmp"], ["-fopenmp"]


omp_compile, omp_link = openmp_flags()

ext_modules = [
    Pybind11Extension(
        "phylotriplet",
        ["bindings.cpp"],
        include_dirs=["../include", "../third_party/LBFGSpp", find_eigen3()],
        extra_compile_args=omp_compile,
        extra_link_args=omp_link,
        cxx_std=14,
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy"],
)
