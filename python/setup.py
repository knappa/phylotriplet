from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os
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


ext_modules = [
    Pybind11Extension(
        "phylotriplet",
        ["bindings.cpp"],
        include_dirs=["../include", find_eigen3()],
        cxx_std=14,
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy"],
)
