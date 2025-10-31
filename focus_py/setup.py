from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11
import os

# adjust these if your system installs Qhull elsewhere
QHULL_INCLUDE = "/usr/include/libqhullcpp"
QHULL_LIB_DIR = None  # e.g. "/usr/lib" or "/usr/lib/x86_64-linux-gnu" if needed

extra_kwargs = {}
if QHULL_LIB_DIR:
    extra_kwargs["library_dirs"] = [QHULL_LIB_DIR]
    # ensure runtime rpath so the shared object can find qhull at runtime
    extra_kwargs["runtime_library_dirs"] = [QHULL_LIB_DIR]

ext_modules = [
    Pybind11Extension(
        "focus",
        ["pybind11_focus.cpp"],            # make sure filename matches
        include_dirs=[
            pybind11.get_include(),
            ".", 
            QHULL_INCLUDE,                 # <---- add Qhull headers
        ],
        libraries=["qhullcpp", "qhull_r"], # <---- link these libraries
        extra_compile_args=["-O3"],
        cxx_std=17,
        **extra_kwargs
    ),
]

setup(
    name="focus",
    version="0.1.0",
    author="Your Name",
    description="FOCuS changepoint detection library",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=["numpy>=1.19.0"],
)
