from glob import glob
from pybind11.setup_helpers import Pybind11Extension
from pathlib import Path
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os
import sys
import re
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = os.path.join(lib_folder, '/requirements.txt')
install_requires = []
if os.path.isfile(requirement_path):
    with open(requirement_path, 'r', encoding='utf-8') as f:
        install_requires = f.read().splitlines()

CMAKE_FLAGS = []
if "CMAKE_ARGS" in os.environ:
    CMAKE_FLAGS = [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

INSTALL_ARGS = []
if "CMAKE_INSTALL_PREFIX" in os.environ:
    INSTALL_ARGS = ['--prefix', os.environ["CMAKE_INSTALL_PREFIX"]]

class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())
        print(self.sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext) -> None:
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)
        
        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg   = "Debug" if debug else "Release"
        
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ] + CMAKE_FLAGS

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        subprocess.run(
            ["cmake", ext.sourcedir, *cmake_args], cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", '--build', '.'], cwd=build_temp, check=True
        )

setup(
    name='fierro_voxelizer_py',
    python_requires='>=3.8',
    author='Kevin Welsh',
    author_email='kwelsh@lanl.gov',
    version='0.2',
    description='',
    include_package_data=True,
    install_requires=install_requires,
    ext_modules=[CMakeExtension('voxerlizer_cpp')],
    cmdclass={"build_ext" : CMakeBuild},
    zip_safe=False,
)
