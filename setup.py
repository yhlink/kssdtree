from setuptools import setup, Extension, find_packages
from os import environ
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
import sys, os, subprocess


class CustomBuildExt(build_ext):
    def run(self):
        super().run()
        self.fix_rpath()

    def fix_rpath(self):
        if sys.platform != "darwin":
            return

        conda_lib = os.path.join(sys.prefix, "lib")

        for ext in self.extensions:
            so_path = self.get_ext_fullpath(ext.name)
            if not os.path.exists(so_path):
                continue

            subprocess.run(
                ["install_name_tool", "-delete_rpath", conda_lib, so_path],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["install_name_tool", "-delete_rpath", "@loader_path", so_path],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["install_name_tool", "-delete_rpath", "@loader_path", so_path],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )


extra_compile_args = []
extra_link_args = []
if 'darwin' in sys.platform:
    extra_compile_args = ['-fopenmp']
    extra_link_args = ['-fopenmp']
    # clang
    # os.environ["CC"] = '/usr/bin/gcc'
    os.environ["CC"] = '/opt/homebrew/bin/gcc-12'

else:
    if environ.get('CC') and 'clang' in environ['CC']:
        # clang
        extra_compile_args = ['-fopenmp=libomp']
        extra_link_args = ['-fopenmp=libomp']
    else:
        # GNU
        extra_compile_args = ['-fopenmp']
        extra_link_args = ['-fopenmp']

MOD1 = 'kssd'
MOD2 = 'nj'
MOD3 = 'dnj'
sources1 = ['co2mco.c',
            'iseq2comem.c',
            'command_dist_wrapper.c',
            'mytime.c',
            'global_basic.c',
            'command_set.c',
            'command_dist.c',
            'command_shuffle.c',
            'command_composite.c',
            'mman.c',
            'pykssd.c']
sources2 = ['align.c',
            'cluster.c',
            'distancemat.c',
            'util.c',
            'tree.c',
            'buildtree.c',
            'sequence.c',
            'pynj.c']

sources3 = ['bytescale.c',
            'dnj.c',
            'str.c',
            'tmp.c',
            'phy.c',
            'filebuff.c',
            'nj.c',
            'qseqs.c',
            'vector.c',
            'matrix.c',
            'mman.c',
            'hclust.c',
            'nwck.c',
            'pherror.c',
            'pydnj.c']
include_dirs1 = ['kssdheaders']
include_dirs2 = ['njheaders']
include_dirs3 = ['dnjheaders']

require_pakages = [
    'pyqt5',
    'ete3',
    'requests',
    'pandas'
]

setup(
    name='kssdtree',
    version='2.1.0',
    author='Hang Yang',
    author_email='yhlink1207@gmail.com',
    description="Kssdtree is a versatile Python package for phylogenetic analysis. It also provides one-stop tree construction and visualization. It can handle DNA sequences of both fasta or fastq format, whether gzipped or not. ",
    url='https://github.com/yhlink/kssdtree',
    download_url='https://pypi.org/project/kssdtree',
    ext_modules=[
        Extension(MOD1, sources=sources1, include_dirs=include_dirs1, libraries=['z'],
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args),
        Extension(MOD2, sources=sources2, include_dirs=include_dirs2),
        Extension(MOD3, sources=sources3, include_dirs=include_dirs3, libraries=['z'],
                  extra_compile_args=extra_compile_args,
                  extra_link_args=extra_link_args)
    ],
    py_modules=['kssdtree', 'toolutils'],
    packages=find_packages(),
    install_requires=require_pakages,
    dependency_links=['https://pypi.python.org/simple/'],
    zip_safe=False,
    include_package_data=True,
    cmdclass={
        "build_ext": CustomBuildExt,
    }
)