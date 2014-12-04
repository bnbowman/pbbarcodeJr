from setuptools import setup, Extension, find_packages
import os
import sys

vFile = 'src/python/pbbarcode/_version.py'

if os.path.exists(vFile):
    lines = open(vFile, 'r').read().splitlines()
    for line in lines:
        elts = line.split('=')
        elts = [e.strip() for e in elts]
        if len(elts) == 2 and elts[0] == '__version__':
            _ReadVersion = elts[1].replace('\'', '').replace('\"', '')
            break
else:
    _ReadVersion = '0.0.0'

setup(
    name = 'pbbarcodejr',
    version=_ReadVersion,
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    packages = find_packages('src/python'),
    package_dir = {'':'src/python'},
    ext_modules=[Extension('pbbarcodejr/sw', ['src/C/sw.c'], extra_compile_args=["-O3","-shared"]),
                 Extension('pbbarcodejr/sw_new', ['src/C/sw_new.c'], extra_compile_args=["-O3","-shared"])],
    zip_safe = False,
    entry_points={
        'console_scripts': [
            'pbbarcodejr = pbbarcodejr.main:main']
    },
    install_requires=[
        'pbcore >= 0.6.3',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0'
        ]
    )
