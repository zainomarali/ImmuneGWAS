# !/usr/bin/env python3
from setuptools import setup

setup(
    name='ImmuneGWAS',
    version='0.0.1',
    author='Zain Omar Ali et al',
    author_email='zain.ali@med.lu.se',
    packages=['ImmuneGWAS/'],
    scripts=[],
    url='https://github.com/zainomarali/ImmuneGWAS',
    license='LICENSE',
    description='Lookup utilities for GWAS hits',
    long_description=open('README.md').read(),
    install_requires=[
        'matplotlib',
        'openpyxl',
        'pandas <= 1.4.2',
        'pytabix==0.1',
        'pytest==7.1.1',
        'regex',
        'requests',
        'seaborn',
    ],
)
