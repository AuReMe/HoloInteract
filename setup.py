#!/usr/bin/env python3
from sys import version_info, stderr
from setuptools import setup

NAME = 'HoloInteract'
CURRENT_PYTHON = version_info[:2]
REQUIRED_PYTHON = (3, 8)

if CURRENT_PYTHON < REQUIRED_PYTHON:
    stderr.write(
        f"{NAME} requires Python 3.10 or higher and your current version is {CURRENT_PYTHON}.")
    exit(1)


setup(
    name=NAME,
    version='0.1.0',
    description='',
    url='https://github.com/CoLucas22',
    author='CoLucas22',
    author_email='corentin.lucas22@outlook.fr',
    packages=["holointeract"],
    zip_safe=False,
    license="LICENSE",
    long_description=open("README.md", encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=["bubbletools", "metage2metabo", "kaleido", "ete3", "pandas", "scipy", "matplotlib", "plotly",
                      "statsmodels", "fastcluster", "seaborn", "rich", "emapper2gbk"],
    entry_points={'console_scripts': ['holointeract=holointeract.main:main']}
)
