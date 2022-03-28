#!/usr/bin/python

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SoundForCp",
    version="0.0.2",
    author="JunPei-China",
    author_email="Jun_Pei@163.com",
    description="A software package for calculating specific heat capacity based on sound velocity method.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JunPei-China/SoundForCp",
    project_urls={
        "Bug Tracker": "https://github.com/JunPei-China/SoundForCp/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    
    install_requires=['numpy>=1.20','scipy>=1.7','PyYaml>=6.0'],
    entry_points={
    'console_scripts': [
        'SoundForCp=SoundForCp.SoundForCp:SoundForCp'] },
)