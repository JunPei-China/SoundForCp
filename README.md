# SoundForCp

`SoundForCp` is a  software package for calculating specific heat capacity based on sound velocity method.  A detailed description is available in xxxx.

## Requirement

`SoundForCp` requires Python 3 and several third-party libraries that are used by SoundForCp.

- Numpy \- The fundamental package for scientific computing with Python. (>=1.20)
- Scipy \- Fundamental algorithms for scientific computing in Python for scientific computing with Python. (>=1.7)
- PyYaml \- A full-featured [YAML](http://yaml.org/) framework for the Python programming language.(>=6.0)
- Python (>=3.6)

As an example, on Ubuntu Linux some of the required software can be installed using

```bash
sudo apt-get install python3 python3-numpy python3-scipy python3-yaml
```

Before you go any further, make sure you have Python and that the expected version is available from your command line. You can check this by running:

Unix/MacOS:

```
python3 --version
```

Windows:

```
py --version
```

You should get some output like `Python 3.6.3`. If you do not have Python, please install the latest 3.x version from [python.org](https://www.python.org/) or refer to the [Installing Python](https://docs.python-guide.org/starting/installation/#installation) section of the Hitchhiker’s Guide to Python.



We recommend to use [PyPI](https://pypi.org/) which allows to conveniently install SoundForCp and all its software dependencies with a single command. you’ll need to make sure you have [pip](https://packaging.python.org/en/latest/key_projects/#pip) available. You can check this by running:

Unix/MacOS:

```
python3 -m pip --version
```

Windows:

```
py -m pip --version
```

If you installed Python from source, with an installer from [python.org](https://www.python.org/), or via [Homebrew](https://brew.sh/) you should already have pip. If you’re on Linux and installed using your OS package manager, you may have to install pip separately, see [Installing pip/setuptools/wheel with Linux Package Managers](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/).

If `pip` isn’t already installed, then first try to bootstrap it from the standard library:

Unix/MacOS:

```
python3 -m ensurepip --default-pip
```

Windows:

```
python3 -m ensurepip --default-pip
```

To install the remaining packages see the installation instructions at their respective web pages.

## Installation

The preferred method is to use `pip` command to install from the PyPi official/mirrors site.  

```
pip install SoundForCp
```

SounfForCp can be then started from a terminal (e.g. “Anaconda Prompt” on Windows) by executing the “SoundForCp” program.

If you don’t use `pip` or prefer to install from sources, make sure the required software is all in place and run:

```
python setup.py install
```

By default the files are installed to standard system directories, which may require the use of `sudo` for write privileges. If administrator (root) access is not available, see the output from `python setup.py install --help` for options to install as a regular user to user-writable locations. Note that installation to non-standard directories may require adjustments to the PATH and PYTHONPATH environment variables.

## Using

SoundForCp can be then started from a terminal (“Anaconda Prompt” on Windows) by executing the “SoundForCp.exe” program.  An alternative method on Windows is to start SoundForCp through double-click the executable program “SoundForCP.exe”.

Before running SoundForCp, a input file, named as `input.yaml` should have been prepared in a certain folder `xxxx`. `Input.yaml` follows the rule of [YAML language](https://yaml.org/spec/1.2.2/#chapter-2-language-overview). `input.yaml` contains `SoundForCp` settings which are summarized at setting tags.

```bash
cd xxxx
SoundForCp.exe
```

after running successfully, there is an output file named as “out.csv” written in the same folder.

## Setting tags

### Basic tags

#### `Sample_Name`

`Sample_Name` represents the sample name.

#### `Longitudinal_Sound_Velocity`

`Longitudinal_Sound_Velocity` represents the longitudinal Sound velocity of the specific sample. The unit is m/s.

#### `Transverse_Sound_Velocity`

`Transverse_Sound_Velocity` represents the transverse Sound velocity of the specific sample. The unit is m/s.

#### `Sample_Density`

`Sample_Density` represents the density of the specific sample, which can be measured by Archimedes method. The unit is g/cm3.

#### `Number_Atoms`

`Number_Atoms` represents the atom numbers per unit cell for the specific sample. 

#### `Volume_Cell`

`Volume_Cell` represents the lattice volume per unit cell for the specific sample. The unit is A^3.

#### `Temperature`

`Temperature` represents the calculated temperature range for the specific sample. there are 3 sub-tags, i.e.  `Start_Temperature`.  `End_Temperature`, `Interval_Temperature`. The unit is K.

#### `Relative_Atomic_Mass`

`Relative_Atomic_Mass` represents  relative atomic mass for molecular formula.  The unit is g/mol.

#### `Atomic_Mole_Quantity`

`Relative_Atomic_Mass` represents relative mole quantity for molecular formula.



### Required tags

#### `Elastic_Constants_Condition`

`Elastic_Constants_Condition` represent  elastic constants C12 equal C44 or not.  `Elastic_Constants_Condition`=2 represents  non-preferred orientation polycrystals. `Elastic_Constants_Condition`=1 represents  C12=C44. For cubic phase, C12 is equal to  C44  in preferred-oritention polycrystals or single crystal. For non cubic phase, C12 is unequal to  C44  in preferred-oritention polycrystals or single crystal. 



## Optional tags

#### `Adiabatic_Bulk_Modulus`

`Adiabatic_Bulk_Modulus` represents the adiabatic bulk modulus for the specific sample. It can also be calculated by sound velocity method. The unit is GPa.

#### `Linear_Expansion_Coefficient`

`Linear_Expansion_Coefficient` represents the linear expansion coefficient for the specific sample. It can also be calculated by sound velocity method. The unit is 1/K.

#### `Debye_Temperature`

`Debye_Temperature` represents the Debye temperature for the specific sample. It can also be calculated by sound velocity method. The unit is K. 

#### `Poisson_Ratio`

`Poisson_Ratio`  represents the Poisson ratio for the specific sample. It can also be calculated by sound velocity method. The unit is dimensionless. 

#### `Gruneisen_Constant`

`Gruneisen_Constant`  represents the Gruneisen constant for the specific sample. It can also be calculated by sound velocity method. The unit is dimensionless. 



## Development

PDFgui is an open-source software available on GitHub https://github.com/JunPei-China/SoundForCp. Feel free to fork the project and contribute. 

## Contacts

For more information on SoundForCp  please visit the project web-page or email Dr. Jun Pei at Jun_Pei@163.com.



