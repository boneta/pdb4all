# pdb4all

![python](https://img.shields.io/badge/python-3.7+-red.svg)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

*Convert between common protein pdb formats and names*

### Installation
Requires `python 3.7+`. No additional dependencies are needed.

*pdb4all* is frequently uploaded to [PyPI](https://pypi.org/project/pdb4all/). To install it, just execute the following command:
```
  pip install pdb4all
```

### Usage
```
  pdb4all [options] .pdb
```

#### Supported ForceFields name conventions
  * Amber
  * Charmm
  * OPLS

#### Supported formats (and software standards)
  * FASTA (output only)
  * generic (Maestro's Protein Preparation Wizard / pdb2pqr / ...)
  * GROMACS
  * fDynamo
