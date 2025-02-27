# eqcutil
Flexible methods that support use of the waveform cross correlation package "EQcorrscan" for a range of projects at the [PNSN](https://pnsn.org).

a small, temporary, addition here

# Authors  
This repository is abstracted from workflows written by N. Stevens and B. Johnson in the course of their work at the [PNSN](https://pnsn.org) and related research endeavors.

# License
![image](./docs/images/gplv3-with-text-136x68.png)  
This repository is distributed under the attached GNU General Public License v3

# Installing with Conda

1) Clone this repository  
```
git clone https://github.com/pnsn/eqcutil.git
```
2) Use the provided `environment.yml` to create a `conda` environment  
```
conda env create -f environment.yml
```
3) Install `eqcutil` using `pip` backend from the root directory of this repo  
```
python -m pip install .
```

# Dependencies & Attribution
This repository build on the following open-source python projects:  
 * [eqcorrscan](https://eqcorrscan.readthedocs.io/en/latest/index.html)
 * [obspy](https://obspy.org)
 * [obsplus](https://niosh-mining.github.io/obsplus/versions/latest/index.html)
 * [pyrocko](https://pyrocko.org)  
 * [quakemigrate](https://quakemigrate.readthedocs.io/en/latest/index.html)
 
We thank the development teams and contributor communities of these projects for their dedication to open-source scientific software!