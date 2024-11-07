# eqcorrscan_utils
Flexible methods that support use of the waveform cross correlation package "EQcorrscan" for a range of projects at the [PNSN](https://pnsn.org).

# Authors  
This repository is abstracted from workflows written by N. Stevens and B. Johnson in the course of their work at the [PNSN](https://pnsn.org) and related research endeavors.

# License
![image](./docs/images/gplv3-with-text-136x68.png)  
This repository is distributed under the attached GNU General Public License v3

# Structure  
The python methods in this repository are roughly grouped by their intended purpose
which can spread in a number of different directions depending on what the user wants
to accomplish. The current structure follows:

## eqcorrscan_utils  
 - augment -- methods that modify ObsPy `Catalog` and EQcorrscan `Template` objects  
 - client -- lightweight wrappers around the `EventBank` and `WaveBank` classes from `ObsPlus`  
 - io -- methods for converting other metadata formats into ObsPy `Catalog` objects  
 - process -- methods wrapping core `EQcorrscan` processes  
 - util -- utilities for the utilities! Helpful other tools (e.g., setting up a Logger instance)  
 - visualize -- methods to help visualize `EQcorrscan` objects and outputs  

# Installation with Conda/Pip from GitHub
Current installation (particularly on Apple Silicon) seems to require
a pre-installation of `scipy` and `eqcorrscan` to comple correctly

1) Create a `conda` environment  
```conda env create -n eqc_util```  
2) Activate environment  
```conda activate eqc_util``` 
3) Pre-install dependencies  
```conda install python=3.10 scipy git```  
4) Pre-install `eqcorrscan` per their `conda` instructions  
```conda install -c conda-forge eqcorrscan```  
5) Install `eqcorrscan_utils` from GitHub  
```pip install git+https://github.com/pnsn.eqcorrscan_utils.git```  
6) Install optional development environment(s), e.g.,  
```conda install ipython```  

# Dependencies & Attribution
This repository build on the following open-source python projects:  
 * [eqcorrscan](https://eqcorrscan.readthedocs.io/en/latest/index.html)
 * [obspy](https://obspy.org)
 * [obsplus](https://niosh-mining.github.io/obsplus/versions/latest/index.html)
 * [pyrocko](https://pyrocko.org)  
 * [quakemigrate](https://quakemigrate.readthedocs.io/en/latest/index.html)
 
We thank the development teams and contributor communities of these projects for their dedication to open-source scientific software!