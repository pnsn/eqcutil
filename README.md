# eqcorrscan_utils
Utility classes and methods that support use of the waveform cross correlation package "EQcorrscan" for a range of projects at the [PNSN](https://pnsn.org).

# Authors  
This repository is abstracted from workflows written by N. Stevens and B. Johnson in the course of their work at the [PNSN](https://pnsn.org).

# License
![image](./docs/images/gplv3-with-text-136x68.png)  
This repository is distributed under the attached GNU General Public License v3

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
 * [obspy](https://obspy.org)
 * [eqcorrscan](https://eqcorrscan.readthedocs.io/en/latest/index.html)
 * [obsplus](https://niosh-mining.github.io/obsplus/versions/latest/index.html)
 * [pyrocko](https://pyrocko.org)  
 
We thank the development teams and contributor communities for their dedication to open-source scientific software!