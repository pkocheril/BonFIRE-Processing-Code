# BonFIRE Processing Code
Code for processing BonFIRE data. Continuously a work-in-progress. No guarantees of accuracy, completeness, or functionality.

Written in MATLAB (R2024b or earlier recommended for figure scaling).


## Highlights and main functions
* Data loading and identification by metadata
* Baseline correction
* Lifetime fitting with a convolved Gaussian IRF
* Batch processing over many files in subfolders
* 'Intelligent' fit type selection and data processing
* Streamlined curve fitting and data visualization


## Prerequisites
Several MATLAB toolboxes (tested in v24.2) are required (or at least recommended):
* Communications Toolbox
* Curve Fitting Toolbox
* DSP System Toolbox
* Image Processing Toolbox
* Optimization Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* Symbolic Math Toolbox
* System Identification Toolbox


## Organization
The main code ("bonfire_2d_proc_v#") should contain all dependencies and can be called in a folder organized as:
```
Folder/
├── bonfire_2d_proc_v#.m
└── Subfolder/
    ├── data.raw
    └── data.txt
```


## Usage
If the code is run in MATLAB with the folder structure as above and the default settings, a GUI should guide you through a "test-run" on your data.


# How to cite
If you found any of these functions useful, please consider citing one or more of the following papers:

1. PA Kocheril, H Wang, D Lee, N Naji, and L Wei. *J Phys Chem Lett* **2024**, *15* (19), 5306-5314.
2. PA Kocheril, J Du, H Wang, RE Leighton, D Lee, Z Yang, N Naji, A Colazo, and L Wei. *Chem Sci* **2025**, *16*, 14905-14918.

