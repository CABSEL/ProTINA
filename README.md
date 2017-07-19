# SALMON

SALMON (**S**ystems **A**nalysis and **L**earning for inferring **Mo**difiers of **N**etworks) is a network-based analysis tool for identifying molecular targets (proteins) of drugs from (time-series) gene expression data as well as prior information of protein-gene network. SALMON provides a perturbation score of each protein for a drug treatment sample of interest. The positive (negative) sign of score indicates enhancement (attenuation) of protein activity under a given drug treatment, and the greater magnitude of the score, the greater perturbation was caused by the drug action. Based on the magnitudes of scores, target proteins can be ranked.

Both MATLAB and R versions of SALMON package are available. Please refer to the SALMON manuscript for more detailed information about SALMON algorithm.


### Prerequisites:
For the MATLAB package,
* MATLAB® R2014b or later
* [GLMNET package](http://web.stanford.edu/~hastie/glmnet_matlab/)

For the R package,
* R software (>=3.3.2)
* `devtools` package in R


### Last Update
* MATLAB package: ***salmon_1.0_MAT*** (19.07.2017)
* R package: ***salmon_1.0_R*** (19.07.2017)


### License
Redistribution and use in source and binary forms, with or without modification, are permitted provided agreeing to the Simplified BSD Style License.

[License](https://github.com/CABSEL/SALMON/blob/master/LICENSE) (RTF, 2 KB)

Read about Simplified [BSD Style License](http://opensource.org/licenses/bsd-license.php)

### Download and installation
Please refer to the [instruction for MATLAB](https://github.com/CABSEL/SALMON/tree/master/salmon_MATLAB/salmon_1.0_MAT/readme.md) and the [instruction for R](https://github.com/CABSEL/SALMON/tree/master/salmon_R/salmon_1.0_R/readme.md).

### Acknowledgements
This work has been supported by the ETH Zurich Research Grant.
