<img style = "float: right;" src = "https://github.com/CABSEL/SALMON/blob/master/salmon_logo.png" width="250" height="100" align="right"> 


# SALMON

SALMON (**S**ystems **A**nalysis and **L**earning for inferring **Mo**difiers of **N**etworks) is a network-based analytical tool for identifying protein targets of compounds from gene transcriptional profiles. SALMON scores proteins according to the enhancement or attenuation of the protein-gene regulatory interactions caused by a compound treatment. The magnitude of the scores indicates the degree of compound-induced alterations in the gene regulatory activity of a protein, while the sign of the scores indicates the direction of the alterations (positive: enhancement, negative: attenuation). In SALMON, the protein targets of a compound are ranked according to the magnitudes of the protein scores. 

MATLAB and R versions of SALMON are available. Please refer to SALMON manuscript for more detailed information.


### Prerequisites:
For MATLAB distribution,
* MATLAB® R2014b or later
* [GLMNET package](http://web.stanford.edu/~hastie/glmnet_matlab/)

For R distribution,
* R software (>=3.3.2)
* `devtools` package in R


### Last Update
* MATLAB package: **version 1.0** (19.07.2017)
* R package: **version 1.0** (19.07.2017)


### License
Redistribution and use in source and binary forms, with or without modification, are permitted provided agreeing to the *Simplified BSD Style License* (see [more](http://opensource.org/licenses/bsd-license.php)).

[License](https://github.com/CABSEL/SALMON/blob/master/LICENSE) (RTF, 2 KB)


### Download and installation
Please refer to [README for MATLAB](https://github.com/CABSEL/SALMON/blob/master/salmon_MATLAB/readme.md) and [README for R](https://github.com/CABSEL/SALMON/blob/master/salmon_R/readme.md).

### Acknowledgements
This work has been supported by the ETH Zurich Research Grant.
