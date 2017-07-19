# SALMON
SALMON (Systems Analysis and Learning for inferring Modifiers of Networks ) is a network-based analysis tool for identifying molecular targets (proteins) of drugs from (time-series) gene expression data as well as prior information of protein-gene network. SALMON provides a perturbation score of each protein for a drug treatment sample of interest. The positive (negative) sign of score indicates enhancement (attenuation) of protein activity under a given drug treatment, and the greater magnitude of the score, the greater perturbation was caused by the drug action. Based on the magnitudes of scores, target proteins can be ranked.

Both MATLAB and R versions of SALMON package are available. Please refer to the SALMON manuscript for more detailed information about DeltaNet algorithm.


### Prerequisites:
For MATLAB
DeltaNet: MATLAB® software (version R2014b or later) and SpaSM package from http://www2.imm.dtu.dk/projects/spasm (last version on 24.10.2012) and [GLMNET package](http://web.stanford.edu/~hastie/glmnet_matlab/).

DeltaNet-R:  R software (versions 3.2.2 or 3.2.3 ) and R packages: "lars" for DeltaNet-lar, "glmnet", "cvTools" and "methods" for DeltaNet-lasso, and "doParallel" and "foreach" for parallel computing.
Last Update
