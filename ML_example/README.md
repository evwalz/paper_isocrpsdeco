# Example: ML Benchmark Datasets

CRPS decomposition for uncertainty quantification methods presented in Block 2 of Table 6 in [EasyUQ paper](https://arxiv.org/pdf/2212.08376.pdf) with corresponding [EasyUQ Repo](https://github.com/evwalz/easyuq). 

The considered methods are:

- Conformal prediction (CP)
- smooth CP
- EasyUQ 
- smooth EasyUQ
- Laplace
- MC Dropout
- Single gaussian


The distributional predictive outcomes of the individual method are either normal distributions (*normal_deco.R*), ensembles or ECDFs (*ensemble_deco.R*) or mixture distributions (*mixture_deco.R*):

- Normal: Laplace and Single Gaussian
- Ensemble / ECDF: CP / EasyUQ
- Mixture: smooth EasyUQ, smooth CP and MC Dropout


The MCS-DSC plots are computed with *figure_deco.R* and saved in [figures](https://github.com/evwalz/paper_isocrpsdeco/tree/main/ML_example/figures)