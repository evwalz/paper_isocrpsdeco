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


The distributional predictive outcome of the individual method is either a normal distribution (*normal_deco.R*), an ensemble or ECDF (*ensemble_deco.R*) or a mixture distribution (*mixture_deco.R*). The MCS-DSC plots are computed with *figure_deco.R* and saved in [figures](https://github.com/evwalz/paper_isocrpsdeco/tree/main/ML_example/figures)