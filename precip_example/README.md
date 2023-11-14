# Example: Precipitation 

CRPS decomposition for Case Study "Probabilistic quantitative precipitation forecasts" from [IDR paper](https://academic.oup.com/jrsssb/article/83/5/963/7056107?login=true).

At each of the 4 different meteorological stations:

- Brussels
- Frankfurt
- London
- Zurich

different probabilistic precipitation forecasts are issued:

- ECMWF Ensemble
- BMA
- EMOS
- HCLR
- IDR<sub>cw</sub>
- IDR<sub>sbg</sub>
- IDR<sub>icx</sub>

The CRPS decomposition of the distributional predictive forecasts is computed  with *precip_deco.R*. The MCS-DSC plots are computed with *figure_deco.R* and saved in [figures](https://github.com/evwalz/paper_isocrpsdeco/tree/main/precip_example/figures)
