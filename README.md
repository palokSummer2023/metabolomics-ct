# Metabolomics analysis of serum and cecal metabolites

## Convexity:

## Cecal metabolite - cecal bacteria correlation

## Repo structure

| directory                  | contents                                         |
|----------------------------|--------------------------------------------------|
| /data                      | raw data files                                   |
| /generated/convexities.pdf | top convexity list                               |
| /generated                 | filtered data, correlation score and convexities |
| *.ipynb                    | notebook files with analysis and results         |
| *.py                       | helper libraries                                 |

## Running order

1. `dirs.ipynb` to make required directories,
2. `xonv.ipynb` to convert and consolidate data format,
3. `xonv_ext_meanstd_2reps.ipynb` to filter data,
4. `bac_mb_correlation-genus-experiment.ipynb`, `bac_mb_correlation-phylum-experiment.ipynb`, `convexity-analysis1.ipynb`, `convexity-dist-plots.ipynb` to do various analysis.
