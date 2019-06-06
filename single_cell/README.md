# droplasso single cell experiments: 
This repo contains R code for reproducing real data results and the mentioned figures, tables and selected genes in: 
> Khalfaoui, Beyrem, and Jean-Philippe Vert. "DropLasso: A robust variant of Lasso for single cell RNA-seq data." arXiv preprint arXiv:1802.09381 (2018).


## Data source: 
To evaluate the performance of methods for supervised classification, we collected 4 publicly available scRNA-seq datasets amenable to this setting. These were already processed as described by [Soneson, C. and Robinson, M. D. (2017)](https://www.nature.com/articles/nmeth.4612)  and as we obtained them from the conquer [website](http://imlspenticton.uzh.ch:3838/conquer/)

## Code description: 
- Data additional preparation, cross validation and comparaison results tables and figures generation are performed in main.R 
- Models are retrained with best found parameters on the entire datasets and stored in models.R
- Selected features lists are generated in models_list.R  and stored in model_selected folder.


## Notes

- Datasets are everytime uploaded and prepared and never stored for reproducibility.  
- Selected feature lists are uploaded in [DAVID Bioinformatics Database](https://david.ncifcrf.gov/) against the background gene lists (\_all lists) for functional analysis.

