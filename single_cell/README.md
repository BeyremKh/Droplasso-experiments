## droplasso single cell experiments: 
This repo contains R code for reproducing real data results and the mentioned figures, tables and selected genes in: 
> Khalfaoui, Beyrem, and Jean-Philippe Vert. "DropLasso: A robust variant of Lasso for single cell RNA-seq data." arXiv preprint arXiv:1802.09381 (2018).


# dependencies: 
The expreriments mainly use our developed R package that implements the dropout lasso, as described in the above paper. 

## Installation
```{r}
library(devtools)
install_github("jpvert/droplasso")
```
## Other dependencies
Please install if not already installed
```r
install.packages(c("rmarkdown", "glmnet","ROCR"))
```

# Notes

- Some implementations depend on `parallel::mclapply` which may be inconvenient for MS Windows users. See [here](http://www.r-bloggers.com/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/) for an easy solution, or manually revert all `mclapply` to `lapply`.
- For Jupyter Notebooks, we make use of an R native kernel. See [here](https://github.com/IRkernel/IRkernel) for an installation. Otherwise, R markdowns are also available. 
