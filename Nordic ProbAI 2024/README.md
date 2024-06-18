# Causality - Tutorial for Nordic Prob AI 2024

The tutorial at main.Rmd covers two topics:

### Part I. Effect identification from causal diagrams via:

-   Adjustment Criterion -- pcalg R package
-   ID Algorithm -- causaleffect R package

### Part II. Effect identification from PAGs via:

-   Generalized Adjustment Criterion (GAG) -- pcalg R package
-   CIDP Algorithm -- PAGId R package

## Setup

The following packages are required:

1.  pcalg and causaleffect R packages:

```         
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
install.packages(c("igraph", "pcalg", "dagitty", "causaleffect"), dependencies=TRUE)
```

2.  PAGId R package, available at <https://github.com/adele/PAGId>:

```         
install.packages("devtools")
devtools::install_github("adele/PAGId", dependencies=TRUE)
```
