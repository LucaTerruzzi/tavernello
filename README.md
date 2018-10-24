# tavernello

## NESSRAquery.R 

### Description 
This script contains a function that downloads every expressionSet coming from the *same microArray* platform, filters them according to the platform standards and concatenates them.

### Usage 
`NESSRAquery(platform, cores = 1)`

### Arguments 

```
platform        the ID of the platform to be considered. An example: "GPL198"
cores           number of cores to be used, few functions run in parallel.
```
### Note 
In the case of our project, the platform used was **"GPL198"** containing data only from *Arabidopsis thaliana*. \
The final output matrix can be used as input for [NESSRA](https://github.com/lucamasera/NESSRA) program. The function requires 
the followin packages: `parallel`and `data.table`availabe in CRAN, `GEOquery`available on Bioconductor. 
