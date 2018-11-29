# NESSRA input 
## be sure to setwd() to an empty new directory before running the script 

# run this function to download files ------------------------------------------

.downloadData <- function(platform,  cores = 1){
    
    plt <- getGEO(platform) # get platform data
    GSEids <- Meta(plt)$series_id # list of experiments
    plt_specs <- Table(plt) # platform specifications
    geneIds <- plt_specs$ID # platform genes 
    print(paste0("Total number of experiments: ", length(GSEids)))
    
    getExpr <- function(i){ 
        expr <- exprs(getGEO(i)[[1]])
        IDs <<- geneIds 
        msg <- F
        if (all(rownames(expr) == IDs)) {  
            fwrite(as.data.table(expr), file = as.character(i), row.names = T)
            msg <- T
        }
        msg
    }
    
    all_data <- mclapply(GSEids, getExpr, mc.cores = cores)
        
    print(cat("Total number of eset passed: ", 
                 sum(as.logical(all_data)), "on" , length(all_data)))
}

# run this function to concatenate matrices ------------------------------------

.concatMatrices <- function(dir, cores = 1) {

    files <- list.files(path = ".", all.files = F)
    
    fwrite(x = do.call(cbind, mclapply(files, function(x) fread(x)[,-1], 
                                      mc.cores = cores)), file = "output file")
}

# wrapper ----------------------------------------------------------------------

NESSRAquery <- function(platform, cores = 1) {
    
    require(data.table)
    require(parallel)
    require(GEOquery)
    
    print(paste0("Starting downloading data at ", Sys.time()))
    .downloadData(platform = platform, cores = cores)
    print("Concatenating many data, will require time and memory.")
    .concatMatrices(cores = cores)
    print(paste0("Done! At ", Sys.time()))
}










