library("GEOquery")
require("Biobase")
library("parallel")

# .inputCheck ------------------------------------------------------------------

.inputCheck <- function(platform, organism, expSubset){
    if (!is.character(platform)) {
        stop("platform is not character.")
    } 
    
    if(length(platform) != 1){
        stop("platform has not length 1.")
    }
    
    if (!is.null(organism)){
        if (!is.character(organism)) {
            stop("Organism is not character")
        }
        if (length(organism) != 1){
            stop("Organism has not length 1")
        }
    }
    
    if (!is.numeric(expSubset)){
        stop("expSubset is not numeric")
    }
    
}

# Get GPLdata ------------------------------------------------------------------

.getGPL <- function(platform, organism) {
    plt <- getGEO(platform)
    plt_specs <- Table(plt) # platform specifications
    
    if (!is.null(organism)) {
        if (!organism%in%plt_specs$`Species Scientific Name`) {
            stop("No such organism for this platform.")
        }
        plt_specs <- as.list(plt_specs[
            plt_specs$`Species Scientific Name` == organism,])
    } 
    list(plt, plt_specs)
}

# Constructor ------------------------------------------------------------------

.constructor <- function(...){
    
    nessraInput <- setClass("nessraInput", 
                            slots = list(
                                exprMatrix = "matrix", 
                                platformInfo = "list"
                                #samplePheno = "factor",
                                #batchVect = "factor",
                                #experimentVect = "integer"
                            ))
    
    out <- nessraInput(exprMatrix = ..1,
                        platformInfo = ..2)
    
    out
}

# Main Funciton ----------------------------------------------------------------
  
NESSRAquery <- function(platform, organism = NULL, cores = 1, expSubset = 5) {
    
    .inputCheck(platform, organism, expSubset)
    
    plts <- .getGPL(platform, organism) # get platform data
    plt <- plts[[1]]
    plt_specs <- plts[[2]]
    
    # creating Matrix ----------------------------------------------------------
    gseids <- Meta(plt)$series_id
    
    if (length(gseids > expSubset)) {
        gseids <- sample(gseids, expSubset)
        print("Random Sampled Experiments:")
        print(gseids)
    }
    
    final_exprs <- matrix(NA, length(plt_specs$ID), 1) 
    geneIds <- rownames(final_exprs) <- plt_specs$ID
    #phenoVect <- 1 
    #batchVect <- factor(0)
    #experimentVect <- integer(0)
    
    # getting exprsMatrices ----------------------------------------------------
    
    getExpr <- function(i){ #shoud add normalization from here
        expr <- exprs(getGEO(i)[[1]])
        IDs <<- geneIds 
        if (all(rownames(expr) == IDs)) {  
            ex_out <- expr 
        } else {
            ex_out <- matrix(NA, length(IDs))
        }
        ex_out
    }
    
    data_list <- mclapply(gseids, getExpr, mc.cores = cores)
    final_exprs <- as.matrix(lapply(data_list, cbind))
    final_exprs[(is.na(final_exprs))] <- NULL
    
    # Creating final object ----------------------------------------------------
    
    nInput <- .constructor(final_exprs, 
                           plt_specs)
    nInput
}

# finish to add other slots informations 


