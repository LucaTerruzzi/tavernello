# DATA EXPLORATION

library(data.table)
library(GEOquery)
library(ggplot2)
library(parallel)

# Samples platform info --------------------------------------------------------

platform <- "GPL198"
plt <- getGEO(platform) # get platform data
GSEids <- Meta(plt)$series_id # list of experiments
GSMids <- Meta(plt)$sample_id
plt_specs <- Table(plt) # platform specifications
geneIds <- plt_specs$ID # platform genes 
print(paste0("Total number of experiments: ", length(GSEids)))

# Beloved data -----------------------------------------------------------------

NESSRAdata <- fread(input = "input_NESSRA.csv", sep = ",", header = T) 

# Data cleaning ----------------------------------------------------------------

# Cleaning from columns duplicates 
dups_bool <- duplicated(colnames(NESSRAdata)) # boolean of duplicated columns 
select_cols <- colnames(NESSRAdata)[dups_bool]
keep_cols <- colnames(NESSRAdata)[!dups_bool]
NESSRAdata <- NESSRAdata[ , ..keep_cols] # cleaning 

# About the NAs 
na_count <- sapply(NESSRAdata, function(y) sum(length(which(is.na(y)))))
na_ind <- which(na_count!=0)
na_cols <- colnames(NESSRAdata)[na_ind]
not_na <- colnames(NESSRAdata)[-na_ind]
NESSRAdata <- NESSRAdata[ , ..not_na] # cleaning 

# About characheter or whatever 
char_count <- sapply(NESSRAdata, function(y) sum(length(which(is.character(y)))))
char_ind <- which(char_count != 0)
not_char <- colnames(NESSRAdata)[-char_ind]
NESSRAdata <- NESSRAdata[, ..not_char]

# Normalization or whatever you want to call it --------------------------------

data_log <- log(NESSRAdata)

# suggestion 
mean_vect <- as.numeric(mclapply(NESSRAdata, mean, mc.cores = 2))
hist(log(mean_vect+1), col = "blue", border = "white") # he wants to cut out tail 
summary(mean_vect)

# cutting 250 mean(exp)
mean_ind <- which(mean_vect <= 25) # I think we could cut also less than 4...
NESSRAdata <- NESSRAdata[, ..mean_ind]

# new_mean_vect 
new_mean_vect <- as.numeric(mclapply(NESSRAdata, mean, mc.cores = 2))
hist(new_mean_vect) # should adjust the cut over 50...
boxplot(log(mean_vect + 1), log(new_mean_vect + 1))
boxplot(mean_vect)
boxplot(new_mean_vect)

fwrite(NESSRAdata, file = "./Norm_data")

# Check log tranform 

checkLog <- function(sample) {
    qx <- as.numeric(stats::quantile(sample,
                                     c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
    logTransformed <- !((qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) ||
                            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
    logTransformed 
}

logScaled <- apply(NESSRAdata, 2, checkLog)
sum(logScaled)

NESSRAdata_log <- NESSRAdata <- NESSRAdata[, ..logScaled]
mean_log_vect <- apply(NESSRAdata, 2, mean)
hist(mean_log_vect) # plot the mean log vect 



