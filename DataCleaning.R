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
na_only <- data.frame(na_only, row.names = na_cols)
barplot(log(na_only[, "na_only"]), col = "darkblue", main = "NA histogram", 
     border = "white", xlab = "Samples", ylab = "NAs")
NESSRAdata <- NESSRAdata[ , ..not_na] # cleaning 

# About characheter or whatever 
char_count <- sapply(NESSRAdata, function(y) sum(length(which(is.character(y)))))
char_ind <- which(char_count != 0)
not_char <- colnames(NESSRAdata)[-na_ind]
NESSRAdata <- NESSRAdata[, ..not_char]

# Normalization or whatever you want to call it --------------------------------

# Blanzieri suggestion 
mean_vect <- as.numeric(mclapply(NESSRAdata, mean, mc.cores = 2))
hist(mean_vect) # he wants to cut out tail 
summary(mean_vect)

# cutting 250 mean(exp)
mean_ind <- which(mean_vect <= 25) # I think we could cut also less than 4...
NESSRAdata <- NESSRAdata[, ..mean_ind]

# new_mean_vect 
new_mean_vect <- as.numeric(mclapply(NE SSRAdata, mean, mc.cores = 2))
hist(new_mean_vect) # should adjust the cut over 50...



