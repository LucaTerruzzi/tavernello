require(data.table)

###### INPUT ######
# Name of the remote as defined in Rclone
remote <- 'DM_drive'
# Result ids mapping for vitis (download csv at http://gene.disi.unitn.it/test/gene_h.php)
result_ids <- fread('results_ids.csv')
# List of wanted genes 
wanted <- fread('example_list.csv')
###################

# Select only new dataset
mapping <- result_ids[result_ids$exp == 'vv_exprdata_2.csv',c(1,3)]

list <- mapping[mapping$lgn %in% wanted$V1,id]

getList <- function(i){
  system2('./rclone',c('copy',paste0(remote,':experiments_results/',i,'_Vv.expansion'),'Vv_results/expansion'))
  system2('./rclone',c('copy',paste0(remote,':experiments_results/',i,'_Vv.interactions'),'Vv_results/interactions'))  
}

lapply(list,getList)
