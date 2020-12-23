IN_DIR <- "data/CITE5p/all_batches/differential_expression_cell_freq/2020_07_24/expressionsets/"

files <- list.files(IN_DIR)

eset_list <- lapply(files, function(path){
  path <- paste0(IN_DIR, path)
  readRDS(path)
})

sapply(eset_list, colnames)

