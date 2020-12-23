diversity <- function(x, x_type = c("clonotypes", "counts"), method){
  #clonotype_vec
  if(x_type == "clonotypes"){
    x <- table(x)
  }else{
    stopifnot(is.numeric(x))
  }

  p <- x / sum(x)

  #Based on wikipedia page 
  #https://en.wikipedia.org/wiki/Diversity_index#Berger%E2%80%93Parker_index
  #Accessed on June 9, 2020

  switch(method,
         shannon = -sum(p * log(p)),
         inv_simpson = 1 / sum(p^2),
         simpson = sum(p^2),
         berger_parker = max(p),
         n_unique = length(p)
         )
}

diversity_all <- function(x, x_type = c("clonotypes", "counts")){
  #clonotype_vec
  if(x_type == "clonotypes"){
    x <- table(x)
  }else{
    stopifnot(is.numeric(x))
  }

  methods <- c("shannon", "simpson", "berger_parker", "n_unique")
  names(methods) <- methods
  sapply(methods, function(method){
        diversity(x, x_type = "counts", method = method)
  })

}
