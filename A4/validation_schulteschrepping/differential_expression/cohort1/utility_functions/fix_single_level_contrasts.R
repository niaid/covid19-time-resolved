#takes a formula and a data.frame that will be put into model.matrix and removes any terms in formula that have only one level in data
# returns a formula object 
update_formula_remove_single_level_contrasts <- function(object, data){
  my_terms <- terms(object)
  
  cols <- rownames(attr(my_terms, "factors"))
  
  rm_cols <- NULL
  for(col in cols){
    if(!col %in% colnames(data)){
      stop(paste("column:", col, "not in data!"))
    }
    unique_items <- unique(data[[col]])
    unique_items <- unique_items[!is.na(unique_items)]
    if(length(unique_items) < 2){
      rm_cols <- c(rm_cols, col)
      print(paste("removing terms involving variable '", col, "': Only one unique value:", unique_items))
    }
  }
  
  if(is.null(rm_cols)){
    return(object)
  }else{
    terms_mat <- attr(my_terms, "factors") 
    all_terms <- colnames(terms)

    terms_mat_rm <- terms_mat[rm_cols, , drop = FALSE]
    terms_mat_rm <- terms_mat[, apply(terms_mat_rm, 2, sum) > 0, drop = FALSE]

    rm_terms <- colnames(terms_mat_rm)

    keep_terms <- setdiff(colnames(terms_mat), rm_terms)
    print("Removing following terms from formula:")
    print(rm_terms)

    out_form <- paste(keep_terms, collapse = " + ")
    out_form <- paste0("~", out_form)
    out_form <- as.formula(out_form)
    return(out_form)
  }
  
}

#Usage
#meta <- mtcars[mtcars$cyl ==4, ]
#FORMULA <- ~cyl + mpg + mpg:cyl
#FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
#design <- model.matrix(FORMULA, meta)

# to test
#tmp <- update_formula_remove_single_level_contrasts(~cyl + cyl:mpg, mtcars)
#tmp <- update_formula_remove_single_level_contrasts(~sdfcyl + cyl:mpg, mtcars)
#tmp <- update_formula_remove_single_level_contrasts(~cyl + cyl:mpg + mpg, mtcars[mtcars$cyl ==4, ])
