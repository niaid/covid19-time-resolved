# DSB Normalization
# isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT" )

# Seurat v 2.3.4 e.g. : cell.columns.protein.matrix = object@assay$CITE@raw.data
# input matrix format for control and stained protein data.
#              AAACCTGCACGGCCAT_H1B1ln1 AAACCTGCATAGTAAG_H1B1ln1 AAACGGGAGTTAGCGG_H1B1ln1
#AnnexinV_PROT                       11                       20                       13
#BTLA_PROT                            1                        2                        1
#CD117_PROT                           .                        .                        1
#CD123_PROT                           2                        3                        2

# control.protein.matrix are the negative cells from cell hashing demultiplexing with multiseq deMULTIplex or HTODemux
# for example below, negative object would ha been object %>% SubsetData(ident.use = "Negative")
# control.protein.matrix = negative_cell_object@assay$CITE@raw.data


DSBNormalizeProtein = function(cell.columns.protein.matrix, control.protein.matrix, define.pseudocount = TRUE, pseudocount.use,  denoise_counts = TRUE, isotype.control.name.vec = NULL){

	adt = cell.columns.protein.matrix %>% as.matrix()
	adtu = control.protein.matrix %>% as.matrix()

	if(define.pseudocount == TRUE) {
		adtu_log = log(adtu + pseudocount.use)
		adt_log = log(adt + pseudocount.use)
	} else {
		# use +1 pseudocount for normalization
		adtu_log = log1p(adtu)
		adt_log = log1p(adt)
	}
	# apply scaling by background (subtract control means from counts and divide by control sd)
	mu_u = apply(adtu_log, 1 , mean)
	sd_u = apply(adtu_log, 1 , sd)
	norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)
	# transpose
	if(denoise_counts == FALSE){
	return(norm_adt) # add back to SCE or Seurat object .
	} else {
		suppressMessages(require(mclust))
		# fit mixture of 2 gaussians to each cell's protein data
		cellwise_background_mean = apply(norm_adt, 2, function(x) {
			g = Mclust(x, G=2, warn = F , verbose = F)
			return(g$parameters$mean[1])
		})

		# define pc1 through isotypes and background protein as a latent variable
		noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
		get_noise_vector = function(noise_matrix) {
			g = prcomp(t(noise_matrix), scale = TRUE)
			return(g$x[ ,1])
		}
		noise_vector = get_noise_vector(noise_matrix)

		# suppressMessages(library(limma))
		denoised_adt.list = list(denoised_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector), cellwise_background_mean = cellwise_background_mean)

	}
	return(denoised_adt.list)
}

# usage:
# normalize and add protein data back to Seurat object
# norm =  DSBNormalizeProtein(cell.columns.protein.matrix = object@assay$CITE@raw.data ... )
# h1 = SetAssayData(h1, new.data = norm, assay.type = "CITE", slot = "data")
# Full function call example
#adt_norm1_TEST = DSBNormalizeProtein(cell.columns.protein.matrix = object@assay$CITE@raw.data,
#                                     control.protein.matrix = negative_cell_object@assay$CITE@raw.data,
# 									  define.pseudocount = TRUE, pseudocount.use = 10, denoise_counts = TRUE,
#                                     isotype.control.name.vec =  c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",
#                                                                  "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT"))
