

readSignatureFile = function(geneList) {
  load("shared_scratch/group1/msigdb/MSigDBSignatures.RData")
  signatureNames = grep("^REACTOME_", names(msigdbSignatures), value = TRUE)
  msigdbSignatures_mat = matrix(0, nrow = length(geneList), ncol = length(signatureNames),
                                dimnames = list(geneList, signatureNames))
  
  for(sig in signatureNames) {
    #printf("%s:\n", sig) #for debug
    #msigdb signatures are all undirected
    genesInSig = intersect(msigdbSignatures[[sig]], geneList)
    #the intersection is not empty
    if(length(genesInSig) > 0) msigdbSignatures_mat[genesInSig, sig] = 1
  }
  
  return(msigdbSignatures_mat)
}




#' Title
#'
#' M - size of the population
#' K - number of items with the desired characteristic in the population
#' N - number of the samples drawn
#' X - number of samples drawn with desired characteristic
hygerich <- function(population, desired, picked, successes, alternative="pos") {
  
  fold = successes / (picked * (desired / population))
  
  if(alternative == "pos") {
    #check whether picked samples are positively enriched within the special
    p = phyper(successes-1, desired, population-desired, picked, lower.tail = FALSE)
    
  } else if(alternative == "neg") {
    p = phyper(successes, desired, population-desired, picked, lower.tail = TRUE)
  } else {
    stop("unknown alternative")
  }
  
  return(list(p = p, fold = fold))
}


### Looks cool: mix colors for hexagons
### http://stats.stackexchange.com/questions/30788/whats-a-good-way-to-use-r-to-make-a-scatterplot-that-separates-the-data-by-trea


#' Title
#'
#' @param reduced_dim_expMat
#' @param anno
#' @param x_axis
#' @param y_axis
#'
#' @return
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
#' @examples
plotSignatures = function(reduced_dim_expMat, anno, x_axis, y_axis,
                          plot_type=c("regular", "hex"), summary_func = mean,
                          num_colors = 100, clip_ends = 0.01,
                          #default shape - all the same
                          shapes = as.factor(1), open_new_window = FALSE,
                          x_label, y_label, tit) {
  
  clip_ends = min(clip_ends, 1-clip_ends)
  plot_type = match.arg(plot_type)
  
  my_palette = colorRampPalette( brewer.pal( 11 , "RdBu" ) ) #11 is the max of RdBu
  my_colors = rev(my_palette(num_colors)) #rev (added 11/2) --> Chao asked me to make blue low and red high by default
  
  for(md_field in colnames(anno)) {
    if(open_new_window) windowIfNotInKnitr()
    
    if(plot_type == "regular") {
      g1 = ggplot(, aes(x = reduced_dim_expMat[, x_axis], y = reduced_dim_expMat[, y_axis])) +
        geom_point(size=2, aes(color = anno[,md_field], shape = shapes))
      
      scale_discrete_fun = scale_color_brewer
      scale_continuous_fun = scale_colour_gradientn
      
    } else if(plot_type == "hex") {
      g1 = ggplot(, aes(x = reduced_dim_expMat[, x_axis],
                        y = reduced_dim_expMat[, y_axis], z = anno[,md_field]))
      
      scale_discrete_fun = scale_fill_brewer
      scale_continuous_fun = scale_fill_gradientn
    }
    
    g1 = g1 + ggtitle(md_field)
    if(is.factor(anno[, md_field])) {
      #discrete color scale
      
      g1 = g1 + scale_discrete_fun(type = "qual", palette = "Paired")
      
    } else  {
      #continuous color scale
      
      # g1 = g1 + scale_color_distiller(type = "div", palette = "RdBu", direction = 1)
      
      #interpolate even more colors - not really distinguishable from just color brewer = color distiller
      g1 = g1 + scale_continuous_fun(colours = my_colors,
                                     values = c(0, seq(from = clip_ends, to = 1 - clip_ends, length.out = num_colors-2), 1))
      
    }
    
    if(length(unique(shapes)) > 1) {
      #user supplied shapes
      g1 = g1 + scale_shape_discrete(solid = TRUE)
    } else {
      #user didn't supply shapes (argument shapes kept its default value) or supplied one unchanging shape - drop the legend
      #I checked that this code has no effect when using hex signatures
      g1 = g1 + scale_shape_discrete(solid = TRUE, guide = FALSE)
    }
    
    if(plot_type == "hex") {
      g1 = g1 + stat_summary_hex(fun = summary_func)
      #can also do: g1+stat_summary_2d(fun = function(z) sum(z))
    }
    
    print(g1 + xlab(x_label) + ylab(y_label) + ggtitle(tit))
  }
}


#convenience function for backwards compatability
#' Title
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plotHexSignatures = function(...) {
  return( plotSignatures(plot_type = "hex", ...) )
}


# First implementation:

#'
#' #' Title
#' #'
#' #' @param reduced_dim_expMat
#' #' @param anno
#' #' @param x_axis
#' #' @param y_axis
#' #'
#' #' @return
#' #' @import ggplot2
#' #' @importFrom RColorBrewer colorRampPalette
#' #' @export
#' #'
#' #' @examples
#' plotSignatures = function(reduced_dim_expMat, anno, x_axis, y_axis) {
#'
#'   for(md_field in colnames(anno)) {
#'     windowIfNotInKnitr()
#'     g1<-ggplot(, aes(x = reduced_dim_expMat[, x_axis], y = reduced_dim_expMat[, y_axis])) +
#'       geom_point(size=2, aes(color = anno[,md_field])) +
#'       ggtitle(md_field)
#'     if(is.factor(anno[, md_field])) {
#'       #discrete color scale
#'
#'       #g1 = g1 + scale_color_discrete()
#'       g1 = g1 + scale_color_brewer(type = "qual", palette = "Paired")
#'
#'     } else  {
#'       #continuous color scale
#'
#'       #g1 = g1 + scale_colour_gradientn(colours=rainbow(2))
#'       #g1 = g1 + scale_colour_gradient2(low="blue", mid="white", high="red")
#'       #see also: http://stackoverflow.com/questions/16295440/r-ggplot-better-gradient-color for a way to interpolate colors between color brewer's limited values
#'
#'       #colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
#'       #g1 = g1 + scale_colour_gradientn(colours = colfunc(50))
#'
#'       #not good, returns discrete scale
#'       #g1 = g1 + scale_color_brewer(type = "div", palette = "RdBu")
#'       #this returns gradient based on the brewer palette:
#'       g1 = g1 + scale_color_distiller(type = "div", palette = "RdBu", direction = 1)
#'       #g1 = g1 + scale_color_distiller(type = "seq", palette = "YlGn", direction = 1)
#'
#'       #interpolate even more colors - not really distinguishable from just color brewer = color distiller
#'       mypal = colorRampPalette( brewer.pal( 11 , "RdBu" ) )
#'       g1 = g1 + scale_colour_gradientn(colours = mypal(100))
#'
#'     }
#'     print(g1)
#'   }
#' }
#'
#'
#' #' Title
#' #'
#' #' @param reduced_dim_expMat
#' #' @param anno
#' #' @param x_axis
#' #' @param y_axis
#' #' @param summary_func
#' #'
#' #' @return
#' #' @import ggplot2
#' #' @export
#' #'
#' #' @examples
#' plotHexSignatures = function(reduced_dim_expMat, anno, x_axis, y_axis,
#'                              summary_func = mean) {
#'
#'   for(md_field in colnames(anno)) {
#'     windowIfNotInKnitr()
#'     g1 = ggplot(, aes(x = reduced_dim_expMat[, x_axis],
#'                      y = reduced_dim_expMat[, y_axis], z = anno[,md_field])) +
#'       ggtitle(md_field)
#'
#'     if(is.factor(anno[, md_field])) {
#'       #discrete color scale
#'       g1 = g1 + scale_fill_brewer(type = "qual", palette = "Paired")
#'     } else  {
#'       #continuous color scale
#'       #g1 = g1 + scale_fill_distiller(type = "div", palette = "RdBu", direction = 1)
#'
#'
#'     }
#'
#'     g1 = g1 + stat_summary_hex(fun = summary_func)
#'     #can also do: g1+stat_summary_2d(fun = function(z) sum(z))
#'
#'     print(g1)
#'   }
#' }
#'
#'
#'
#' # #also useful: stat_density2d
#' # #(the z in the aes is the extra variable to be summarized)
#' # g1 = ggplot(, aes(x = res_pca$x[, "PC1"],
#' #                   y = res_pca$x[, "PC2"], z = anno[,md_field]) ) +
#' #   stat_density2d(aes(fill = ..level..), geom="polygon")
#' # g1 = g1 + scale_fill_distiller(type = "div", palette = "RdBu", direction = 1)
#' # print(g1)

