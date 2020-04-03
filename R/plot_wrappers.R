#' PCA plot wrapper
#'
#' @param count_set A count_set of mRNA or miRNA counts generated using count_set
#' @param norm_method Name of the count_set assay to be plotted. Options are
#' "counts", "background_corrected","positive_control_scaled",
#' "housekeeping_scaled", "geNorm_housekeeping", "all_endogenous_scaled", "loess", "vsn", "quantile", "ruvIII".
#' Default = "housekeeping_scaled"
#'
#' @param title Plot title
#' @param comp1 PCA compenent to plot on the x axis. Default = Comp 1
#' @param comp2 PCA compenent to plot on the y axis Default = Comp 1
#' @param legend Option for plot legend. Default = TRUE
#' @param label Option to label points with sample id. Default= TRUE
#' @param colors Option to define colour scheme. Default = NA
#'
#' @return pca plots
#'
#' @importFrom pcaMethods prep pca R2cum scores
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot
#' @importFrom stats na.omit
#'
pca_plot_wrap <- function(count_set = NULL,
                          norm_method = "housekeeping_scaled",
                          title = NA,
                          comp1 = 1,
                          comp2 = 2,
                          legend = TRUE,
                          label = TRUE,
                          colors = NA) {

  get_counts <- paste0("assays(count_set)$", norm_method)
  data_in <- na.omit(eval(parse(text=get_counts)))

  points <- seq(15, (15+(length(unique(count_set$batch))-1)),
                length = length(unique(count_set$batch)))

  # data transformations
  md <- pcaMethods::prep(t(data_in), scale = "none", centre = FALSE)
  pc <- pca(md, method="svd", center=FALSE, nPcs=ncol(data_in))
  var_3 <- R2cum(pc)[3] # accumulated variance
  pc_1 <- round(pc@R2[comp1]*100, 2)
  pc_2 <- round(pc@R2[comp2]*100, 2)

  pc_scores <- as.data.frame(scores(pc))
  pc_scores <- data.frame(pc_scores, "group"=count_set$group, "samp_id"=count_set$samp_id)

  # initialise plot:
  graphics::plot(1, type="n", xlim=c(min(pc_scores[comp1])-5,
                                     max(pc_scores[comp1])+5),
                 ylim=c(min(pc_scores[comp2])-5,
                        max(pc_scores[comp2])+5),
                 axes=TRUE,
                 xlab=paste("PC", comp1, " - ", pc_1, "%", sep=""),
                 ylab=paste("PC", comp2, " - ", pc_2, "%", sep=""),
                 main=paste("PCA - ", title, sep="")
  )

  # add labels:
  graphics::abline(h = 0, v = 0, col = "gray", lty = 2)
  # add legend:
  if(legend==TRUE){
    legend("bottomright", c(as.character(unique(pc_scores$group))),
           col=colors[unique(count_set$group)],
           pch = c(rep(19, length(unique(pc_scores$group)))),
           title = "SAMPLE GROUPS", cex= 1)
    legend("bottomleft", c(as.character(unique(count_set$batch))),
           pch = points[unique(count_set$batch)],
           title = "BATCH", cex= 1)
  }
  if(label==TRUE){
    graphics::text(pc_scores[,comp1], pc_scores[,comp2], pc_scores$samp_id,
                   cex=1.1, pos=3, col="black")
  }
  graphics::points(pc_scores[,comp1], pc_scores[,comp2], cex = 2,
                   col = colors[count_set$group], pch= points[count_set$batch])
}

#' RLE plot wrapper
#'
#' @param count_set A count_set of mRNA or miRNA counts generated using count_set
#' @param norm_method Name of the count_set assay to be plotted. Options are
#' "counts", "background_corrected","positive_control_scaled",
#' "housekeeping_scaled","geNorm", "all_endogenous_scaled", "loess", "vsn", "quantile", "ruvIII".
#' Default = "housekeeping_scaled"
#' @param title Plot title
#' @param legend Option for plot legend. Default = TRUE
#' @param colors Option to define colour scheme. Default = NA
#'
#' @return rle plots
#'
#' @importFrom EDASeq plotRLE
#' @importFrom RColorBrewer brewer.pal

rle_plot_wrap <- function(count_set = NULL,
                          norm_method = "housekeeping_scaled",
                          title = NA,
                          legend = TRUE,
                          colors = NA) {

  get_counts <- paste0("assays(count_set)$", norm_method)
  data_in <- 2^(eval(parse(text=get_counts))) #EDASeq::plotRLE expects non log data

  EDASeq::plotRLE(data_in,
          main = paste("RLE - ", title, sep=""),
          xaxt = "n",
          las = 2,
          col = colors[count_set$group],
          outline = FALSE)

  axis(side = 1, at = 1:length(count_set$samp_id),
       labels = paste(count_set$samp_id, count_set$batch, sep = "_"), las = 2)

  # add legend:
  if(legend==TRUE){
    legend("bottomright", c(as.character(unique(count_set$group))),
           col=colors[unique(count_set$group)],
           pch = c(rep(19, length(unique(count_set$group)))),
           title = "SAMPLE GROUPS", cex=1)
  }

}

#' hclust plot wrapper
#'
#' @param count_set A count_set of mRNA or miRNA counts generated using count_set
#' @param norm_method Name of the count_set assay to be plotted. Options are
#' "counts", "background_corrected","positive_control_scaled",
#' "housekeeping_scaled", "geNorm", "all_endogenous_scaled", "loess", "vsn", "quantile", "ruvIII".
#' Default = "housekeeping_scaled"
#' @param title Plot title
#' @param legend Option for plot legend. Default = TRUE
#' @param colors Option to define colour scheme. Default = NA
#'
#' @return hclust plots
#'
#' @importFrom dplyr %>%
#' @importFrom dendextend colored_bars set %>%
#' @importFrom stats as.dendrogram na.omit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.omit

hclust_plot_wrap <- function(count_set = NULL,
                          norm_method = "housekeeping_scaled",
                          title = NA,
                          legend = TRUE,
                          colors = NA) {

  get_counts <- paste0("assays(count_set)$", norm_method)
  data_in <- na.omit(eval(parse(text=get_counts)))
  data_in <- t(data_in)

  # relabel hclust by samp_id
  rownames(data_in) <- paste(count_set$samp_id, count_set$batch, sep = "_")
  dd <- stats::dist(scale(data_in), method = "euclidean")
  hc <- stats::hclust(dd, method = "ward.D2")
  # labelling of nodes, by samples
  colors_hclust <- colors[count_set$group]

  # adding coloured nodes
  colours_order <- 1:(length(count_set$group))
  colours_tab <- as.data.frame(cbind(colours_order, colors_hclust))
  colors_hclust_ordered <- colours_tab[sort(colours_tab$colours_order)[hc$order],]

  # GENERATE DENDROGRAM
  dend <- data_in %>%
    scale %>%
    stats::dist(method = "euclidean") %>%
    stats::hclust(method = "ward.D2") %>%
    as.dendrogram
  dend <- dend %>%
    set("labels_col", value=as.character(colors_hclust_ordered$colors_hclust))
  plot(dend,
       main = paste("hclust - ", title, sep = ""),
       sub ="euclidean + ward")

  if(legend == TRUE){
    legend("topright",
          legend = as.character(unique(count_set$group)),
          col = colors[unique(count_set$group)],
          pch = c(rep(19, length(unique(count_set$group)))),
          title = "SAMPLE GROUPS",
          cex= 1)
  }

}

#' density plot wrapper
#'
#' @param count_set A count_set of mRNA or miRNA counts generated using count_set
#' @param norm_method Name of the count_set assay to be plotted. Options are
#' "counts", "background_corrected","positive_control_scaled",
#' "housekeeping_scaled", "geNorm", "all_endogenous_scaled", "loess", "vsn", "quantile", "ruvIII".
#' Default = "housekeeping_scaled"
#' @param title Plot title
#' @param legend Option for plot legend. Default = TRUE
#' @param colors Option to define colour scheme. Default = NA
#'
#' @return density plots
#'
#' @importFrom dplyr %>%
#' @importFrom dendextend colored_bars set %>%
#' @importFrom stats as.dendrogram na.omit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats na.omit

density_plot_wrap <- function(count_set = NULL,
                             norm_method = "housekeeping_scaled",
                             title = NA,
                             legend = TRUE,
                             colors = NA) {

  get_counts <- paste0("assays(count_set)$", norm_method)
  data_in <- na.omit(eval(parse(text=get_counts)))
  #data_in <- t(data_in)

  #group colour
  density_colours <- colors[count_set$group]

  # determine densities
  # everything:
  dens <- stats::density(data_in)

  all_dens <- lapply(seq_len(ncol(data_in)), function(i)
    stats::density(data_in[,i]))


  y_range <- range(sapply(seq_len(length(all_dens)), function(i)
    all_dens[[i]]$y))
  # to include everything in the range:
  y_range <- range(c(y_range, dens$y))
  x_range <- range(sapply(seq_len(length(all_dens)), function(i)
    all_dens[[i]]$x))

  # initialise plot
  graphics::plot(all_dens[[1]], xlim = x_range, ylim = y_range,
                 main=paste("Density - ", title, sep=""),
                 col = density_colours[1],
                 xlab="log2(counts)")
  # add samples
  sapply(2:length(all_dens), function(i)
    graphics::lines(all_dens[[i]], col = density_colours[i])
  )
  # add total density of all samples
  graphics::lines(dens, col ="black", lwd=2)
  # add legend
  if(legend==TRUE){
    legend("topright", c(as.character(unique(count_set$group)), "ALL SAMPLES"),
           col=c(colors[unique(count_set$group)], "black"),
           pch = c(rep(19, length(unique(count_set$group)))),
           title = "SAMPLE GROUPS", cex=1)
  }
}











