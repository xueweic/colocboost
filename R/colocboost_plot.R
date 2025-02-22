

cb_plot <- function(cb_output, y = "log10p", 
                    gene_name = NULL,
                    trait_idx = NULL, 
                    trait_names = NULL,
                    plot_cols = 2,
                    pos = NULL,
                    plot_cs_idx = NULL,
                    variant_coord = FALSE,
                    show_coloc = TRUE,
                    show_hits = FALSE,
                    show_cos_to_uncoloc = FALSE,
                    show_cos_to_uncoloc_idx = NULL,
                    show_cos_to_uncoloc_trait = NULL,
                    points_color = "grey80", cos_color = NULL,
                    ylim_each = TRUE, 
                    trait_legend_pos = "top",
                    trait_legend_size = 1.2,
                    cos_legend_pos = "bottomleft",
                    show_snp = FALSE,
                    lab_style = c(2, 1),
                    axis_style = c(2, 1),
                    title_style = c(2.5, 2), 
                    add_vertical = FALSE, add_vertical_idx = NULL, 
                    add_genetrack = FALSE,
                    ...){
  
    # get cb_plot_input data from colocboost results
    cb_plot_input <- get_input_plot(cb_output, plot_cs_idx = plot_cs_idx, 
                                    variant_coord = variant_coord,
                                    trait_names = trait_names,
                                    show_cos_to_uncoloc = show_cos_to_uncoloc,
                                    show_cos_to_uncoloc_idx = show_cos_to_uncoloc_idx,
                                    show_cos_to_uncoloc_trait = show_cos_to_uncoloc_trait)
    # get initial set up of plot
    cb_plot_init <- plot_initial(cb_plot_input, y = y, points_color = points_color, cos_color = cos_color,
                                 ylim_each = ylim_each, gene_name = gene_name,
                                 trait_legend_pos = trait_legend_pos, trait_legend_size = trait_legend_size,
                                 cos_legend_pos = cos_legend_pos,
                                 show_snp = show_snp, lab_style = lab_style, axis_style = axis_style,
                                 title_style = title_style, ... )
    
    colocboost_plot_basic(cb_plot_input, cb_plot_init, pos = pos,
                          trait_idx = trait_idx, plot_cols = plot_cols, 
                          add_vertical = add_vertical, add_vertical_idx = add_vertical_idx, 
                          show_hits = show_hits, 
                          ...)
  
}



# get input data for cb_plot
get_input_plot <- function(cb_output, plot_cs_idx = NULL, variant_coord = FALSE,
                           trait_names = NULL,
                           show_cos_to_uncoloc = FALSE,
                           show_cos_to_uncoloc_idx = NULL,
                           show_cos_to_uncoloc_trait = NULL){
  # redefined trait names
  if (!is.null(trait_names)){
      cb_output$data_info$traits_info$traits_names <- trait_names
      names(cb_output$data_info$z) <- trait_names
  }
  
  # extract results from colocboost
  analysis_trait <- cb_output$data_info$traits_info$traits_names
  target_idx <- which(cb_output$data_info$traits_info$is_target)
  if ( length(target_idx)!=0 ){
      target_trait <- analysis_trait[target_idx]
  } else { target_trait <- NULL }
  # extract z-scores
  snps <- cb_output$data_info$variants
  Z <- cb_output$data_info$z
    
  # if finemapping
  if (cb_output$data_info$n_trait==1){
      cb_output$cos_details$cos$cos_variants <- cb_output$ucos_details$ucos$ucos_variants
      cb_output$cos_details$cos$cos_index <- cb_output$ucos_details$ucos$ucos_index
      cb_output$cos_details$cos_top_variants <- cb_output$ucos_details$ucos_top_variants
      cb_output$cos_details$cos_traits <- cb_output$ucos_details$ucos_traits
      if (!is.null(cb_output$cos_details$cos$cos_variants)){
          cb_output$cos_details$cos_vcp <- cb_output$ucos_details$ucos_weight
      }
      
  }
  # extract coloc_cos
  coloc_snps <- cb_output$cos_details$cos$cos_variants
  coloc_cos <- cb_output$cos_details$cos$cos_index
  coloc_index <- cb_output$cos_details$cos_traits$trait_index
  # top_variants
  coloc_hits <- lapply(names(coloc_cos), function(cn){
    p <- grep(cn, rownames(cb_output$cos_details$cos_top_variants))
    cb_output$cos_details$cos_top_variants$top_index[p]
  })
  names(coloc_hits) <- names(coloc_cos)
  if (!is.null(coloc_cos)){
    # extract vcp each trait
    vcp <- lapply(1:length(analysis_trait), function(iy){
        pos <- which(sapply(coloc_index, function(idx) iy %in% idx))
        if (length(pos)!=0){
          w <- do.call(cbind, cb_output$cos_details$cos_vcp[pos])
          return(1-apply(1-w, 1, prod))
        } else {
          return(rep(0, length(snps)))
        }
    })
    ncos <- length(cb_output$cos_details$cos$cos_index)
    if (is.null(plot_cs_idx)){
      select_cs <- 1:ncos
    } else {
      if (length(setdiff(plot_cs_idx, c(1:ncos)))!=0){
        stop("please check plot_cs_idx!")
      }
      select_cs <- plot_cs_idx
    }
    coloc_snps <- coloc_snps[select_cs]
    coloc_cos <- coloc_cos[select_cs]
    coloc_index <- coloc_index[select_cs]
    coloc_hits <- coloc_hits[select_cs]
  } else {
    if (cb_output$data_info$n_trait==1){
      warnings("No fine-mapped causal effects in this region!")
    } else {
      warnings("No colocalized effects in this region!")
    }
    coloc_index <- NULL
    vcp <- rep(0, cb_output$data_info$n_variants)
  }
  # extract x axis
  if (variant_coord){
    snp.info <- do.call(rbind, lapply(snps, function(snp)strsplit(snp, ":")[[1]]))
    chrom <- as.numeric(sapply(snp.info[,1], function(ss) strsplit(ss, "chr")[[1]][2]))
    x <- data.frame("chrom" = chrom,
                    "pos" = as.numeric(snp.info[,2]))
    coloc_cos <- lapply(coloc_cos, function(cos) x$pos[cos])
    coloc_hits <- lapply(coloc_hits, function(cos) x$pos[cos])
  } else {
    x <- data.frame("chrom" = NA,
                    "pos" = c(1:length(snps)))
  }
  plot_input <- list("traits" = analysis_trait,
                     "target_trait" = target_trait,
                     "variants" = snps,
                     "x" = x,
                     "Zscores" = Z,
                     "vcp" = vcp,
                     "cos" = coloc_cos,
                     "cos_hits" = coloc_hits,
                     "coloc_index" = coloc_index)
  
  # check if plot cos to uncolocalized trait
  if (show_cos_to_uncoloc & !is.null(coloc_cos)){
      if (is.null(show_cos_to_uncoloc_idx)){
        cos_to_uncoloc <- coloc_cos
        cos_idx_to_uncoloc <- 1:length(coloc_index)
        if (is.null(show_cos_to_uncoloc_trait)){
          warning("Show all CoSs to uncolocalized traits.")
          trait_to_uncoloc <- lapply(coloc_index, function(cidx){ setdiff(1:length(analysis_trait), cidx) })
        } else {
          warning("Show all CoSs to uncolocalized traits ", paste(show_cos_to_uncoloc_trait, collapse = ","))
          trait_to_uncoloc <- lapply(coloc_index, function(cidx){ setdiff(show_cos_to_uncoloc_trait, cidx) })
        }
      } else {
        if (show_cos_to_uncoloc_idx > length(coloc_cos)){
          warning("There are only ", length(coloc_cos), " CoS in this region. ",
                  "Cannot show the ordered ", paste(show_cos_to_uncoloc_idx, collapse = ","), 
                  " CoS, which does not exist")
          trait_to_uncoloc <- cos_to_uncoloc <- cos_idx_to_uncoloc <- NULL
        } else {
          cos_idx_to_uncoloc <- show_cos_to_uncoloc_idx
          cos_to_uncoloc <- coloc_cos[show_cos_to_uncoloc_idx]
          if (is.null(show_cos_to_uncoloc_trait)){
            warning("Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","), 
                    " CoS for all uncolocalized traits.")
            trait_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx){
              l <- list(setdiff(1:length(analysis_trait), coloc_index[[idx]]))
              names(l) <- names(coloc_index[idx])
              return(l)
            })
          } else {
            warning("Show the ordered ", paste(cos_idx_to_uncoloc, collapse = ","), 
                    " CoS for traits ", paste(show_cos_to_uncoloc_trait, collapse = ","))
            trait_to_uncoloc <- sapply(show_cos_to_uncoloc_idx, function(idx){
              l <- list(setdiff(show_cos_to_uncoloc_trait, coloc_index[[idx]]))
              names(l) <- names(coloc_index[idx])
              return(l)
            })
          }
        }
      }
      if (!is.null(trait_to_uncoloc)){
        cos_uncoloc_texts <- rep("Uncolocalized effect", length(trait_to_uncoloc))
      } else {
        cos_uncoloc_texts <- NULL
      }
      draw_uncoloc <- sapply(trait_to_uncoloc, length)!=0
      if (any(draw_uncoloc)){
        pos <- which(draw_uncoloc)
        trait_to_uncoloc <- trait_to_uncoloc[pos]
        cos_to_uncoloc <- cos_to_uncoloc[pos]
        cos_idx_to_uncoloc <- cos_idx_to_uncoloc[pos]
        cos_uncoloc_texts <- cos_uncoloc_texts[pos]
      }
      uncoloc <- list("trait_to_uncoloc" = trait_to_uncoloc,
                      "cos_to_uncoloc" = cos_to_uncoloc,
                      "cos_idx_to_uncoloc" = cos_idx_to_uncoloc,
                      "cos_uncoloc_texts" = cos_uncoloc_texts)
      plot_input$uncoloc <- uncoloc
  }
  class(plot_input) <- "colocboost"
  return(plot_input)
}



plot_initial <- function(cb_plot_input, y = "log10p", 
                         points_color = "grey80", cos_color = NULL,
                         ylim_each = TRUE, gene_name = NULL,
                         trait_legend_size = 1.5,
                         trait_legend_pos = "right",
                         cos_legend_pos = "bottomleft",
                         show_snp = FALSE,
                         lab_style = c(2, 1),
                         axis_style = c(1, 1),
                         title_style = c(2.5, 2), 
                         ...){
  
  args = list(...)
  if (!exists("pch", args)) args$pch = 16
  
  # - set background point color and cos color pools
  args$bg <- points_color
  if (is.null(cos_color)){
    cos_color <- c("dodgerblue2", "#6A3D9A", "#FF7F00", "#FB9A99", "#33A02C",
                   "#A6CEE3",  "gold1", "#01665E","#FDBF6F", "#CAB2D6", "#B2DF8A", 
                   "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#01665E",
                   "#35978F", "#80CDC1", "#C7EAE5", "#003C30")
    # cos_color <- c("#1F70A9", "#33A02C", "#CAB2D6", "#EA7827")
  }
  args$col <- cos_color
  
  # - set data and x-lab and y-lab
  if (y == "log10p") {
    plot_data <- lapply(cb_plot_input$Zscores, function(z) {
      -log10(2 * pnorm(-abs(z)))
    })
    ylab = "-log10(p)"
  } else if (y == "z_original") {
    plot_data <- cb_plot_input$Zscores
    ylab = "Z score"
  } else if (y == "vcp"){
    plot_data <- cb_plot_input$vcp
    ylab = "VCP"
    if (length(cb_plot_input$traits)==1){ ylab = "PIP" }
    args$ylim <- c(0,1)
  } else {
    stop("Invalid y value! Choose from 'z' or 'z_original'")
  }
  if (!exists("xlab", args)) args$xlab = "variants"
  if (!exists("ylab", args)) args$ylab = ylab
  args$lab_size <- as.numeric(lab_style[1])
  args$lab_face <- lab_style[2]
  
  # - set title format
  args$title <- gene_name
  args$title_size <- as.numeric(title_style[1])
  args$title_face <- title_style[2]
  
  # - set x-axis and y-axis
  args$x <- cb_plot_input$x$pos
  args$y <- plot_data
  if (show_snp){
    args$xaxt = "n"
    args$xtext = cb_plot_input$variants
  } 
  args$axis_size <- as.numeric(axis_style[1])
  args$axis_face <- axis_style[2]
  
  # - set ylim for each subfigure
  if (exists("ylim", args)) ymax = rep(args$ylim[2],length(args$y)) else ymax = NULL
  if (ylim_each & is.null(ymax)) {
    ymax <- sapply(plot_data, function(p){
      valid_values <- p[is.finite(p)]
      if (length(valid_values) == 0) {
        ymax <- 330  # Default if no valid values
      } else {
        ymax <- max(valid_values)*1.05
      }
      return(ymax)
    })
  } 
  args$ymax = ymax
  
  # - set legend text position and format
  args$trait_legend_pos <- trait_legend_pos
  args$trait_legend_size <- trait_legend_size
  if (trait_legend_pos == "right"){
    args$trait_legend_angle = 90
  } else if (trait_legend_pos == "left") {
    args$trait_legend_angle = 270
  } else {
    args$trait_legend_angle = 0
  }
  
  if (!(cos_legend_pos %in% c("bottomright", "bottom", "bottomleft", "left",
                              "topleft", "top", "topright", "right", "center"))){
    cos_legend_pos = "bottomleft"
  }
  args$cos_legend_pos <- cos_legend_pos
  
  return(args)
}


colocboost_plot_basic = function (cb_plot_input, cb_plot_init,
                                  trait_idx = NULL, pos = NULL,
                                  plot_cols = 2, 
                                  add_vertical = FALSE, add_vertical_idx = NULL, 
                                  show_hits = TRUE, 
                                  ...) {
  
  args <- list(...)
  args <- c(args, cb_plot_init[c("xlab","ylab")])
  args$col = cb_plot_init$bg
  if (is.null(pos)){
    args$x = cb_plot_init$x
    y = cb_plot_init$y
  } else {
    args$x = cb_plot_init$x[pos]
    y = lapply(cb_plot_init$y, function(yy) yy[pos])
  }
  args$pch = cb_plot_init$pch
  args$cex.axis = cb_plot_init$axis_size
  args$cex.lab = cb_plot_init$lab_size
  args$font.lab = cb_plot_init$lab_face
  # - change position
  cb_plot_init$trait_legend_pos <- switch(cb_plot_init$trait_legend_pos,
                                          "right" = 4, "left" = 2, "top" = 3,"bottom" = 1)
  
  # - begin plotting
  coloc_cos <- cb_plot_input$cos
  traits <- cb_plot_input$traits
  if (is.null(trait_idx)){
    if (is.null(coloc_cos)){
      # - no colocalized effects, draw all traits in this region
      if (length(cb_plot_input$traits)==1){
        message("There is no fine-mapped causal effect in this region!. Showing margianl for this trait!")
      } else {
        message("There is no colocalization in this region!. Showing margianl for all traits!")
      }
      trait_idx <- 1:length(y)
    } else {
      n.coloc <- length(coloc_cos)
      coloc_index <- cb_plot_input$coloc_index
      trait_idx <- Reduce(union, coloc_index)
    }
    if (!is.null(cb_plot_input$target_trait)){
      p_target <- grep(cb_plot_input$target_trait, traits)
      include_target <- sapply(cb_plot_input$coloc_index, function(ci){ p_target %in% ci })
      if (any(include_target)){
          coloc_index <- cb_plot_input$coloc_index[order(include_target == "FALSE")]
          coloc_index <- Reduce(union, coloc_index)
          trait_idx <- c(p_target, setdiff(coloc_index, p_target))
      }
    }
  }
  if(length(trait_idx)==1){plot_cols=1}
  nrow <- ceiling( length(trait_idx) / plot_cols )
  if (!is.null(cb_plot_init$xtext)){bottom = 6} else {bottom=2}
  if (!is.null(cb_plot_init$title)){
    par(mfrow=c(nrow, plot_cols), mar=c(bottom,5,2,1), oma = c(0, 0, 3, 0))
  } else {
    par(mfrow=c(nrow, plot_cols), mar=c(bottom,5,2,1), oma = c(0, 0, 1, 0))
  }
  for (iy in trait_idx){
    args$y = y[[iy]]
    args$ylim = c(0, cb_plot_init$ymax[iy])
    if (!is.null(cb_plot_init$xtext)){
      args$xaxt = "n"
      do.call(plot, args)
      axis(1, at = args$x, labels = FALSE)
      text(x = args$x, y = par("usr")[3] - 0.1, labels = cb_plot_init$xtext, srt = 45, adj = 1, xpd = TRUE)
    } else {
      do.call(plot, args)
    }
    mtext(traits[iy], side = cb_plot_init$trait_legend_pos, line = 0.2, adj = 0.5,
          cex = cb_plot_init$trait_legend_size, font = 1)
    if (add_vertical){
      for (iii in 1:length(add_vertical_idx)){
        abline(v = add_vertical_idx[iii], col = '#E31A1C', lwd = 1.5, lty = 'dashed')
      }
    }
    
    # mark variants in CoS to colocalized traits
    if (!is.null(coloc_cos)){
      n.coloc <- length(coloc_cos)
      coloc_index <- cb_plot_input$coloc_index
      legend_text = list(col = vector())
      legend_text$col = head(cb_plot_init$col, n.coloc)
      
      # check which coloc set for this trait
      p.coloc <- sapply(coloc_index, function(idx) sum(idx==iy)!=0 )
      p.coloc <- which(p.coloc)
      coloc_cos.idx <- coloc_cos[p.coloc]
      for (i.cs in p.coloc){
        # add the points with specific color
        cs <- as.numeric(coloc_cos[[i.cs]])
        x0 <- intersect(args$x, cs)
        y1 = args$y[match(x0, args$x)]
        points(x0,y1, pch = 21, bg = legend_text$col[i.cs], col = NA, cex = 1.5,lwd = 2.5)
        if (show_hits){
          # add the hits points with "red"
          cs_hits <- as.numeric(cb_plot_input$cos_hits[[i.cs]])
          x_hits <- intersect(args$x, cs_hits)
          y_hits = args$y[match(x_hits, args$x)]
          points(x_hits, y_hits, pch = 21, bg=legend_text$col[i.cs], col = "#E31A1C", cex = 2, lwd = 3)
        }
      }
      
      # mark variants in CoS to uncolocalized traits
      uncoloc <- cb_plot_input$uncoloc
      if (!is.null(uncoloc)){
        p.uncoloc <- sapply(uncoloc$trait_to_uncoloc, function(idx) sum(idx==iy)!=0 )
        p.uncoloc <- which(p.uncoloc)
        texts <- shape_col <- texts_col <- c()
        for (i.uncoloc in p.uncoloc){
          uncoloc_trait <- uncoloc$trait_to_uncoloc[[i.uncoloc]]
          if (iy %in% uncoloc_trait){
            # add the points with specific color
            cs <- as.numeric(uncoloc$cos_to_uncoloc[[i.uncoloc]])
            x0 <- intersect(args$x, cs)
            y1 = args$y[match(x0, args$x)]
            points(x0,y1, pch = 4, col = adjustcolor(legend_text$col[i.uncoloc], alpha.f = 0.3), 
                   cex = 1.5, lwd = 1.5)
            texts <- c(texts, uncoloc$cos_uncoloc_texts[i.cs])
            shape_col <- c(shape_col, adjustcolor(legend_text$col[i.uncoloc], alpha.f = 1))
            texts_col <- c(texts_col, adjustcolor(legend_text$col[i.uncoloc], alpha.f = 0.8))
          }
        }
        if (length(texts)==0){next}
        legend(cb_plot_init$cos_legend_pos, texts, bty = "n", col = shape_col, text.col = texts_col,
               cex = 1.5, pt.cex = 1.5, pch = 4, x.intersp = 0.1, y.intersp = 0.3)
      }
    }
  }
  if (!is.null(cb_plot_init$title)){
    mtext(cb_plot_init$title, side = 3, line = 0, outer = TRUE, 
          cex = cb_plot_init$title_size, font = cb_plot_init$title_face)
  }
  return(invisible())
  
}
