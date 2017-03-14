space_time_plot <- function(the_brick,
                            out_file,
                            timelim,
                            timeaxis,
                            timelab = "Year AD",
                            zlim = NULL,
                            zlim_mid_range = NULL,
                            zlab,
                            zaxis,
                            zlim_colors = brewer.pal(9, "YlGn"),
                            fig_width = 5,
                            graph_height = 1.5,
                            margin = 0.1,
                            smoother = dnorm(seq(-10,10,1), sd=5)){
  
  mean.all <- mean(the_brick[], na.rm = T) %>% round()
  mean.spatial <- mean(the_brick, na.rm = T)
  mean.temporal <- cellStats(the_brick, mean, na.rm = T)
  mean.temporal.smooth <- stats::filter(mean.temporal, filter = smoother)
  
  aspect <- ifelse(raster::isLonLat(the_brick),
                   1/cos((mean(ymin(the_brick),ymax(the_brick)) * pi)/180),
                   nrow(the_brick)/ncol(the_brick))
  
  plot_width <- fig_width - (margin * 2)
  plot_height <- plot_width * aspect
  fig_height <- plot_height + (margin * 3) + graph_height
  
  if(is.null(zlim)){
    zlim <- c(min(the_brick[], na.rm = T),max(the_brick[], na.rm = T))
  }
  
  if(is.null(zlim_mid_range)){
    zlim_mid_range <- zlim
  } 
  
  colors.begin.range <- c(zlim[1],max(zlim[1],zlim_mid_range[1]-1))
  colors.end.range <- c(zlim_mid_range[2]+1,zlim[2])
  
  colors.begin <- rep(head(zlim_colors,1),
                      length(colors.begin.range[1]:colors.begin.range[2]))
  colors.mid <- colorRampPalette(zlim_colors)(length(zlim_mid_range[1]:zlim_mid_range[2]))
  colors.end <- rep(tail(zlim_colors,1),
                    length(colors.end.range[1]:colors.end.range[2]))
  colors <- c(colors.begin,
              colors.mid,
              colors.end)
  
  cairo_pdf(filename = out_file,
         width = fig_width,
         height = fig_height,
         antialias = "none", 
         bg = "white",
         pointsize = 8,
         fallback_resolution = 600)
  par(mai=c(graph_height + (margin * 2),
            margin,
            margin,
            margin),
      xpd=F)
  
  plot(1,
       type='n',
       xlab="",
       ylab="", 
       xlim=c(extent(the_brick)@xmin,extent(the_brick)@xmax),
       ylim=c(extent(the_brick)@ymin,extent(the_brick)@ymax), 
       xaxs="i",
       yaxs="i",
       axes=FALSE,
       main='')
  
  plot(mean.spatial,
       maxpixels=1000000,
       zlim=zlim,
       add=T,
       col=colors,
       useRaster=TRUE, legend=FALSE)
  
  par(mai=c(margin,
            margin,
            (margin * 2) + plot_height,
            margin), xpd=T, new=T)
  
  plot(1,
       type='n',
       xlab="",
       ylab="",
       xlim=c(0,fig_width),
       ylim=zlim_mid_range,
       xaxs="i",
       yaxs="i",
       axes=FALSE,
       main='')
  
  legend.breaks <- seq(from=zlim_mid_range[1], to=zlim_mid_range[2], length.out=(length(colors[zlim_mid_range[1]:zlim_mid_range[2]])+1))
  rect(col=colors[zlim_mid_range[1]:zlim_mid_range[2]],
       border=NA,
       ybottom=zlim_mid_range[1]:(zlim_mid_range[2]-1),
       ytop=(zlim_mid_range[1]+1):zlim_mid_range[2],
       xleft=0.15,
       xright=0.35,
       xpd=T)
  text(x = 0,
       y=mean(zlim_mid_range),
       labels=zlab,
       adj=c(0.5,1),
       cex=1,
       srt = 90,
       family='Helvetica Bold')
  text(x = 0.5, 
       y = c(mean.all,
             zlim_mid_range[1],
             zlim_mid_range[2],
             zaxis), 
       labels = c(mean.all,
                  ifelse(zlim_mid_range[1] > zlim[1],paste0("< ",zlim_mid_range[1]),zlim[1]),
                  ifelse(zlim_mid_range[2] < zlim[2],paste0("> ",zlim_mid_range[2]),zlim[2]),
                  zaxis),
       adj=c(0.5,0.5),
       cex=0.8,
       family='Helvetica Bold')
  
  par(mai=c(margin,
            margin * 8,
            (margin * 2) + plot_height,
            margin * 2),
      xpd=T,
      new=T)
  plot(1, type='n', xlab="", ylab="", xlim=timelim, ylim=zlim_mid_range, xaxs="i", yaxs="i", axes=FALSE, main='')
  abline(h = mean.all,col = "gray50", lty = 1, xpd = F)
  lines(y = mean.temporal, x = timelim[1]:timelim[2])
  lines(y = mean.temporal.smooth, x = timelim[1]:timelim[2], col = "red")
  axis(2,
       at = c(mean.all,
              zlim_mid_range[1],
              zlim_mid_range[2],
              zaxis),
       labels = F)
  
  par(mai=c(margin,
            margin * 8,
            (margin * 2) + plot_height,
            margin * 2),
      xpd=T,
      new=T)
  plot(1, type='n', xlab="", ylab="", xlim=timelim, ylim=c(0,graph_height), xaxs="i", yaxs="i", axes=FALSE, main='')
  segments(x0 = timeaxis, 
           x1 = timeaxis, 
           y0 = 0, 
           y1 = graph_height - (margin*2),
           col = "gray50",
           lty = 3)
  text(x = timeaxis,
       y = graph_height - (margin*1),
       labels = timeaxis,
       adj = c(0.5,0.5),
       cex=0.8,
       family='Helvetica Bold')
  text(x = mean(timelim),
       y = graph_height,
       labels = timelab,
       adj = c(0.5,0.5),
       cex=1,
       family='Helvetica Bold')
  
  dev.off()
  distill(out_file)
}