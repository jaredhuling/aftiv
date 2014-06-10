

genParamIndices <- function(res, X = NULL, Y = NULL, IS = NULL, sample.size = NULL) {
  #given a value for 2 of (X, Y, IS), this function returns indices
  #for parameters that meet the conditions specified for the two
  #given values
  param.names.actual <- c("Conf.Corr.X", "Conf.Corr.Y", "Instrument.Strength", "sample.size")
  
  stopifnot(sum(c(is.null(X), is.null(Y), is.null(IS), is.null(sample.size))) == 1)
  stopifnot(all(length(X) < 2, length(Y) < 2, length(IS) < 2))
  X.range <- unique(res[1, , 1])
  Y.range <- unique(res[1, , 2])
  IS.range <- unique(res[1, , 3])
  #n.data.range <- unique(res[1, , 4])
  sample.size.range <- unique(res[1, , 5])
  #lambda.range <- unique(res[1, , 6])
  
  if (!is.null(X)){if (!is.element(round(X, 4), round(X.range, 4))){stop(paste("X must be in:", X.range))}}
  if (!is.null(Y)){if (!is.element(round(Y, 4), round(Y.range, 4))){stop(paste("Y must be in:", Y.range))}}
  if (!is.null(IS)){if (!is.element(round(IS, 4), round(IS.range, 4))){stop(paste("IS must be in:", IS.range))}}
  if (!is.null(sample.size)){if (!is.element(round(sample.size, 4), round(sample.size.range, 4))){ 
    stop(paste("sample.size must be in:", sample.size.range))} 
  }
  nonnull.idx <- which(unlist(lapply(c("X", "Y", "IS", "sample.size"), function(x) !is.null(eval(parse(text = x))))))
  nonnull.idx <- match(param.names.actual[nonnull.idx], dimnames(res)[[3]])
  if (length(nonnull.idx) > 3) {stop("Only supply 3 variables please")}
  values <- unlist(lapply(c("X", "Y", "IS", "sample.size"), function(x) eval(parse(text = x))))
  
  idx <- which(round(res[1, , nonnull.idx[1]], 4) == round(values[1], 4) & 
                 round(res[1, , nonnull.idx[2]], 4) == round(values[2], 4) & 
                 round(res[1, , nonnull.idx[3]], 4) == round(values[3], 4))
  vary <- which(!(1:4 %in% nonnull.idx))
  list(idx = idx, vary = vary, held = list(idx = nonnull.idx, value = values))
}

meltResults <- function(res, idx.list, ss.static = F) {
  #function takes results array and 
  #idx given by 'genParamIndices' and
  #returns appropriate melted data
  param.names <- c("Conf.Effect.X", "Conf.Effect.Y", "Instrument.Strength", "sample.size")
  param.names.actual <- c("Conf.Corr.X", "Conf.Corr.Y", "Instrument.Strength", "sample.size")
  if (ss.static) {param.names[4] <- "placeholder"; param.names.actual[4] <- "placeholder"}
  vary.param.name <- param.names.actual[idx.list$vary]; vary.param.name <- vary.param.name[which(vary.param.name != "placeholder")]
  held.names <- param.names.actual[idx.list$held$idx]
  vary.loc <- match(vary.param.name, dimnames(res)[[3]]); vary.loc <- vary.loc[which(!is.na(vary.loc))]
  held.loc <- match(held.names, dimnames(res)[[3]])
  melted <- data.frame(array(0, dim = c(length(idx.list$idx) * dim(res)[1], 14)))
  for (i in 1:dim(res)[1]){
    est.name <- dimnames(res)[[1]][i]
    if (length(idx.list$idx) == 1) {
      replacement <- c(res[i, idx.list$idx, vary.loc],
                       rep(est.name, length(idx.list$idx)),
                       res[i, idx.list$idx, (length(param.names) + 3):18])
    } else {
      replacement <- cbind(res[i, idx.list$idx, vary.loc],
                           rep(est.name, length(idx.list$idx)),
                           res[i, idx.list$idx, (length(param.names) + 3):18])
    }
    melted[(1 + length(idx.list$idx) * (i - 1)):(length(idx.list$idx) * i),] <- replacement
  }
  for (i in c(1, 3:14)){melted[,i] <- as.numeric(melted[,i])}
  colnames(melted) <- c(vary.param.name, "variable", dimnames(res)[[3]][(length(param.names) + 3):18])
  melted
}



internalPlotBetaAndCI <- function(melted, truth, faceted = F, held, multi = F, CI) {
  cbPalette <- c("springgreen4", "blue", "darkorange", "red")
  var.idx <- ifelse(multi, 2, 1)
  footMessages <- c("Cor(U,Y):", "Cor(U,X):", "Cor(Z,X):")
  #values <- c(round(mean(unique(melted[,2 + held$idx[1]])), 3), round(mean(unique(melted[,2 + held$idx[2]])), 3))
  correlations <- c(round(mean(unique(melted$Cor.U.Y)), 3), round(mean(unique(melted$Cor.U.X)), 3), round(mean(unique(melted$Cor.Z.X)), 3))
  subtitle <- paste(paste(footMessages, correlations), collapse=", ")
  if (colnames(melted)[var.idx] == "Conf.Effect.X" | colnames(melted)[var.idx] == "Conf.Corr.X") {
    xlabtext <- "Confounder Effect Size on X"
  } else if (colnames(melted)[var.idx] == "Conf.Effect.Y" | colnames(melted)[var.idx] == "Conf.Corr.Y") {
    xlabtext <- "Confounder Effect Size on Y"
  } else if (colnames(melted)[var.idx] == "Instrument.Strength") {
    xlabtext <- "Instrument Strength"
  } else if (colnames(melted)[var.idx] == "sample.size") {
    xlabtext <- "Sample Size"
  } else {stop("Variable name is messed up. Check that out")}
  if (CI == "quantile90") {
    y.value.name <- colnames(melted)[which(colnames(melted) == "Med")]
    upper.int.name <- colnames(melted)[which(colnames(melted) == "Q95")]
    lower.int.name <- colnames(melted)[which(colnames(melted) == "Q5")]
  } else if (CI == "quantile95") {
    y.value.name <- colnames(melted)[which(colnames(melted) == "Med")]
    upper.int.name <- colnames(melted)[which(colnames(melted) == "Q97.5")]
    lower.int.name <- colnames(melted)[which(colnames(melted) == "Q2.5")]
  } else if (CI == "ci") {
    y.value.name <- colnames(melted)[which(colnames(melted) == "Mean")]
    upper.int.name <- colnames(melted)[which(colnames(melted) == "UCI")]
    lower.int.name <- colnames(melted)[which(colnames(melted) == "LCI")]
  } else {stop("wrong CI name")}
  
  p <- ggplot(data=melted, aes_string(x=colnames(melted)[var.idx], y=y.value.name, colour = "variable")) + 
    geom_line(size = 2) + 
    scale_fill_hue(c=65, l=80) + scale_colour_manual(values=cbPalette, name="Estimator  ") + ylab("Beta Estimate") + 
    xlab(xlabtext) + 
    ggtitle(as.expression(bquote(atop("Comparison of Estimators", atop(italic(.(subtitle))))))) + 
    theme_set(theme_complete_bw(24)) + 
    theme(plot.title = element_text(size = 26, face = "bold", colour = "black", vjust = -1),
          axis.title.x=element_text(size=22),
          axis.title.y=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22),
          legend.position="bottom") + 
    geom_hline(size = 2, yintercept = truth) +
    geom_ribbon(data = melted, aes_string(ymin = lower.int.name, ymax = upper.int.name), alpha = 0.3)
  if (faceted){
    if (multi) { 
      return(p + facet_grid(as.formula(paste(colnames(melted)[1], paste("~", colnames(melted)[1 + var.idx]))))
             )
    } else {
      return(p + facet_wrap(as.formula(paste("~", colnames(melted)[1+var.idx])), nrow = 1))
    }
  } else {
    if (multi) { 
      return(p + facet_grid(as.formula(paste(colnames(melted)[1], "~."))))
    } else {
      return (p)
    }
  }

}

plotBetaAndCI <- function(res, X = NULL, Y = NULL, IS = NULL, sample.size = NULL, truth,
                          faceted = F, CI = c("ci", "quantile90", "quantile95")){
  CI <- match.arg(CI)
  ss.loc <- match("sample.size", dimnames(res)[[3]])
  ss.static <- F
  if (length(unique(res[1,,ss.loc])) == 1) {sample.size <- unique(res[1,,ss.loc]); ss.static = T}
  nonnull.idx <- which(unlist(lapply(c("X", "Y", "IS", "sample.size"), function(x) !is.null(eval(parse(text = x))))))
  values <- lapply(c("X", "Y", "IS", "sample.size"), function(x) eval(parse(text = x)))
  vec.loc <- which(unlist(lapply(values, function(x) length(x) > 1)))
  if (length(vec.loc) > 1) {stop("Please provide no more than one vector")}
  if (length(vec.loc) == 0) { 
    idx <- genParamIndices(res, X, Y, IS, sample.size)
    melted <- meltResults(res, idx, ss.static = ss.static)
    internalPlotBetaAndCI(melted, truth, faceted = faceted, held = idx$held, multi = F, CI = CI)
  } else {
    for (i in 1:length(values[[vec.loc]])) {
      values.tmp <- values
      values.tmp[[vec.loc]] <- values[[vec.loc]][i]
      idx <- genParamIndices(res, X = values.tmp[[1]], Y = values.tmp[[2]], IS = values.tmp[[3]], sample.size = values.tmp[[4]])
      if (i == 1){
        melted <- meltResults(res, idx, ss.static = ss.static)
        melted <- cbind(rep(values[[vec.loc]][i], nrow(melted)), melted)
        colnames(melted)[1] <- c("X", "Y", "IS", "sample.size")[vec.loc]
      } else {
        melted.tmp <- meltResults(res, idx, ss.static = ss.static)
        melted.tmp <- cbind(rep(values[[vec.loc]][i], nrow(melted.tmp)), melted.tmp)
        colnames(melted.tmp)[1] <- c("X", "Y", "IS", "sample.size")[vec.loc]
        melted <- rbind(melted, melted.tmp)
      }
    }
    melted$variable <- as.factor(melted$variable)
    internalPlotBetaAndCI(melted, truth, faceted = faceted, held = idx$held, multi = T, CI = CI)
  }
  
}

meltParamVec <- function(res, X = NULL, Y = NULL, IS = NULL){
  idx <- genParamIndices(res, X, Y, IS)
  melted <- meltResults(res, idx)
  melted
}

plot.AFTsim <- function(res) {
  stopifnot(class(res) == "AFTsim")
  cbPalette <- c("springgreen4", "blue", "darkorange", "red")
  means <- vector()
  names <- vector()
  for (i in 1:length(res)) {
    means <- c(means, res[[i]])
    names <- c(names, rep(names(res)[i], length(res[[i]])))
  }
  dat2plot <- data.frame(variable = names, value = means)
  ggplot(dat2plot, aes(value, fill=variable, colour=variable)) + 
    geom_density(size = 1.25, alpha = 0.3) + 
    xlab("Beta Estimate") +
    scale_colour_manual(values=cbPalette, name="Estimator  ") +
    guides(fill=FALSE) + 
    scale_fill_manual(values=cbPalette) + theme_set(theme_complete_bw(24)) + 
    theme(plot.title = element_text(size = 26),
          axis.title.x=element_text(size=22),
          axis.title.y=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22),
          legend.position="bottom") + 
    geom_vline(size = 0.75, xintercept = attr(res, "truth"), colour = "grey50")
}

plotEstimatingEquations <- function(type = "AFT", sample.size, conf.corr.X = 0.0, conf.corr.Y = 0.0, instrument.strength, lambda, beta0, beta1) {
  # use this function to generate a plot of the values of an estimating equation (or vector of est eqns)
  # in order to diagnose potential issues
  # type can be set to a vector for multiple est. eqns (ie type = c("AFT", "AFT-IV"))
  types <- c("AFT", "AFT-IV", "AFT-2SLS", "AFT-IPCW", "AFT-2SLS-xhat")
  funcs <- c("vEvalAFTScore", "vEvalAFTivScore", "vEvalAFT2SLSScore", "vEvalAFTivIPCWScore", "vEvalAFT2SLSxhatScore")
  for (i in length(type)) {if (!is.element(type[i], types)) {stop("'type' must only contain 'AFT', 'AFT-IV',' AFT-2SLS' or 'AFT-IPCW'")}}
  Data.simu <- simIVData(sample.size, conf.corr.X = conf.corr.X, conf.corr.Y = conf.corr.Y, 
                         instrument.strength = instrument.strength, beta0 = beta0, beta1 = beta1, lambda = lambda)
  betas <- seq(beta1 - 2, beta1 + 2, by = 0.05)
  results <- data.frame(array(0, dim = c(0, 3)))
  for (i in 1:length(type)) {
    est.eqn <- match.fun(funcs[[match(type[i], types)]])
    est.eqn.values <- est.eqn(Data.simu, betas)
    results <- rbind(results, data.frame(betas, est.eqn.values, rep(type[i], length(betas))))
  }
  names(results) <- c("Beta", "Value", "Est.Eqn")
  results$Est.Eqn <- levels(results$Est.Eqn)[results$Est.Eqn]
  require(ggplot2)
  p <- ggplot(data = results, aes(x = Beta, y = Value, colour = Est.Eqn)) + 
    geom_line(size=2) + geom_vline(size=2, xintercept = beta1)
  plot(p)
  list(plot = p, data = results)
}

theme_complete_bw <- function(base_size = 12, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
      axis.ticks =        element_line(colour = "black"),
      axis.title.x =      element_text(size = base_size, vjust = 0.5),
      axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
      axis.ticks.length = theme_bw()$axis.ticks.length,
      axis.ticks.margin = theme_bw()$axis.ticks.margin,
      
      legend.background = element_rect(colour=NA), 
      legend.key =        element_rect(fill = NA, colour = "black", size = 0.25),
      legend.key.size =   theme_bw()$legend.key.size,
      legend.text =       element_text(size = base_size * 0.8),
      legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
      legend.position =   "right",
      
      panel.background = element_rect(fill = "white", colour = NA), 
      panel.border =     element_rect(fill = NA, colour = "grey50"), 
      panel.grid.major = element_line(colour = "grey90", size = 0.3), 
      panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
      panel.margin =     theme_bw()$panel.margin,
      
      strip.background = element_rect(fill = NA, colour = NA), 
      strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
      strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = -90),
      
      plot.background =  element_rect(colour = NA, fill = "white"),
      plot.title =       element_text(size = base_size * 1.2),
      plot.margin =      theme_bw()$plot.margin)
}




# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


