## The file stan_utilities.R contains the following functions:

########################################

## stan

## This is a replacement for the stan function in the rstan package that
##      caches already compiled models and runs in parallel if possible.
## Returns a stanfit object.

## stan(stanProgram, optimize = FALSE, cores = parallel::detectCores(),
##      pedantic = FALSE, sound = 1, ...)

# stanProgram:      Path to .stan file or program string.
#                       (Not required if model_code, file, or fit is present.)
# optimize:         Use sampling() (default) or optimizing().
# cores:            Number of cores to use if sampling().
# pedantic:         Whether to print status messages.
#                       (These messages are dumb. Don't use this.)
# sound:            Sound when sampling complete (if beepr is installed).
#                       (This is fun! Sounds like a kitchen timer.)
# ...               Other arguments passed to Stan:
#                       In this document, we will generally call stan by
#                       passing it a data list and a model code string.
#                       stan(data = data_model, model_code = stan_model)

########################################

## cleanup

## Cleans up the objects left over after fitting a Stan model:
##      data_model:     The data list passed to stan.
##      stan_model:     The model code as a character string.
##      fit_model:      The stanfit object output.
##      samples_model:  The samples as a data frame.
##      params_model:   The parameters to track for plotting.
## Returns nothing.

## cleanup(model)

# model:                Root name of the Stan model as character string.

########################################

## multiplot (winston Chang)

## Places ggplot objects in a grid.
## Returns nothing, but displays desired graphics.

## multiplot(..., plotlist=NULL, file, cols=1, layout=NULL)

# ...               One can pass ggplot objects directly.
#                           (But it's easier to pass them as a plotlist.)
# plotlist:         A list of ggplot objects.
# cols:             Number of columns in layout.
#                       If layout is not specified, multiplot will calculate
#                       how many rows to create using cols.
# layout:           A layout matrix.
#                       For example:
#                           matrix(c(1,2,3,3), nrow=2, byrow=TRUE)
#                       Plot 1 will go in the upper left,
#                       plot 2 will go in the upper right, and
#                       plot 3 will go all the way across the bottom.
#                       If present, cols is ignored.

########################################

## HDIofMCMC (John Kruschke)

## Computes highest density interval from a sample of representative values,
##      estimated as shortest credible interval.

## HDIofMCMC(sampleVec , cred_mass = 0.95)

# sampleVec:        A vector of representative values from a probability distribution.
# cred_mass:        A scalar between 0 and 1, indicating the mass within the credible
#                       interval that is to be estimated.

########################################

## HDI_calc

## Calculates HDI and prepares graphical parameters to pass to ggplot.
## Returns a list of three objects:
##      breaks:     The breaks to use for the histogram.
##      line_pos:   Parameters for plotting the HDI line segments.
##      text_pos:   Parameters for plotting the text annotations.

## HDI_calc(sample, cred_mass = 0.95, digits = 2)

# sample:           A vector of representative values from a probability distribution.
# cred_mass:        A scalar between 0 and 1, indicating the mass within the credible
#                       interval that is to be estimated.
# digits:           The number of decimal places when displaying the HDI.

########################################

## param_plot

## Plots a sampled parameter along with HDI.
## Returns a ggplot object.

## param_plot(sample, param, cred_mass = 0.95)

# sample:           A vector of representative values from a probability distribution.
# param:            A string with the name of the parameters to be plotted.
# cred_mass:        A scalar between 0 and 1, indicating the mass within the credible
#                       interval that is to be estimated.

########################################

## sample_plots

## Creates a grid of histograms of sampled parameters with HDIs.
## Returns nothing, but displays the result of multiplot.

## sample_plots(samples, params, cred_mass = 0.95, layout = NULL)

# samples:          A data frame extracted from a stanfit object.
# params:           A character vector of parameters to be plotted.
# cred_mass:        A scalar between 0 and 1, indicating the mass within the credible
#                       interval that is to be estimated.
# layout:           A layout matrix.
#                       For example:
#                           matrix(c(1,2,3,3), nrow=2, byrow=TRUE)
#                       Plot 1 will go in the upper left,
#                       plot 2 will go in the upper right, and
#                       plot 3 will go all the way across the bottom.
#                       If present, cols is ignored.

########################################

library(rstan)
library(parallel)
library(dplyr)
library(ggplot2)

stan <- function(stanProgram,      # path to .stan file or program string
                 optimize = FALSE, # use sampling() (default) or optimizing()
                 cores = parallel::detectCores(), # number of cores to use if sampling()
                 pedantic = FALSE, # whether to print status messages
                 sound = 1,        # sound when sampling complete (if beepr is installed)
                 ...) {            # further arguments to sampling() or optimizing()
    dots <- list(...)
    if (missing(stanProgram)) {
        if (!is.null(dots$fit)) {
            stanProgram <- rstan::get_stancode(dots$fit)
        }
        else if (!is.null(dots$file)) {
            stanProgram <- dots$file
            dots$file <- NULL
        }
        else {
            stanProgram <- dots$model_code
            dots$model_code <- NULL
        }
    }
    if (length(stanProgram) > 1 || !grepl("stan$", stanProgram)) { # program string
        tf <- tempfile()
        writeLines(stanProgram, con = tf)
        stanProgram <- file.path(dirname(tf), paste0(tools::md5sum(tf), ".stan"))
        if (!file.exists(stanProgram)) file.rename(from = tf, to = stanProgram)
    }
    else if (!file.exists(stanProgram)) stop(paste(stanProgram, "does not exist"))

    mtime <- file.info(stanProgram)$mtime
    chains <- dots$chains
    if (is.null(chains)) chains <- 4L
    stanProgram.rda <- gsub("stan$", "rda", stanProgram)
    if (!file.exists(stanProgram.rda) |
        file.info(stanProgram.rda)$mtime <  mtime |
        mtime < as.POSIXct(packageDescription("rstan")$Date) ) {

        if (pedantic) cat("Model needs compilation.\n")
        dots$chains <- 0L
        dots$file <- stanProgram
        stanExe <- suppressMessages(do.call(rstan::stan, args = dots))
        saveRDS(stanExe, file = stanProgram.rda)
        dots$file <- NULL
    }
    else {
        if (pedantic) cat("Loading cached model.\n")
        stanExe <- readRDS(stanProgram.rda)
    }

    if (optimize) {
        dots$object <- get_stanmodel(stanExe)
        dots$chains <- NULL
        out <- do.call(rstan::optimizing, args = dots)
        return(out)
    }

    # sampling()
    dots$chains <- 1L
    dots$fit <- stanExe
    if (chains == 0) return(stanExe)
    if (chains == 1) return(do.call(rstan::stan, args = dots))

    dots$refresh <- 500L
    sinkfile <- paste0(tempfile(), "_StanProgress.txt")
    cat("Refresh to see progress\n", file = sinkfile)
    #   if(interactive()) browseURL(sinkfile)
    callFun <- function(i) {
        dots$chain_id <- i
        sink(sinkfile, append = TRUE)
        on.exit(sink(NULL))
        return(do.call(rstan::stan, args = dots))
    }
    stopifnot(require(parallel))
    if (.Platform$OS.type != "windows") {
        out <- mclapply(1:chains, FUN = callFun, mc.cores = cores, mc.silent = TRUE)
    }
    else { # Windows
        cl <- makePSOCKcluster(cores)
        on.exit(stopCluster(cl))
        clusterExport(cl, envir = environment(), varlist = c("dots", "sinkfile"))
        clusterEvalQ(cl, expr = require(Rcpp))
        out <- parLapply(cl, X = 1:chains, fun = callFun)
    }
    if (!is.na(sound) && suppressWarnings(require(beepr, quietly = TRUE))) beep(sound)
    if (all(sapply(out, is, class2 = "stanfit"))) {
        out <- rstan::sflist2stanfit(out)
    }
    return(out)
}


cleanup <- function(model) {
    file1 <- paste("data", model, sep = "_")
    file2 <- paste("stan", model, sep = "_")
    file3 <- paste("fit", model, sep = "_")
    file4 <- paste("samples", model, sep = "_")
    file5 <- paste("params", model, sep = "_")
    cleanup_list <- c(file1, file2, file3, file4, file5)
    rm(list = cleanup_list, pos = ".GlobalEnv")
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

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

    if (numPlots == 1) {
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


HDIofMCMC = function(sampleVec , cred_mass = 0.95) {
    # Computes highest density interval from a sample of representative values,
    #   estimated as shortest credible interval.
    # Arguments:
    #   sampleVec
    #     is a vector of representative values from a probability distribution.
    #   cred_mass
    #     is a scalar between 0 and 1, indicating the mass within the credible
    #     interval that is to be estimated.
    # Value:
    #   HDIlim is a vector containing the limits of the HDI
    sortedPts = sort( sampleVec )
    ciIdxInc = floor( cred_mass * length( sortedPts ) )
    nCIs = length( sortedPts ) - ciIdxInc
    ciWidth = rep( 0 , nCIs )
    for (i in 1:nCIs) {
        ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
    }
    HDImin = sortedPts[ which.min( ciWidth ) ]
    HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
    HDIlim = c( HDImin , HDImax )
    return( HDIlim )
}


HDI_calc <- function(sample, cred_mass = 0.95, digits = 2) {
    # Get the breaks the way the hist function does.
    breaks <- pretty(range(sample),
                     n = nclass.Sturges(sample),
                     min.n = 1)
    # Get histogram output to calculate tallest bar
    h <- hist(sample, plot = FALSE)
    top <- max(h$counts)

    #Calculate endpoints of HDI
    HDI <- HDIofMCMC(sample, cred_mass)

    # Calculate position of line segments
    line_pos <- data_frame(x = c(HDI[1], HDI[1], HDI[2]),
                           xend = c(HDI[1], HDI[2], HDI[2]),
                           y = c(0.7*top, 0, 0.7*top),
                           yend = c(0, 0, 0),
                           size = factor(c(1, 2, 1)),
                           color = factor(c("blue", "black", "blue")),
                           linetype = factor(c("dashed", "solid", "dashed")))

    #Calculuate position of text annotations
    text_pos <- data_frame(x = c(HDI[1], HDI[2], (HDI[1] + HDI[2])/2),
                           y = c(0.75*top, 0.75*top, 0.05*top),
                           label = c(round(HDI[1], digits = digits),
                                     round(HDI[2], digits = digits),
                                     paste(100*cred_mass, "% HDI", sep = ""))) %>%
        mutate(width = strwidth(label, "inches") + 0.25,
               height = strheight(label, "inches") + 0.25)

    return(list(breaks = breaks, line_pos = line_pos, text_pos = text_pos))
}


param_plot <- function(sample, param, cred_mass = 0.95) {

    # Calculate HDI
    HDI <- HDI_calc(sample, cred_mass)
    df_sample <- data_frame(sample)

    # Plot
    ggplot(df_sample, aes(x = sample)) +
        geom_histogram(color = "gray", fill = "white", breaks = HDI$breaks) +
        geom_segment(data = HDI$line_pos, aes(x = x,
                                              xend = xend,
                                              y = y,
                                              yend = yend,
                                              size = size,
                                              color = color,
                                              linetype = linetype),
                     lineend = "round") +
        scale_size_manual(values = c(1, 2, 1)) +
        scale_color_manual(values = c("blue" = "blue",
                                      "black" = "black",
                                      "blue" = "blue")) +
        scale_linetype_manual(values = c("dashed" = "dashed",
                                         "solid" = "solid",
                                         "dashed" = "dashed")) +
        geom_text(data = HDI$text_pos, aes(x = x, y = y, label = label),
                  size = 7, vjust = 0) +
        ggtitle(param) +
        theme(plot.title = element_text(size = rel(2)),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = rel(1.5)),
              axis.title.y = element_blank(),
              legend.position = "none")
}



sample_plots <- function(samples, params, cred_mass = 0.95, layout = NULL) {
    # Subset samples by only the desired parameters
    samples <- samples[params]

    plots <- mapply(param_plot,
                    sample = samples,
                    param = params,
                    cred_mass = cred_mass,
                    SIMPLIFY = FALSE)

    multiplot(plotlist = plots, layout = layout)
}
