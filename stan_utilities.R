## The file stan_utilities.R contains the following functions:

########################################

## cleanup

## Cleans up the objects left over after fitting a Stan model:
##      data_model:     The data list passed to stan.
##      fit_model:      The stanfit object output.
##      samples_model:  The samples as a data frame.
##      params_model:   The parameters to track for plotting.
## Returns nothing.

## cleanup(model)

# model:                Root name of the Stan model as character string.

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

## HDIofMCMC (John Kruschke)

## Computes highest density interval from a sample of representative values,
##      estimated as shortest credible interval.

## HDIofMCMC(sampleVec , cred_mass = 0.95)

# sampleVec:        A vector of representative values from a probability distribution.
# cred_mass:        A scalar between 0 and 1, indicating the mass within the credible
#                       interval that is to be estimated.

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

## params_per_program

## Takes a list of sampled parameters for a collection of programs and refactors
## it so that each row is a program and the columns contain summaries of the
## parameters that correspond to that program.
## Returns a data frame.

## params_per_program(df)

# df:               A data frame of MCMC samples. For this to work, program-level
#                   parameters should be of the form 'param[#]' where # is the
#                   corresponding program new_ID.

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


library(dplyr)
library(ggplot2)


##### cleanup #####

cleanup <- function(model) {
    file1 <- paste("data", model, sep = "_")
    file2 <- paste("stan", model, sep = "_")
    file3 <- paste("fit", model, sep = "_")
    file4 <- paste("samples", model, sep = "_")
    file5 <- paste("params", model, sep = "_")
    cleanup_list <- c(file1, file2, file3, file4, file5)
    rm(list = cleanup_list, pos = ".GlobalEnv")
}


##### HDI_calc #####

HDI_calc <- function(sample, cred_mass = 0.95, digits = 2) {
    # Get the breaks the way the hist function does.
    breaks <- pretty(range(sample),
                     n = nclass.Sturges(sample),
                     min.n = 1)
    # Get histogram output to calculate tallest bar
    h <- hist(sample, plot = FALSE)
    top <- max(h$counts)


    #Get median
    sample_median <- median(sample)

    #Calculate endpoints of HDI
    HDI <- HDIofMCMC(sample, cred_mass)

    # Calculate position of line segments
    line_pos <- data_frame(x = c(HDI[1], HDI[1], HDI[2], sample_median),
                           xend = c(HDI[1], HDI[2], HDI[2], sample_median),
                           y = c(0.7*top, 0, 0.7*top, 0.9*top),
                           yend = c(0, 0, 0, 0),
                           size = factor(c(1, 2, 1, 1)),
                           color = factor(c("blue", "black", "blue", "red")),
                           linetype = factor(c("dashed", "solid", "dashed", "dashed")))

    #Calculuate position of text annotations
    text_pos <- data_frame(x = c(HDI[1],
                                 HDI[2],
                                 0.25*HDI[1] + 0.75*HDI[2],
                                 sample_median),
                           y = c(0.75*top, 0.75*top, 0.05*top, 0.91*top),
                           label = c(round(HDI[1], digits = digits),
                                     round(HDI[2], digits = digits),
                                     paste(100*cred_mass, "% HDI", sep = ""),
                                     paste("Median = ", round(sample_median,
                                                              digits = digits))))

    return(list(breaks = breaks, line_pos = line_pos, text_pos = text_pos))
}


##### HDIofMCMC #####

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


##### multiplot #####

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


##### params_per_program #####

params_per_program <- function(df) {
    df %>%
        summarise_each(funs(q025 = quantile(., probs = 0.025),
                            median = quantile(., probs = 0.5),
                            q975 = quantile(., probs = 0.975))) %>%
        lapply(unname) %>%      # Necessary because this is a list of named numbers
        as_data_frame() %>%
        gather(key, value) %>%
        separate(key, into = c("Param", "Program", "Stat"),
                 sep = "[\\[\\]]",
                 extra = "drop",
                 convert = TRUE) %>%    # Separate parameter from program from stat
        filter(Program != "") %>%       # Keep only program-specific values
        mutate(Program = factor(Program),
               Stat = str_sub(Stat, start = 2)) %>%     # Remove extra _
        unite(Parameter, Param, Stat) %>%
        spread(Parameter, value) %>%
        arrange(Program)
}


##### param_plot #####

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
        scale_size_manual(values = c(1, 2, 1, 1)) +
        scale_color_manual(values = c("blue" = "blue",
                                      "black" = "black",
                                      "blue" = "blue",
                                      "red" = "red")) +
        scale_linetype_manual(values = c("dashed" = "dashed",
                                         "solid" = "solid",
                                         "dashed" = "dashed",
                                         "dashed" = "dashed")) +
        geom_text(data = HDI$text_pos, aes(x = x, y = y, label = label),
                  size = 4, vjust = 0) +
        ggtitle(param) +
        theme(plot.title = element_text(size = rel(1.5)),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = rel(1.25)),
              axis.title.y = element_blank(),
              legend.position = "none")
}


##### sample_plots #####

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
