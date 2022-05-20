# For the first time, you may need to run:
# install.packages("directlabels")
# install.packages("tidyverse")

library(directlabels)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("long_read_enhancement", "min_qual", "metric", "value"))
precision <- read_csv(args[2], col_names = c("long_read_enhancement", "TP", "FP"))

res <- res %>%
    filter(metric %in% c("TP", "FP")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(TP!=0 | FP!=0)

# Creating an empty column:
precision <- add_column(precision, min_qual = NA, .after="long_read_enhancement")

# Combine res & precision
total <- rbind(res, precision)
total$long_read_enhancement = factor(total$long_read_enhancement,
                                     levels=c('no_enhancement', 'SL2x1', 'SL2x2', 'SL2x3', 'L2x1', 'L2x2', 'L2x3', 'L2',
                                              '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'),
                                     labels=c('iGenVar 0.0.3: S',
                                              'iGenVar SL with 1x coverage', 'iGenVar SL with 2x coverage', 'iGenVar SL with 3x coverage',
                                              'iGenVar L with 1x coverage', 'iGenVar L with 2x coverage', 'iGenVar L with 3x coverage',
                                              'iGenVar 0.0.3: L',
                                              '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'))

scale_custom <- list(
      # https://stackoverflow.com/questions/46803260/assigning-40-shapes-or-more-in-scale-shape-manual
      scale_shape_manual(values = c(15, 16, 16, 16, 17, 17, 17, 18, 46, 46, 46, 46, 46, 46, 46, 46, 46)),
      # https://www.farb-tabelle.de/de/farbtabelle.htm
      scale_color_manual(breaks = c('iGenVar 0.0.3: S',
                                    'iGenVar SL with 1x coverage', 'iGenVar SL with 2x coverage', 'iGenVar SL with 3x coverage',
                                    'iGenVar L with 1x coverage', 'iGenVar L with 2x coverage', 'iGenVar L with 3x coverage',
                                    'iGenVar 0.0.3: L'),
                         values = c("chartreuse3",
                                    "deepskyblue", "dodgerblue1", "dodgerblue2",# "dodgerblue3", "dodgerblue4", "darkblue",
                                    "firebrick1", "firebrick1", "firebrick1", "firebrick3",
                                    "tan", "tan", "tan", "tan", "tan", "tan", "tan", "tan", "tan"),
                         na.value = "gray80"))

ggplot(total, aes(TP, FP, color = long_read_enhancement)) +
      geom_point(size=0.5) +
      scale_custom +
      # Add labels to curves
      # cex: character expansion; rot: rotation
      geom_dl(aes(label = long_read_enhancement), method = list("last.points", rot=-20, cex = 0.5, dl.trans(x=x+0.1))) +
      geom_path() +
      labs(x = "TP", y = "FP", color = "Tool",
            title="Long Read Enhancement - iGenVar 0.0.3 - Precision",
            subtitle="Short reads combined with low coverage (sampled) long reads") +
      lims(x=c(0,9641), y=c(0,9641)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) +
      theme(plot.title = element_text(size=20, hjust=0.5, face="bold", color="black")) +
      theme(plot.subtitle = element_text(size=10, hjust=0.5, face="italic", color="black"))

ggsave(args[3], width=10, height=8)
