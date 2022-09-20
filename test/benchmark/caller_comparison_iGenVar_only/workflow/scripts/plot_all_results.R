# For the first time, you may need to run:
# install.packages("directlabels")
# install.packages("tidyverse")

library(directlabels)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

plotTitle <- args[1]
res <- read_tsv(args[2], col_names = c("input_combination", "min_qual", "metric", "value"))
f1 <- read_csv(args[3], col_names = c("input_combination", "precision", "recall"))

res <- res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0)

# Creating an empty column:
f1 <- add_column(f1, min_qual = NA, .after="input_combination")

# Combine res & f1
total <- rbind(res, f1)
total$input_combination = factor(total$input_combination,
                                 levels=c('S1', 'S2', 'S1L1', 'S2L1', 'S1L2', 'S2L2', 'S1L3', 'S2L3', 'L1', 'L2', 'L3',
                                          '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'),
                                 labels=c('Illumina Paired End', 'Illumina Mate Pair',
                                          'Illumina Paired End & MtSinai PacBio', 'Illumina Mate Pair & MtSinai PacBio',
                                          'Illumina Paired End & PacBio CCS', 'Illumina Mate Pair & PacBio CCS',
                                          'Illumina Paired End & 10X Genomics', 'Illumina Mate Pair & 10X Genomics',
                                          'MtSinai PacBio', 'PacBio CCS', '10X Genomics',
                                          '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'))

scale_custom <- list(
      # https://stackoverflow.com/questions/46803260/assigning-40-shapes-or-more-in-scale-shape-manual
      scale_shape_manual(values = c(15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 17, 46, 46, 46, 46, 46, 46, 46, 46, 46)),
      # https://www.farb-tabelle.de/de/farbtabelle.htm
      scale_color_manual(breaks = c('Illumina Paired End', 'Illumina Mate Pair',
                                    'Illumina Paired End & MtSinai PacBio', 'Illumina Mate Pair & MtSinai PacBio',
                                    'Illumina Paired End & PacBio CCS', 'Illumina Mate Pair & PacBio CCS',
                                    'Illumina Paired End & 10X Genomics', 'Illumina Mate Pair & 10X Genomics',
                                    'MtSinai PacBio', 'PacBio CCS', '10X Genomics'),
                         values = c("chartreuse3", "chartreuse1",
                                    "deepskyblue", "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4", "darkblue",
                                    "firebrick1", "firebrick3", "firebrick4",
                                    "tan", "tan", "tan", "tan", "tan", "tan", "tan", "tan", "tan"),
                         na.value = "gray80"))

ggplot(total, aes(recall, precision, color = input_combination)) +
      geom_point(size=0.5) +
      scale_custom +
      # Add labels to curves
      # cex: character expansion; rot: rotation
      geom_dl(aes(label = input_combination), method = list("last.points", rot=-20, cex = 0.5, dl.trans(x=x+0.1))) +
      geom_path() +
      labs(y = "Precision", x = "Recall", color = "Tool",
            title=plotTitle, subtitle="All combinations of short and long reads") +
      lims(x=c(0,1), y=c(0,1)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) +
      theme(plot.title = element_text(size=20, hjust=0.5, face="bold", color="black")) +
      theme(plot.subtitle = element_text(size=10, hjust=0.5, face="italic", color="black"))

ggsave(args[4], width=10, height=8)
