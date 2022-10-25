# For the first time, you may need to run:
# install.packages("directlabels")
# install.packages("tidyverse")

library(directlabels)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

res <- read_tsv(args[1], col_names = c("input_combination", "min_qual", "metric", "value"))
f1 <- read_csv(args[2], col_names = c("input_combination", "precision", "recall"))

res <- res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0)

# Creating an empty column:
f1 <- add_column(f1, min_qual = NA, .after="input_combination")

# Combine res & f1
total <- rbind(res, f1)
total$input_combination = factor(total$input_combination,
                                 levels=c('hg38_Sim_default', 'hg38_Sim_InDel', 'hg38_Sim_noSNP', 'hg38_Sim_SNPandSV',
                                          '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'),
                                 labels=c('hg38_Sim_default', 'hg38_Sim_InDel', 'hg38_Sim_noSNP', 'hg38_Sim_SNPandSV',
                                          '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'))

scale_custom <- list(
      # https://stackoverflow.com/questions/46803260/assigning-40-shapes-or-more-in-scale-shape-manual
      scale_shape_manual(values = c(15, 15, 16, 16,
                                    46, 46, 46, 46, 46, 46, 46, 46, 46)),
      # https://www.farb-tabelle.de/de/farbtabelle.htm
      scale_color_manual(breaks = c('hg38_Sim_default', 'hg38_Sim_InDel', 'hg38_Sim_noSNP', 'hg38_Sim_SNPandSV'),
                         values = c("chartreuse1", "deepskyblue", "darkorchid1", "maroon1",
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
            title="Caller Comparison - Vaquita LR", subtitle="With Mason simulated Illumina short reads") +
      lims(x=c(0,1), y=c(0,1)) +
      theme_bw() +
      theme(panel.spacing = unit(0.75, "lines")) +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) +
      theme(plot.title = element_text(size=20, hjust=0.5, face="bold", color="black")) +
      theme(plot.subtitle = element_text(size=10, hjust=0.5, face="italic", color="black"))

ggsave(args[3], width=10, height=8)
