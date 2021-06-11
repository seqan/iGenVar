library(tidyverse)

# Rscript --vanilla scripts/plot_all.R {input} {parameter_name} {output}

args = commandArgs(trailingOnly=TRUE)

parameter_name = args[2]

if (parameter_name == "min_var_length") {
  parameter_values = c(10,30,50,100,500,1000,2000)
} else if (parameter_name == "max_var_length") {
  parameter_values = c(500000,1000000,2000000)
} else if (parameter_name == "max_tol_inserted_length") {
  parameter_values = c(1,5,10)
} else if (parameter_name == "max_overlap") {
  parameter_values = c(5,10,20)
} else {                    # hierarchical_clustering_cutoff
  parameter_values = c(20,50,100,200,500,5000)
}

# res <- read_tsv("results/truvari/min_var_length/all_results.txt", col_names = c(parameter_name, "min_qual", "metric", "value"))
res <- read_tsv(args[1], col_names = c(parameter_name, "min_qual", "metric", "value"))
res[[parameter_name]] = factor(res[[parameter_name]], levels=parameter_values, labels=as.character(parameter_values))

res %>%
    filter(metric %in% c("recall", "precision")) %>%
    pivot_wider(names_from=metric, values_from=value) %>%
    filter(recall!=0 | precision!=0) %>%
    mutate(precision = 100*precision, recall = 100*recall) %>%
    ggplot(aes(recall, precision, color = get(parameter_name), pch = get(parameter_name))) +
      geom_point(size=0.5) +
      geom_path() +
      labs(y = "Precision", x = "Recall", color = parameter_name, pch = parameter_name) +
      lims(x=c(0,100), y=c(0,100)) +
      theme_bw() +
      theme(text = element_text(size=14), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))

ggsave(args[3], width=7, height=6)
