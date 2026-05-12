#dir
setwd("work/dir/path")

# Packages
library(readr)
library(dplyr)
library(ggplot2)
#Import

summary <- read_csv("summary.csv")
View(summary)

#Arguments 

ligands <- c("IHP", "2_3_BPG")
modern <- c("1A3N", "2DN3", "6BB5")
ancestral <- c("AncA", "AncB", "AncMH_nat")
mutant <- c("mut_D99N", "mut_K82N", "mut_N102T", "mut_Y35F")


#Filter
filtered_summary <- subset(summary, ligand %in% ligands)
View(filtered_summary)

ihp_summ <- subset(summary, ligand %in% ligands[1])
View(ihp_summ)

bpg_summ <- subset(summary, ligand %in% ligands[2])
View(bpg_summ)


#Output
write.csv(ihp_summ, "ihp_summ.csv", row.names = FALSE)
write.csv(bpg_summ, "bpg_summ.csv", row.names = FALSE)

#Box Plots


make_boxplot <- function(data, ligand_name) {

  ggplot(data, aes(x = protein, y = mean_r)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
    labs(
      title = paste0("Centroid dispersion (mean_r) — ", ligand_name),
      x = "Protein",
      y = "Mean radial centroid dispersion (Å)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_ihp <- make_boxplot(ihp_summ, "IHP")
plot_bpg <- make_boxplot(bpg_summ, "2,3-BPG")

ggsave(file.path(".", "boxplot_IHP_mean_r.png"), plot_ihp, width = 10, height = 5, dpi = 300)
ggsave(file.path(".", "boxplot_2_3_BPG_mean_r.png"), plot_bpg, width = 10, height = 5, dpi = 300)











