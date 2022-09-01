#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

##### Read in Libraries #####
library(tidyverse)
library(data.table)
library(scales)
library(RColorBrewer)
library(colorspace)
library(grid)


args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
samples_file <- args[2]




print("Reading in samples dataframe.")
samples <- fread(samples_file, sep = "\t")
colnames(samples) <- c("Pool", "Individual")


##### Read in the ancestry prediction dataframes #####
dt_list <- list()

for (row in 1:nrow(samples)){
    dt_list[[paste0(samples[row,"Pool"], "_", samples[row,"Individual"])]] <- fread(paste0(outdir, "/", samples[row,"Pool"], "/individual_", samples[row,"Individual"], "/pca_sex_checks_original/ancestry_assignments.tsv"), sep = "\t")
    dt_list[[paste0(samples[row,"Pool"], "_", samples[row,"Individual"])]]$Pool <- samples[row,"Pool"]
    dt_list[[paste0(samples[row,"Pool"], "_", samples[row,"Individual"])]]$IID <- samples[row,"Individual"]
}



##### Account for possibility that there may only be one individual in one pool at one location
message("Combining predictions together.")
if (nrow(samples) > 1 ){

    dt <- do.call(rbind, dt_list)

} else {

    dt <- dt_list[[1]]

}


fwrite(dt, paste0(outdir, "/ancestry_assignments.tsv"), sep = "\t")




##### Make a PCA figure with all data combined #####
message("Reading in reference 1000 genomes data.")
ref_1000g <- fread(paste0(outdir, "/", samples[1,"Pool"], "/individual_", samples[1,"Individual"], "/pca_sex_checks_original/ancestry_assignments_w_ref.tsv"), sep = "\t")
ref_1000g$Pool <- NA 
ref_1000g$Plotting_Ancestry <- ref_1000g$Final_Assignment

print(ref_1000g)

###### Make different dts for each predicted ancestry #####
# dt$Plotting_Ancestry <- paste0("Predicted ", dt$Final_Assignment)
dt$Plotting_Ancestry <- "Predicted"


message("Creating comparison facets.")

dt_ref_list <- list()

dt_ref_list[["all"]] <- rbind(ref_1000g, dt)
dt_ref_list[["all"]]$Group <- "1000G Reference +\nAll Predicted Data"

for (anc in unique(dt$Final_Assignment)){
    dt_ref_list[[anc]] <- rbind(ref_1000g, dt[Final_Assignment == anc])
    dt_ref_list[[anc]]$Group <- paste0("1000G Reference +\n", anc, " Predicted Data")
}

print(dt_ref_list)

##### Combine datatables #####
message("Binding data.tables.")
dt_ref <- do.call(rbind, dt_ref_list)

print(dt_ref)


##### Plot Results #####
##### Set up population colors #####
base_colors <- brewer.pal((length(unique(dt_ref$Plotting_Ancestry)) - 1), "Dark2")
names(base_colors) <- unique(dt_ref$Plotting_Ancestry)[!(unique(dt_ref$Plotting_Ancestry) %in% "Predicted")]

pop_colors <- c(alpha(base_colors, 0.2), "black")
names(pop_colors) <- c(unique(ref_1000g$Plotting_Ancestry), "Predicted")


plot_PCs_medoids <- ggplot(dt_ref, aes(PC1, PC2, color = Plotting_Ancestry)) +
  geom_point(size = 0.5) +
  theme_bw() +
  facet_wrap(vars(factor(Group, levels =  c("1000G Reference +\nAll Predicted Data", paste0("1000G Reference +\n", unique(dt$Final_Assignment), " Predicted Data")))), nrow = 1) +
  scale_color_manual(values = pop_colors)



pred_colors <- c("white",base_colors[unique(dt$Final_Assignment)])
names(pred_colors) <- c("1000G Reference +\nAll Predicted Data", paste0("1000G Reference +\n", unique(dt$Final_Assignment), " Predicted Data"))



# Generate the ggplot2 plot grob
g <- grid.force(ggplotGrob(plot_PCs_medoids))
# Get the names of grobs and their gPaths into a data.frame structure
grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))
# Build optimal gPaths that will be later used to identify grobs and edit them
grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")
grobs_df$gPath_full <- gsub(pattern = "layout::", 
                            replacement = "", 
                            x = grobs_df$gPath_full, 
                            fixed = TRUE)


# Get the gPaths of the strip background grobs
strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = ".*strip\\.background.x.*", 
                                            x = grobs_df$gPath_full)]


# Edit the grobs
for (i in 1:length(strip_bg_gpath)){

  g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = pred_colors[i]))
}


ggsave(g, filename = paste0(outdir,"/Ancestry_PCAs.png"), height = 3, width = 1.5 + 2.5 * (length(unique(dt$Final_Assignment)) + 1))
