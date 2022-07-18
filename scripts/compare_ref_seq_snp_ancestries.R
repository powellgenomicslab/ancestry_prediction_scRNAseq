library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(caret)
library(colorspace)
library(ggpubr)
library(grid)


### Set up directories ###
args <- commandArgs(TRUE)
arguments <- read.table(args, header = F, sep p= "\t")
datadir <- arguments[1,]
ref_predictions <- arguments[2,]
outdir <- arguments[3,]
pools <- c(arguments[4,])



### Read in reference SNP-based ancestry predictions ###
ref_df <- fread(ref_predictions)



## Read in predictions ##
## Souporcell ##
souporcell_anc_list <- lapply(pools, function(pool){
    fread(paste0(datadir, pool, "/souporcell/pca_sex_checks_original/ancestry_assignments.tsv"))
})
names(souporcell_anc_list) <- pools


## Add pool name to dataframes ##
## Souporcell ##
souporcell_anc_list <- lapply(names(souporcell_anc_list), function(pool){
    souporcell_anc_list[[pool]]$Pool <- pool
    souporcell_anc_list[[pool]]$Software <- "souporcell"
    return(souporcell_anc_list[[pool]])
})
names(souporcell_anc_list) <- pools


## read in cluster to individual keys ##
souporcell_id_conversions <- lapply(pools, function(pool){
    fread(paste0(datadir, pool, "/souporcell/Genotype_ID_key.txt"), sep = "\t")
})
names(souporcell_id_conversions) <- pools


## Add IDs to ancestry assignments ##
## Souporcell ##
souporcell_anc_list <- lapply(names(souporcell_anc_list), function(pool){
    souporcell_anc_list[[pool]]$IID <- gsub("donor", "", souporcell_anc_list[[pool]]$IID)
    souporcell_anc_list[[pool]] <- data.table(souporcell_id_conversions[[pool]])[souporcell_anc_list[[pool]], on = c("Cluster_ID" = "IID")]
    return(souporcell_anc_list[[pool]])
})


## combine all dataframes together ##
## Souporcell ##
souporcell_anc <- do.call(rbind, souporcell_anc_list)

souporcell_anc_sub <- souporcell_anc[,c("Genotype_ID", "Cluster_ID", "Pool", "AFR", "AMR", "EAS", "EUR", "SAS", "Final_Assignment")]
colnames(souporcell_anc_sub) <- c("IID", "Cluster_ID", "Pool", "AFR", "AMR", "EAS", "EUR", "SAS", "seq_snp_ancestry")
souporcell_anc_sub_joined <- merge(ref_df, souporcell_anc_sub, by = "IID", all = TRUE)
souporcell_anc_sub_joined$seq_snp_ancestry <- ifelse(is.na(souporcell_anc_sub_joined$seq_snp_ancestry), "unidentified_in_pool", souporcell_anc_sub_joined$seq_snp_ancestry)

souporcell_anc_sub_joined$Software <- "Souporcell"


anc_joined <- souporcell_anc_sub_joined


anc_joined$seq_snp_ancestry <- factor(anc_joined$seq_snp_ancestry, levels = c("EUR", "EAS", "AMR", "SAS", "AFR"))
anc_joined$ref_snp_ancestry <- factor(anc_joined$ref_snp_ancestry, levels = c("EUR", "EAS", "AMR", "SAS", "AFR"))


fwrite(anc_joined, paste0(outdir, "snp_ancestry_predictions.tsv"), sep = "\t")


### Make Figures for Presentation ###
colors <- c(brewer.pal(5, "Dark2"), "grey70")
names(colors) = c("EUR", "EAS", "AMR", "SAS", "AFR","unidentified_in_pool")

ref_plot <- ggplot(unique(anc_joined[,c("IID", "Pool", "ref_snp_ancestry")]), aes(ref_snp_ancestry, fill = ref_snp_ancestry)) +
    geom_bar() +
    theme_classic() +
    scale_fill_manual(values = colors)+
    stat_count(geom = "text", colour = "black", size = 3.5,
        aes(label = ..count..),vjust=0) +
        xlab("Ancestry")

ggsave(ref_plot, filename = paste0(outdir, "reference_ancestry_numbers.png"), width =3, height = 3)



correct_colors <- c(Correct = "grey30", Incorrect = "firebrick3")


seq_plot_both_correct <- ggplot(anc_joined, aes(seq_snp_ancestry, fill = Annotation)) +
    geom_bar() +
    theme_classic() +
    scale_fill_manual(values = correct_colors) +
    facet_grid(~Software) +
    stat_count(geom = "text", colour = "black", size = 3.5,
        aes(label = ..count..),vjust=0)

ggsave(seq_plot_both_correct, filename = paste0(outdir, "predicted_ancestry_numbers_correct.png"), width = 7.5, height = 3)



stat_metrics <- list()

for (soft in unique(anc_joined$Software)){
    temp <- anc_joined[Software == soft]

    ##### Calculate TPR, TNR, FPR, FNR, accuracy, balanced accuracy, MCC #####
    stat_metrics[[soft]] <- data.table(Group = unique(c(as.character(temp$seq_snp_ancestry), as.character(temp$ref_snp_ancestry))), 
                TP = as.numeric(NA),
                TN = as.numeric(NA),
                FP = as.numeric(NA),
                FN = as.numeric(NA),
                TPR = as.numeric(NA),
                TNR = as.numeric(NA),
                FPR = as.numeric(NA),
                FNR = as.numeric(NA),
                PPV = as.numeric(NA),
                NPV = as.numeric(NA),
                Accuracy = as.numeric(NA),
                Balanced_Accuracy = as.numeric(NA),
                MCC = as.numeric(NA)
                )


    for (group in stat_metrics[[soft]]$Group){
        

        stat_metrics[[soft]][Group == group, "TP"] <- nrow(temp[seq_snp_ancestry == group & Annotation == "Correct"])
        stat_metrics[[soft]][Group == group, "TN"] <- nrow(temp[seq_snp_ancestry != group & Annotation == "Correct"])
        stat_metrics[[soft]][Group == group, "FP"] <- nrow(temp[seq_snp_ancestry == group & Annotation == "Incorrect"])
        stat_metrics[[soft]][Group == group, "FN"] <- nrow(temp[seq_snp_ancestry != group & Annotation == "Incorrect"])


        stat_metrics[[soft]][Group == group, "TPR"] <- stat_metrics[[soft]][Group == group]$TP/(stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$FN) ## Sensitivity
        stat_metrics[[soft]][Group == group, "TNR"] <- stat_metrics[[soft]][Group == group]$TN/(stat_metrics[[soft]][Group == group]$TN + stat_metrics[[soft]][Group == group]$FP) ## Specificity
        stat_metrics[[soft]][Group == group, "FPR"] <- stat_metrics[[soft]][Group == group]$FP/(stat_metrics[[soft]][Group == group]$TN + stat_metrics[[soft]][Group == group]$FP)
        stat_metrics[[soft]][Group == group, "FNR"] <- stat_metrics[[soft]][Group == group]$FN/(stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$FN)
        stat_metrics[[soft]][Group == group, "PPV"] <- stat_metrics[[soft]][Group == group]$TP/(stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$FP)
        stat_metrics[[soft]][Group == group, "NPV"] <- stat_metrics[[soft]][Group == group]$TN/(stat_metrics[[soft]][Group == group]$TN + stat_metrics[[soft]][Group == group]$FN)
        stat_metrics[[soft]][Group == group, "Accuracy"] <- (stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$TN)/nrow(temp)
        stat_metrics[[soft]][Group == group, "Balanced_Accuracy"] <- (stat_metrics[[soft]][Group == group]$TPR + stat_metrics[[soft]][Group == group]$TNR)/2
        stat_metrics[[soft]][Group == group, "MCC"] <- (stat_metrics[[soft]][Group == group]$TP * stat_metrics[[soft]][Group == group]$TN - stat_metrics[[soft]][Group == group]$FP * stat_metrics[[soft]][Group == group]$FN)/sqrt((stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$FP)*(stat_metrics[[soft]][Group == group]$TP + stat_metrics[[soft]][Group == group]$FN)*(stat_metrics[[soft]][Group == group]$TN + stat_metrics[[soft]][Group == group]$FP)*(stat_metrics[[soft]][Group == group]$TN + stat_metrics[[soft]][Group == group]$FN))

    }
    stat_metrics[[soft]]$Software <- soft

}



stat_metrics_dt <- do.call(rbind, stat_metrics)

stat_metrics_long_dt <- melt(stat_metrics_dt, id.vars = c("Group", "Software"))



stat_metrics_long_dt$Expected_Direction <- ifelse(stat_metrics_long_dt$variable %in% c("FP", "FN", "FNR", "FPR"), "Low", "High")
stat_metrics_long_dt$Group <- factor(stat_metrics_long_dt$Group, levels = c("EUR", "EAS", "AMR", "SAS", "AFR"))


stat_heat <- ggplot(stat_metrics_long_dt[!(variable %in% c("TP", "TN", "FP", "FN"))], aes(Group, variable, fill = value)) +
    geom_tile() +
    theme_classic() +
    geom_text(aes(label=round(value, 2)), color = "white", size = 3) +
    scale_fill_continuous_sequential(palette = "Terrain 2", na.value="white") +
    facet_grid(Expected_Direction~Software, scales = "free_y", space='free') +
    xlab("Ancestry")


ggsave(stat_heat, filename = paste0(outdir, "statistics_heatmap.png"), width = 7.5, height = 3.5)



combined_plot <- cowplot::plot_grid(seq_plot_both_correct + theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle ("MoST Sample Ancestry Annotations"), stat_heat + theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
          ncol = 1, nrow = 2, align = "v", axis = "l", rel_heights = c(1,0.75))

ggsave(combined_plot, filename = paste0(outdir, "predicted_numbers_statistics_heatmap_combined.png"), width = 7, height = 5)


##### Plot the incorrect vs correct per reference group predictions #####
anc_joined$AFR <- as.numeric(anc_joined$AFR)
anc_joined$AMR <- as.numeric(anc_joined$AMR)
anc_joined$EAS <- as.numeric(anc_joined$EAS)
anc_joined$EUR <- as.numeric(anc_joined$EUR)
anc_joined$SAS <- as.numeric(anc_joined$SAS)

anc_joined_long <- melt(anc_joined, measure.vars = c("AFR", "AMR", "EAS", "EUR", "SAS"))




### Add reference data ###
ref_anc_sub <- ref_df[,c("IID","AFR", "AMR", "EAS", "EUR", "SAS", "combined_assignment")]
colnames(ref_anc_sub) <- c("IID","AFR", "AMR", "EAS", "EUR", "SAS", "seq_snp_ancestry")


ref_anc_sub$AFR <- as.numeric(ref_anc_sub$AFR)
ref_anc_sub$AMR <- as.numeric(ref_anc_sub$AMR)
ref_anc_sub$EAS <- as.numeric(ref_anc_sub$EAS)
ref_anc_sub$EUR <- as.numeric(ref_anc_sub$EUR)
ref_anc_sub$SAS <- as.numeric(ref_anc_sub$SAS)
ref_anc_sub$ref_snp_ancestry <- ref_anc_sub$seq_snp_ancestry

ref_anc_sub <- ref_anc_sub[unique(anc_joined[,c("IID", "Pool")]), on = "IID"]


ref_anc_sub_long <- melt(ref_anc_sub, measure.vars = c("AFR", "AMR", "EAS", "EUR", "SAS"))
ref_anc_sub_long$Software <- "Reference"
ref_anc_sub_long$Annotation <- "Correct"


anc_joined_long_ref <- rbind(anc_joined_long[,c("IID", "seq_snp_ancestry", "ref_snp_ancestry", "Pool", "variable", "value", "Software", "Annotation")], ref_anc_sub_long)


### reference annotation was upadted by new method of annotation so need to update group ###
anc_joined_long_ref$ref_snp_ancestry <- factor(anc_joined_long_ref$ref_snp_ancestry, levels = c("EUR", "EAS", "AMR", "SAS", "AFR"))

anc_joined_long_ref$Annotation <- ifelse(anc_joined_long_ref$ref_snp_ancestry == anc_joined_long_ref$seq_snp_ancestry, "Correct", "Incorrect")


### Make figure of stacked combined with ref ### 
p_probabilities_ref <- ggplot(anc_joined_long_ref, aes(paste0(gsub("-\\d", "", IID), "-", gsub("RZ731_", "", Pool)), value, fill = variable)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    facet_grid(factor(gsub(" ", "\n", Software), levels =  c("Reference", "Freemuxlet", "Freemuxlet\nImpute", "Souporcell", "Souporcell\nImpute")) ~ ref_snp_ancestry, scales="free", space="free_x") +
    scale_fill_manual(values = colors[1:5]) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.x=element_blank()) +
    geom_text(aes(label = ifelse(Annotation == "Incorrect", "X", ""), y = 0.9),  vjust = 1) +
    ylab("Probability")

# Generate the ggplot2 plot grob
g_ref <- grid.force(ggplotGrob(p_probabilities_ref))
# Get the names of grobs and their gPaths into a data.frame structure
grobs_df_ref <- do.call(cbind.data.frame, grid.ls(g_ref, print = FALSE))
# Build optimal gPaths that will be later used to identify grobs and edit them
grobs_df_ref$gPath_full <- paste(grobs_df_ref$gPath, grobs_df_ref$name, sep = "::")
grobs_df_ref$gPath_full <- gsub(pattern = "layout::", 
                            replacement = "", 
                            x = grobs_df_ref$gPath_full, 
                            fixed = TRUE)


# Get the gPaths of the strip background grobs
strip_bg_gpath_ref <- grobs_df_ref$gPath_full[grepl(pattern = ".*strip\\.background.x.*", 
                                            x = grobs_df_ref$gPath_full)]


# Edit the grobs
for (i in 1:length(strip_bg_gpath_ref)){
  g_ref <- editGrob(grob = g_ref, gPath = strip_bg_gpath_ref[i], gp = gpar(fill = colors[levels(anc_joined_long_ref$ref_snp_ancestry)[i]]))
}

ggsave(g_ref, filename = paste0(outdir,"assignments_probabilities_w_ref.png"), width = 10, height = 7)

