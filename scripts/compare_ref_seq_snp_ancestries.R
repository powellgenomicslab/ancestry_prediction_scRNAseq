#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

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
arguments <- read.table(args, header = F, sep = "\t")
datadir <- arguments[1,]
ref_predictions <- arguments[2,]
outdir <- arguments[3,]
samples_file <- c(arguments[4,])


print("Reading in samples dataframe.")
samples <- fread(samples_file, sep = "\t")
colnames(samples) <- c("Pool", "Individual")
pools <- unique(samples$Pool)
indiv <- unique(samples$Individual)



individuals <- lapply(pools, function(pool){
    dir(paste0(datadir,"/",pool), pattern = "individual")
})
names(individuals) <- pools



### Read in reference SNP-based ancestry predictions ###
ref_df <- fread(ref_predictions)
ref_df <- ref_df[IID %in% indiv]


## Read in predictions ##
print("Reading in results.")

freebayes_anc_list <- lapply(pools, function(pool){
    tmp <- lapply(individuals[[pool]], function(indiv){
        if (file.exists(paste0(datadir, "/", pool,"/",indiv, "/pca_sex_checks_original/ancestry_assignments.tsv"))){
            tmp2 <- fread(paste0(datadir, "/", pool,"/",indiv, "/pca_sex_checks_original/ancestry_assignments.tsv"))
            tmp2$Pool <- pool
            tmp2$IID <- gsub("individual_", "", indiv)
            return(tmp2)
        } else {
            NULL
        }
    })
    names(tmp) <- individuals[[pool]]
    return(tmp)
})
names(freebayes_anc_list) <- pools


## combine all dataframes together ##
print("Updating data.")
freebayes_anc <- do.call(rbind, lapply(freebayes_anc_list, function(x) do.call(rbind, x)))


freebayes_anc_sub <- freebayes_anc[,c("IID", "Pool", "AFR", "AMR", "EAS", "EUR", "SAS", "Final_Assignment")]
colnames(freebayes_anc_sub) <- c("IID", "Pool", "AFR", "AMR", "EAS", "EUR", "SAS", "seq_snp_ancestry")


print("Combining data with reference.")
anc_joined <- merge(ref_df[,c("IID", "Final_Assignment")], freebayes_anc_sub, by = "IID", all = TRUE)
colnames(anc_joined) <- gsub("Final_Assignment", "ref_snp_ancestry", colnames(anc_joined))

anc_joined$Software <- "Freebayes"
anc_joined$Annotation <- ifelse(as.character(anc_joined$seq_snp_ancestry) == as.character(anc_joined$ref_snp_ancestry), "Correct", "Incorrect")


anc_joined$seq_snp_ancestry <- ifelse(is.na(anc_joined$seq_snp_ancestry), "unidentified_in_pool", as.character(anc_joined$seq_snp_ancestry))
anc_joined$seq_snp_ancestry <- factor(anc_joined$seq_snp_ancestry, levels = c("EUR", "EAS", "AMR", "SAS", "AFR", "unidentified_in_pool"))
anc_joined$ref_snp_ancestry <- factor(anc_joined$ref_snp_ancestry, levels = c("EUR", "EAS", "AMR", "SAS", "AFR"))


fwrite(anc_joined, paste0(outdir, "snp_ancestry_predictions.tsv"), sep = "\t")


### Make Figures for Presentation ###
colors <- c(brewer.pal(5, "Dark2"), "grey70")
names(colors) = c("EUR", "EAS", "AMR", "SAS", "AFR","unidentified_in_pool")

ref_plot <- ggplot(unique(anc_joined[!is.na(ref_snp_ancestry),c("IID", "Pool", "ref_snp_ancestry")]), aes(ref_snp_ancestry, fill = ref_snp_ancestry)) +
    geom_bar() +
    theme_classic() +
    scale_fill_manual(values = colors[1:5])+
    stat_count(geom = "text", colour = "black", size = 3.5,
        aes(label = ..count..),vjust=0) +
        xlab("Ancestry")

ggsave(ref_plot, filename = paste0(outdir, "reference_ancestry_numbers.png"), width =2.5 + length(unique(unique(anc_joined[!is.na(ref_snp_ancestry),c("IID", "Pool", "ref_snp_ancestry")])$ref_snp_ancestry))*0.25, height = 3)



correct_colors <- c(Correct = "grey30", Incorrect = "firebrick3")


seq_plot_both_correct <- ggplot(anc_joined, aes(seq_snp_ancestry, fill = Annotation)) +
    geom_bar() +
    theme_classic() +
    scale_fill_manual(values = correct_colors) +
    stat_count(geom = "text", colour = "black", size = 3.5,
        aes(label = ..count..),vjust=0) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.x=element_blank())

ggsave(seq_plot_both_correct, filename = paste0(outdir, "predicted_ancestry_numbers_correct.png"), width =2 + length(unique(anc_joined$seq_snp_ancestry))*0.25, height = 4)


seq_plot_both_correct_identified <- ggplot(anc_joined[seq_snp_ancestry != "unidentified_in_pool"], aes(seq_snp_ancestry, fill = Annotation)) +
    geom_bar() +
    theme_classic() +
    scale_fill_manual(values = correct_colors) +
    stat_count(geom = "text", colour = "black", size = 3.5,
        aes(label = ..count..),vjust=0) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.x=element_blank())

ggsave(seq_plot_both_correct_identified, filename = paste0(outdir, "predicted_ancestry_numbers_correct_identified.png"), width =2 + length(unique(anc_joined[seq_snp_ancestry != "unidentified_in_pool"]$seq_snp_ancestry))*0.25, height = 4)




anc_joined_noNA <- anc_joined[!is.na(Pool) & !is.na(IID)]


stat_metrics <- list()

for (soft in unique(anc_joined_noNA$Software)){
    temp <- anc_joined_noNA[Software == soft]

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


stat_heat <- ggplot(stat_metrics_long_dt[!(variable %in% c("TP", "TN", "FP", "FN")) & !is.na(Group)], aes(Group, variable, fill = value)) +
    geom_tile() +
    theme_classic() +
    geom_text(aes(label=round(value, 2)), color = "white", size = 3) +
    scale_fill_continuous_sequential(palette = "Terrain 2", na.value="white") +
    facet_grid(vars(Expected_Direction), scales = "free_y", space='free') +
    xlab("Ancestry")


ggsave(stat_heat, filename = paste0(outdir, "statistics_heatmap.png"), width = 4 + length(unique(stat_metrics_long_dt[!(variable %in% c("TP", "TN", "FP", "FN")) & !is.na(Group)]$Group))*0.25, height = 3.5)



combined_plot <- cowplot::plot_grid(seq_plot_both_correct_identified + theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), stat_heat + theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
          ncol = 1, nrow = 2, align = "v", axis = "l", rel_heights = c(1,0.75))

ggsave(combined_plot, filename = paste0(outdir, "predicted_numbers_statistics_heatmap_combined.png"), width = 3 + length(unique(stat_metrics_long_dt[!(variable %in% c("TP", "TN", "FP", "FN")) & !is.na(Group)]$Group))*0.25, height = 5)


##### Plot the incorrect vs correct per reference group predictions #####
anc_joined$AFR <- as.numeric(anc_joined$AFR)
anc_joined$AMR <- as.numeric(anc_joined$AMR)
anc_joined$EAS <- as.numeric(anc_joined$EAS)
anc_joined$EUR <- as.numeric(anc_joined$EUR)
anc_joined$SAS <- as.numeric(anc_joined$SAS)

anc_joined_long <- melt(anc_joined, measure.vars = c("AFR", "AMR", "EAS", "EUR", "SAS"))




### Add reference data ###
ref_anc_sub <- ref_df[,c("IID","AFR", "AMR", "EAS", "EUR", "SAS", "Final_Assignment")]
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

anc_joined_long_ref$Annotation <- ifelse(as.character(anc_joined_long_ref$ref_snp_ancestry) == as.character(anc_joined_long_ref$seq_snp_ancestry), "Correct", "Incorrect")


### Make figure of stacked combined with ref ### 
base_colors <- colors[1:5]

p_probabilities_ref <- ggplot(anc_joined_long_ref, aes(paste0(IID, "-",Pool), value, fill = variable)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_classic() +
    facet_grid(factor(gsub(" ", "\n", Software), levels =  c("Reference", "Freebayes")) ~ ref_snp_ancestry, scales="free", space="free_x") +
    scale_fill_manual(values = base_colors) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.x=element_blank()) +
    geom_text(aes(label = ifelse(Annotation == "Incorrect", "X", ""), y = 0.9),  vjust = 1) +
    ylab("Probability")


pred_colors <- base_colors[levels(anc_joined_long_ref$ref_snp_ancestry)[levels(anc_joined_long_ref$ref_snp_ancestry) %in% unique(anc_joined_long_ref$ref_snp_ancestry)]]


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
if (any(is.na(anc_joined_long_ref$ref_snp_ancestry))) {
    for (i in 1:(length(strip_bg_gpath_ref) - 1)){
        g_ref <- editGrob(grob = g_ref, gPath = strip_bg_gpath_ref[i], gp = gpar(fill = pred_colors[i]))
    }
        g_ref <- editGrob(grob = g_ref, gPath = strip_bg_gpath_ref[i+1], gp = gpar(fill = "white"))

} else {
    for (i in 1:length(strip_bg_gpath_ref)){
        g_ref <- editGrob(grob = g_ref, gPath = strip_bg_gpath_ref[i], gp = gpar(fill = pred_colors[i]))
    }
}

ggsave(g_ref, filename = paste0(outdir,"assignments_probabilities_w_ref.png"), width = 1 + 0.15 * nrow(unique(anc_joined_long_ref[,c("IID", "Pool")])), height = 4.5)





