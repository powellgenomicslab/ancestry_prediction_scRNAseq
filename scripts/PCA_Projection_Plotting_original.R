#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

##### Read in Libraries #####
library(tidyverse)
library(ggpubr)
library(cluster)
library(RColorBrewer)
library(caret)
library(data.table)


args <- commandArgs(TRUE)
arguments <- read.table(args, header = F)
outdir <- arguments[1,]
data_score <- read_delim(as.character(arguments[2,]), delim = "\t", col_types = cols(.default = "d", "#FID" = "c", "IID" = "c", "Provided_Ancestry" = "c"))
onekg_score <- read_delim(as.character(arguments[3,]), delim = "\t", col_types = cols(.default = "d", "#IID" = "c", "SuperPop" = "c", "Population" = "c"))
onekg_anc <- read_delim(as.character(arguments[4,]), delim = "\t", col_types = cols(.default = "c", "SEX" = "d"))




##### Remove # from colnames #####
data_score <- data_score[!(is.na(data_score$PC1_AVG) | data_score$PC1_AVG == "NaN"),]
colnames(data_score) <- gsub("#", "",colnames(data_score))
colnames(onekg_score) <- gsub("#", "",colnames(onekg_score))
colnames(onekg_anc) <- gsub("#", "",colnames(onekg_anc))


print("data_score")
print(head(as.data.frame(data_score)))
print("onekg_score")
print(head(as.data.frame(onekg_score)))
print("onekg_anc")
print(head(as.data.frame(onekg_anc)))


##### Combine PCA Results Together #####
onekg_score_temp <- onekg_score[,c("IID", colnames(onekg_score)[grep("PC", colnames(onekg_score))])]
colnames(onekg_score_temp) <- c("IID",paste0("PC",1:(ncol(onekg_score_temp)-1)))
onekg_score_temp$FID <- NA
onekg_score_temp <- onekg_score_temp[,c("FID","IID",paste0("PC",1:(ncol(onekg_score_temp)-2)))]
onekg_score_temp <- left_join(onekg_score_temp, onekg_anc[,c("IID", "SuperPop")])
print(head(as.data.frame(onekg_score_temp)))


if (any(grepl("FID", colnames(data_score)))){
  data_score_temp <- data_score[,c(colnames(data_score)[grep("FID", colnames(data_score))], colnames(data_score)[grep("IID", colnames(data_score))], colnames(data_score)[grep("PC", colnames(data_score))])]
  colnames(data_score_temp) <- c("FID","IID",paste0("PC",1:(ncol(data_score_temp)-2)))
} else {
  data_score_temp <- data_score[,c(colnames(data_score)[grep("IID", colnames(data_score))], colnames(data_score)[grep("PC", colnames(data_score))])]
  colnames(data_score_temp) <- c("IID",paste0("PC",1:(ncol(data_score_temp)-1)))
  data_score_temp$FID <- NA
}
data_score_temp$SuperPop <- NA
print(head(as.data.frame(data_score_temp)))



model <- train(SuperPop ~ ., data = onekg_score_temp[,c("SuperPop", paste0("PC", 1:10))], method = "knn")
predictions <- predict(model, newdata = data_score_temp, type = "prob")
predictions$combined_assignment <- colnames(predictions)[max.col(predictions,ties.method="first")]

data_score_temp <- cbind(data_score_temp, predictions)


onekg_score_temp$AFR <- NA
onekg_score_temp$AMR <- NA
onekg_score_temp$EAS <- NA
onekg_score_temp$EUR <- NA
onekg_score_temp$SAS <- NA
onekg_score_temp$combined_assignment <- onekg_score_temp$SuperPop


scores <- rbind(onekg_score_temp, data_score_temp)
print(head(as.data.frame(scores)))

scores$IID <- as.character(scores$IID)
onekg_anc$IID <- as.character(onekg_anc$IID)
print(head(as.data.frame(scores)))



##### Plot Results #####
scores$Plot <- ifelse(is.na(scores$SuperPop), "Projected Data Assignments", "1000G Reference")
scores$Final_Assignment <- scores$combined_assignment


##### Set up population colors #####
pop_colors <- brewer.pal(length(unique(scores$Final_Assignment)), "Dark2")
names(pop_colors) <- unique(scores$Final_Assignment)


plot_PCs_medoids <- ggplot(scores, aes(PC1, PC2, color = Final_Assignment)) +
  geom_point() +
  theme_bw() +
  facet_wrap(vars(Plot)) +
  scale_color_manual(values = pop_colors)


ggsave(plot_PCs_medoids, filename = paste0(outdir,"Ancestry_PCAs.png"), height = 5, width = 12)


fwrite(scores[scores$Plot == "Projected Data Assignments", !(colnames(scores) %in% c("SuperPop", "combined_assignment", "Plot"))], paste0(outdir, "ancestry_assignments.tsv"), sep = "\t", na = "NA")


if (nrow(arguments) > 4){
  samples <- fread(as.character(arguments[5,]), sep = "\t")
  colnames(samples) <- c("Pool", "Individual")
  pool <- as.character(arguments[6,])
  indiv <- as.character(arguments[7,])

  if (samples$Pool[1] == pool & samples$Individual[1] == indiv & !is.na(samples)){
    fwrite(scores[scores$Plot == "1000G Reference", !(colnames(scores) %in% c("SuperPop", "combined_assignment", "Plot"))], paste0(outdir, "ancestry_assignments_w_ref.tsv"), sep = "\t", na = "NA")
  }
}