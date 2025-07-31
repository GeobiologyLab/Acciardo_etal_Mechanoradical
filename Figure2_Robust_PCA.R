library(reshape2)
library(ggplot2)
library(dplyr)
library(vegan)
library(dada2)
library(phyloseq)
library(compositions)
library(rrcov)

load("Data/Mechanoradical_Incubation.RData")

### Remove all ASVs found in controls
controlSamples <- rownames(df)[grepl("Blank", rownames(df))]
largeVolume <- rownames(df)[grepl("1L", rownames(df))]

df <- df[-which(rownames(df)%in%largeVolume),]
contaminants <- which(colSums(df[controlSamples,])>0)
df.clean <- df[,-contaminants]
df.clean <- df.clean[,-which(colSums(df.clean)==0)]
snvTable <- df.clean[-which(rownames(df.clean)%in%controlSamples),]

taxa$SNV <- rownames(taxa)

snvID <- gsub(".*-(\\w+)_.*", "\\1", rownames(snvTable))
snvID <- sub("(_S).*", "", snvID)

tmp <- snvTable
rownames(tmp) <- snvID
tmp <- tmp/rowSums(tmp) #result is same if normalized to relative abundance or counts just put in

### Robust PCA analysis

transformedSNVs <- decostand(tmp,method="rclr")
rpca <- PcaHubert(transformedSNVs)

eigenvalues <- rpca@eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues) * 100  


### Plotting by timpoint

scores <- rpca@scores
scores_df <- data.frame(PC1 = scores[, 1], PC2 = scores[, 2], samples = rownames(tmp))

scores_df$Treatment = "t0"
scores_df$Treatment[which(grepl("W",scores_df$samples))] <- "Water"
scores_df$Treatment[which(grepl("G",scores_df$samples))] <- "Glass"
scores_df$Treatment[which(grepl("TAn",scores_df$samples))] <- "Rock"
scores_df$Timepoint <- as.numeric(sub(".*?(\\d+).*", "\\1", scores_df$samples))

scores_df$Timepoint[which(scores_df$Timepoint==0)]<-5


scores_df$Treatment <- factor(scores_df$Treatment, levels = c("t0", "Water", "Rock", "Glass"))

xlim <- c(-20,20)
ylim <- c(-10,15)

ggplot(scores_df[which(scores_df$Timepoint%in%c(1) | scores_df$Treatment=="t0" ),], aes(x = PC1, y = PC2, group = Treatment,alpha=Timepoint)) +
  geom_point(aes(shape=Treatment,color=Treatment),size = 3) +
  scale_shape_manual(values = c("t0" = 8, "Water" = 15, "Rock" = 16, "Glass" = 17)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values = c("t0"="black","Water" = "blue","Rock"="red","Glass"="orange")) +
  labs(
    x = paste("PC1 (", round(variance_explained[1], 2), "% variance)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2], 2), "% variance)", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16) 
  ) +
  scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)

ggplot(scores_df[which(scores_df$Timepoint%in%c(1,2) | scores_df$Treatment=="t0" ),], aes(x = PC1, y = PC2, group = Treatment,alpha=Timepoint)) +
  geom_point(aes(shape=Treatment,color=Treatment),size = 3) +
  scale_shape_manual(values = c("t0" = 8, "Water" = 15, "Rock" = 16, "Glass" = 17)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values = c("t0"="black","Water" = "blue","Rock"="red","Glass"="orange")) +
  labs(
    x = paste("PC1 (", round(variance_explained[1], 2), "% variance)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2], 2), "% variance)", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16) 
  ) +
  scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)

ggplot(scores_df[which(scores_df$Timepoint%in%c(1,2,3) | scores_df$Treatment=="t0" ),], aes(x = PC1, y = PC2, group = Treatment,alpha=Timepoint)) +
  geom_point(aes(shape=Treatment,color=Treatment),size = 3) +
  scale_shape_manual(values = c("t0" = 8, "Water" = 15, "Rock" = 16, "Glass" = 17)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values = c("t0"="black","Water" = "blue","Rock"="red","Glass"="orange")) +
  labs(
    x = paste("PC1 (", round(variance_explained[1], 2), "% variance)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2], 2), "% variance)", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16)  
  ) +
  scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)

ggplot(scores_df[which(scores_df$Timepoint%in%c(1,2,3,4) | scores_df$Treatment=="t0" ),], aes(x = PC1, y = PC2, group = Treatment,alpha=Timepoint)) +
  geom_point(aes(shape=Treatment,color=Treatment),size = 3) +
  scale_shape_manual(values = c("t0" = 8, "Water" = 15, "Rock" = 16, "Glass" = 17)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values = c("t0"="black","Water" = "blue","Rock"="red","Glass"="orange")) +
  labs(
    x = paste("PC1 (", round(variance_explained[1], 2), "% variance)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2], 2), "% variance)", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16)  
  ) +
  scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)

ggplot(scores_df, aes(x = PC1, y = PC2, group = Treatment,alpha=Timepoint)) +
  geom_point(aes(shape=Treatment,color=Treatment),size = 3) +
  scale_shape_manual(values = c("t0" = 8, "Water" = 15, "Rock" = 16, "Glass" = 17)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_color_manual(values = c("t0"="black","Water" = "blue","Rock"="red","Glass"="orange")) +
  labs(
    x = paste("PC1 (", round(variance_explained[1], 2), "% variance)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2], 2), "% variance)", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16) 
  ) +
  scale_y_continuous(limits=ylim) + scale_x_continuous(limits=xlim)

