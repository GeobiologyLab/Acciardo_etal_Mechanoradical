library(ANCOMBC)
library(igraph)
library(vegan)
library(dplyr)
library(phyloseq)
library(compositions)

load("Data/Mechanoradical_Incubation.RData")

taxa$SNV <- rownames(taxa)

### Remove ASVs in control samples

controlSamples <- rownames(df)[grepl("Blank", rownames(df))]
largeVolume <- rownames(df)[grepl("1L", rownames(df))]
df <- df[-which(rownames(df)%in%largeVolume),]
contaminants <- which(colSums(df[controlSamples,])>0)
df.clean <- df[,-contaminants]
df.clean <- df.clean[,-which(colSums(df.clean)==0)]
snvTable <- df.clean[-which(rownames(df.clean)%in%controlSamples),]

writeANCOM2File <- function(myDF,outputfilename){
  o <- myDF[,c("taxon","lfc_Treatmentzrock")]
  o$Kingdom <- taxa[o$taxon,"Kingdom"]
  o$Phylum <- taxa[o$taxon,"Phylum"]
  o$Class <- taxa[o$taxon,"Class"]
  o$Order <- taxa[o$taxon,"Order"]
  o$Family <- taxa[o$taxon,"Family"]
  o$Genus <- taxa[o$taxon,"Genus"]
  write.csv(o,outputfilename)
}


### make the water and glass matrix

gw <- snvTable[grepl('-W1',rownames(snvTable)),]
tmp <- snvTable[grepl('-W2',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-W3',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-W4',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-W5',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-G1',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-G2',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-G3',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-G4',rownames(snvTable)),]
gw <- rbind(gw,tmp)
tmp <- snvTable[grepl('-G5',rownames(snvTable)),]
gw <- rbind(gw,tmp)

### make the rock, water, and glass table
tmp <- snvTable[grepl('-TAn1',rownames(snvTable)),]
wgr <- rbind(gw,tmp)
tmp <- snvTable[grepl('-TAn2',rownames(snvTable)),]
wgr <- rbind(wgr,tmp)
tmp <- snvTable[grepl('-TAn3',rownames(snvTable)),]
wgr <- rbind(wgr,tmp)
tmp <- snvTable[grepl('-TAn4',rownames(snvTable)),]
wgr <- rbind(wgr,tmp)
tmp <- snvTable[grepl('-TAn5',rownames(snvTable)),]
wgr <- rbind(wgr,tmp)

sample_ids <- rownames(wgr)

short_ids <- gsub(".*-([A-Za-z0-9]+_[0-9]+)_.*", "\\1", sample_ids)

rownames(wgr) <- gsub(".*-([A-Za-z0-9]+_[0-9]+)_.*", "\\1", rownames(wgr))

# Create the metadata object
metadata <- data.frame(
  SampleID = short_ids,                           # Sample IDs
  Treatment = c(rep("water",15), rep("zglass", 15), rep("zrock",15)),  # Treatment group; reference will be the first in alphabetical order
  Time = c(24, 24, 24, 48, 48, 48, 96, 96, 96, 168, 168, 168, 240, 240, 240, 24, 24, 24, 48, 48, 48,96, 96, 96, 168, 168, 168, 240, 240, 240, 24, 24, 24, 48, 48, 48,96, 96, 96, 168, 168, 168, 240, 240, 240),
  row.names = short_ids                           # Set rownames for phyloseq compatibility
)

longmetadata <- metadata
longmetadata <- longmetadata[match(rownames(wgr) , longmetadata$SampleID), ]
rownames(longmetadata) <- longmetadata$SampleID


physeqWGR <- phyloseq(otu_table(wgr, taxa_are_rows = FALSE), tax_table(as.matrix(taxa)),sample_data(longmetadata))


min_prevalence <- 0.1 * nsamples(physeqWGR) 
physeqWGR <- filter_taxa(physeqWGR, function(x) sum(x > 0) >= min_prevalence, prune=TRUE)


outputWGR = ancombc2(data = physeqWGR, tax_level = "SNV",
                     fix_formula = "Time + Treatment", rand_formula = NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE, group = "Treatment",
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, 
                     struc_zero = TRUE, neg_lb = TRUE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                     iter_control = list(tol = 1e-2, max_iter = 20, 
                                         verbose = TRUE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     lme_control = lme4::lmerControl(),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                 nrow = 2, 
                                                                 byrow = TRUE),
                                                          matrix(c(-1, 0, 1, -1),
                                                                 nrow = 2, 
                                                                 byrow = TRUE),
                                                          matrix(c(1, 0, 1, -1),
                                                                 nrow = 2, 
                                                                 byrow = TRUE)),
                                          node = list(2, 2, 1),
                                          solver = "ECOS",
                                          B = 10))

resWGR <- as.data.frame(outputWGR$res)
resWGR <- as.data.frame(resWGR)
rownames(resWGR)<-resWGR$taxon

myOutput <- resWGR
myOutput$Kingdom <- taxa$Kingdom[match(myOutput$taxon, taxa$SNV)]
myOutput$Phylum <- taxa$Phylum[match(myOutput$taxon, taxa$SNV)]
myOutput$Class <- taxa$Class[match(myOutput$taxon, taxa$SNV)]
myOutput$Order <- taxa$Order[match(myOutput$taxon, taxa$SNV)]
myOutput$Family <- taxa$Family[match(myOutput$taxon, taxa$SNV)]
myOutput$Genus <- taxa$Genus[match(myOutput$taxon, taxa$SNV)]


eRock <- myOutput[which(myOutput$diff_Treatmentzrock==TRUE),]
uRock <- eRock[which(eRock$diff_Treatmentzglass==FALSE),]
glassRock <- eRock[which(eRock$diff_Treatmentzglass==TRUE),]

