library(dada2)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggtext)


md <- read.csv("Data/IncubationResults_Metadata.csv", sep=";")
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

s <- snvTable
rownames(s)<- gsub(".*-([A-Za-z0-9]+_[0-9]+)_.*", "\\1", rownames(s))

myOrder <- c("t0_1","t0_2","t0_3","W1_1","W1_2","W1_3","W2_1","W2_2","W2_3","W3_1","W3_2","W3_3","W4_1","W4_2","W4_3",
             "W5_1","W5_2","W5_3","G1_1","G1_2","G1_3","G2_1","G2_2","G2_3","G3_1","G3_2","G3_3","G4_1","G4_2","G4_3",
             "G5_1","G5_2","G5_3","TAn1_1","TAn1_2","TAn1_3","TAn2_1","TAn2_2","TAn2_3","TAn3_1","TAn3_2","TAn3_3",
             "TAn4_1","TAn4_2","TAn4_3","TAn5_1","TAn5_2","TAn5_3")

s <- s[match(myOrder, rownames(s)), ]

s <- s/rowSums(s)

## make table of select taxa

splot <- data.frame(ID = rownames(s))

splot["Ca. Desulforudis"] <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Genus=="Candidatus Desulforudis"),"SNV"])])
splot["Other Bacillota"] <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Bacillota"),"SNV"])]) - splot["Ca. Desulforudis"]
splot$Thermodesulfovibrionia <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Thermodesulfovibrionia"),"SNV"])])#rowSums(s[,intersect(uRock[which(uRock$TaxaOfInterest=="Nitrospirota"),"taxon"],taxa[which(taxa$Class=="Thermodesulfovibrionia"),"SNV"])])
splot$Chloroflexota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Chloroflexota"),"SNV"])])
splot$Ignavibacteriota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Ignavibacteriota"),"SNV"])])
splot$Alphaproteobacteria <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Alphaproteobacteria"),"SNV"])])
splot$Syntrophia <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Syntrophia"),"SNV"])])
splot$Thermodesulfobacteriota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Thermodesulfobacteriota"),"SNV"])]) - splot$Syntrophia
uc <- taxa[which(taxa$Family=="Comamonadaceae"),]
splot["Unclassified Comamonadaceae"] <- rowSums(s[,which(colnames(s)%in%uc[is.na(uc$Genus),"SNV"])])
splot$Hydrogenophaga <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Genus=="Hydrogenophaga"),"SNV"])])
splot$Other <- rowSums(s) - rowSums(splot[,2:ncol(splot)])


## prep for plotting

s.melt <- melt(splot,id.vars = "ID")
s.melt$ID <- factor(s.melt$ID, levels = myOrder)

s.melt$abundance <- s.melt$value * md$qPCR_Count[match(s.melt$ID, md$Sample)]

renameZero <- function(string){
  print(dim(s.melt))
  tmp <- s.melt[grepl("t",s.melt$ID),]
  tmp$ID <- gsub("^t", string, tmp$ID)
  #print(tmp)
  return(bind_rows(s.melt,tmp))
}

s.melt<-renameZero("W")
s.melt<-renameZero("TAn")
s.melt<-renameZero("G")

myReorder <- read.csv("Data/barplot_mapping_file.csv", sep=",")
s.melt$key <- myReorder$new[match(s.melt$ID,myReorder$orig)]

s.melt$normabundance <- s.melt$abundance/10000

line_data <- s.melt %>%
  group_by(key) %>%
  summarise(total_abundance = sum(abundance))

line_data$total_abundance <- log2(line_data$total_abundance)



scale_factor <- 20 # for log2
 



a <- ggplot(s.melt[grepl('W', s.melt$ID), ], aes(x = key, y = value * scale_factor, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_line(data = line_data[grepl("W", line_data$key), ], inherit.aes = FALSE,
            aes(x = key, y = total_abundance, group = 1),
            color = "black", size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.text.y.right=element_blank(),
        axis.title.y = element_markdown(color = "black", size = 14),
        axis.text.y.left = element_markdown(color="black",size=12)) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(
    name =  "**---** log<sub>2</sub>(16S rRNA gene copies mL<sup>-1</sup>)"  ,       
    sec.axis = sec_axis(~ . / scale_factor)) 

b <- ggplot(s.melt[grepl('TAn', s.melt$ID), ], aes(x = key, y = value * scale_factor, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_line(data = line_data[grepl("TAn", line_data$key), ], inherit.aes = FALSE,
            aes(x = key, y = total_abundance, group = 1),
            color = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Paired") +
  labs(y = " ", x = "") +
  scale_y_continuous(
    name = " ",
    sec.axis = sec_axis(~ . / scale_factor, name = " ")
  )

c <- ggplot(s.melt[grepl('G', s.melt$ID), ], aes(x = key, y = value * scale_factor, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_line(data = line_data[grepl("G", line_data$key), ], inherit.aes = FALSE,
            aes(x = key, y = total_abundance, group = 1),
            color = "black", size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y.left=element_blank(), axis.text.y.left=element_blank(),
        axis.title.y.right=element_markdown(color = "black", size = 14),
        axis.text.y.right=element_markdown(color = "black", size = 12)) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "")
  )

# Extract a legend from one of the plots
legend <- get_legend(
  a + labs(fill = "Relative Abundance\nof Assigned Taxa") + theme(legend.position = "right")
)

# Remove legends from individual plots and arrange them in a grid
a <- a + theme(legend.position = "none")
b <- b + theme(legend.position = "none")
c <- c + theme(legend.position = "none")

grid.arrange(a, b, c, legend, ncol = 4)

