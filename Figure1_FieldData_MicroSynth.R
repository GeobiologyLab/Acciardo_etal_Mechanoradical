library(reshape2)
library(dada2)
library(ggplot2)
library(gridExtra)
library(cowplot)



load("Data/dada2_st1_int11_microsynth_annotations_v138.2.RData")

taxa <- data.frame(taxa)
taxa$SNV = rownames(taxa)
controls <- seqtab.nochim[grepl("Blank",rownames(seqtab.nochim)),]

s <- data.frame(seqtab.nochim[-which(rownames(seqtab.nochim)%in%rownames(controls)),-which(colSums(controls)>0)])

s <- s/rowSums(s)#*100

### These are sampling dates: 2023-07-05; 2023-07-27; 2024-02-01; 2024-02-27; 2024-02-01; 2024-02-27; 2024-02-01; 2024-02-27 ie Exp 1 is on 2024-02-01 and Exp2 on 2024-02-27

newID <- c("PreRestim","PostRestim","Exp1_A","Exp2_A","Exp1_B","Exp2_B","Exp1_C","Exp2_C")


rownames(s) <- newID


plottingOrder <- c("PreRestim","PostRestim","Exp1_A","Exp1_C","Exp1_B","Exp2_B","Exp2_A","Exp2_C")

s <- s[match(plottingOrder, rownames(s)),]

dates <- c(as.Date("2023-07-07"),as.Date("2023-07-27"),as.Date("2024-02-01"),as.Date("2024-02-02"),as.Date("2024-02-03"),as.Date("2024-02-27"),as.Date("2024-02-28"),as.Date("2024-02-29"))

splot <- data.frame(ID = rownames(s), Date = dates)

splot$Ca.Desulforudis <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Genus=="Candidatus Desulforudis"),"SNV"])])
splot$Other.Bacillota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Bacillota"),"SNV"])]) - splot$Ca.Desulforudis
splot$Thermodesulfovibrionia <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Thermodesulfovibrionia"),"SNV"])])#rowSums(s[,intersect(uRock[which(uRock$TaxaOfInterest=="Nitrospirota"),"taxon"],taxa[which(taxa$Class=="Thermodesulfovibrionia"),"SNV"])])
splot$Chloroflexota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Chloroflexota"),"SNV"])])
splot$Alphaproteobacteria <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Alphaproteobacteria"),"SNV"])])
splot$Ignavibacteriota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Ignavibacteriota"),"SNV"])])
splot$Syntrophia <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Class=="Syntrophia"),"SNV"])])
splot$Other.Thermodesulfobacteriota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Thermodesulfobacteriota"),"SNV"])]) - splot$Syntrophia
splot$Comamonadaceae <-  rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Family=="Comamonadaceae"),"SNV"])])
splot$Other.Pseudomonadota <- rowSums(s[,which(colnames(s)%in%taxa[which(taxa$Phylum=="Pseudomonadota"),"SNV"])]) - splot$Comamonadaceae -splot$Alphaproteobacteria
splot$Other <- rowSums(s) - rowSums(splot[,3:ncol(splot)])


s.melt <- melt(splot,id.vars = c("ID","Date"))
s.melt$ID <- factor(s.melt$ID, levels = plottingOrder)

a <- ggplot(s.melt[which(s.melt$ID%in%c("PreRestim","PostRestim")),], aes(x = Date, y = value, fill = variable)) +
  geom_area(alpha=0.5) +
  geom_bar(stat = "identity", position = "stack", width = 1,color="black",size=0.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(y = "Relative Abundance") +
  ggtitle("ST1-int11 (2023); Induced Seismicity") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",axis.title.x=element_blank())


b <- ggplot(s.melt[which(s.melt$ID%in%c("Exp1_A","Exp1_B","Exp1_C","Exp2_A","Exp2_B","Exp2_C")),], aes(x = Date, y = value, fill = variable)) +
  geom_area(alpha=0.5) +
  geom_bar(stat = "identity", position = "stack", width = 1,color="black",size=0.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(y = "-- Additional  Seismicity  2023 to 2024 --") +
  ggtitle("ST1-int11 (2024); No Seismicity") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",axis.title.x=element_blank())

legend_plot <- ggplot(s.melt[which(s.melt$ID %in% c("PreRestim", "PostRestim2")),], aes(x = Date, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  theme(legend.position = "top", legend.box = "horizontal",legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 2)) 

# Function to extract the legend manually
extract_legend <- function(plot) {
  g <- ggplotGrob(plot) # Convert plot to grob
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")]
  if (length(legend) == 1) {
    return(legend[[1]]) # Return the legend grob
  } else {
    stop("Legend not found or multiple legends exist.")
  }
}

# Extract the legend
legend <- extract_legend(legend_plot)


# Arrange plots and legend
plot_grid(

  a,      
  b,       
  legend, a,
  ncol = 2, 
  rel_heights = c(1, 1),  
  rel_widths = c(1, 1)    
)


