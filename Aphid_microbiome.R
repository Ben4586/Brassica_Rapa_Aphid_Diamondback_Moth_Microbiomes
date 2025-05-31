##Code for DADA2 processing adapted from https://benjjneb.github.io/dada2/tutorial.html

# The following R packages are needed to run the analysis
library(dada2)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(vegan)
library(phangorn)
library(DECIPHER)
library(stringr)
library(MicrobiotaProcess)
library(patchwork)
library(VennDiagram)
library(DESeq2)
library(ggbeeswarm)
library(ggrepel)
library(tidyverse)

#Remove adapter, barcode and primer sequences in raw reads (This step has been performed by the sequencing company upon receiving the data)
path <- "data/clean_fastq" # CHANGE ME to the directory containing the clean fastq files.
head(list.files(path))

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)) #change file extension pattern if needed
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

#Plot Sequence Quality
# Extract sample names, assuming filenames have format: NUMBER_SAMPLENAME_XXX.fastq <- file names must be in this format
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2]) 
plotQualityProfile(fnRs[1:2])

#Filter and trim sequences
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 200), #trunclen will remove seqs shorter than provided value
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                     compress=TRUE, multithread=FALSE, matchIDs=TRUE, verbose=TRUE) # On Windows set multithread=FALSE
head(out)
colSums(out) 

#Estimate errors
errF <- learnErrors(filtFs, multithread=FALSE, randomize = T) #keep track of how many steps required for convergence
plotErrors(errF, nominalQ=TRUE)
errR <- learnErrors(filtRs, multithread=FALSE, randomize = T)
plotErrors(errR, nominalQ=TRUE)

#Dereplicate (unique) samples
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Create dada object
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
dadaFs[[1]]

#Merge read pairs and make seq table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #gives number of ASVs and samples 
sum(seqtab) #gives total number of seqs 

#Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab))) 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(370,400)] #remove seqs that are too long to be correct.  This parameter can be adjusted as needed. 
table(nchar(getSequences(seqtab2)))
dim(seqtab2) #gives number of ASVs and samples
sum(seqtab2) #gives total number of seqs 

#Find and remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim) #number of samples and ASVs after removing chimeras
sum(seqtab.nochim) #number of seqs after removing chimeras 
sum(seqtab.nochim)/sum(seqtab2) #percent seqs kept after removing chimeras

#Summary of filtered seqs
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy to ASVs
taxa <- assignTaxonomy(seqtab.nochim, "data/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE, tryRC=TRUE) #tryRC = try reverse complement. 
saveRDS(taxa, "C:/Users/bengsoon/Desktop/Plant_microbiome/Plant_microbiome_day21_root/taxa.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Make phyloseq object
mapfile <- Aphid_metadata.csv" #make a csv file with sample metadata included in data folder. 
map <- read.csv(mapfile)
map <- sample_data(map)
colnames(map)[1] <- "SampleID"
rownames(map) <- map$SampleID

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa)) #makes phyloseq object, which contains ASV count table, sample data and taxonomy information.

#Save phyloseq objects and sample files
#Save dada2/phyloseq objects in case you need them again
saveRDS(ps, "ps.rds")
saveRDS(taxa, "taxa.rds")
saveRDS(map, "map.rds")
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#Filter out contaminants
ps16 <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Family != "Mitochondria" &
      Class   != "Chloroplast"
  )

ps
ps16 #filtered out 646 ASVs
sum(otu_table(ps16)) #number of seqs remaining
saveRDS(ps16, "ps16.rds")

##Code for Buchnera filtering processing adapted from https://github.itap.purdue.edu/LaramyEndersGroup/Milkweed-Aphid-Microbiome/blob/master/github%20upload/code/MilkweedAphid_preprocessing_08182021revised.R
#Count OTU abundance
taxa_sum_df <- data.frame(sum = taxa_sums(ps16))
head(taxa_sum_df, 15)
smin <- min(taxa_sums(ps16))
smean <- mean(taxa_sums(ps16))
smax <- max(taxa_sums(ps16))
smin
smean
smax

#Retain taxa found more than 10 times in at least 5% of the samples
ps16.filter <- filter_taxa(ps16, function(x) sum(x > 10) > (0.05*length(x)), TRUE)
ps16.filter
ps16

#Filter out samples with low counts.
#Histogram of sample sizes
sample_sum_df <- data.frame(sum = sample_sums(ps16.filter))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(ps16.filter))
smean <- mean(sample_sums(ps16.filter))
smax <- max(sample_sums(ps16.filter))

smin
smean
smax

#Filter samples with fewer than 2500 reads.
ps16.trim <- prune_samples(sample_sums(ps16.filter)>=2500, ps16.filter)
sum(otu_table(ps16.filter))
sum(otu_table(ps16.trim))
saveRDS(ps16.trim, "C:/Users/bengsoon/Desktop/Aphid_microbiome_second/ps16.trim.rds")
sample_sum_df2 <- data.frame(sum = sample_sums(ps16.trim))

# Histogram of sample read counts
ggplot(sample_sum_df2, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#Sample min/mean/max after removing low/high count samples
smin <- min(sample_sums(ps16.trim))
smean <- mean(sample_sums(ps16.trim))
smax <- max(sample_sums(ps16.trim))

smin
smean
smax

#Remove samples with low count
all.samples <- as.vector(sample_data(ps16.filter)[[1]])
nolow.samples <- as.vector(sample_data(ps16.trim)[[1]])
subs <- all.samples %in% nolow.samples
subset(all.samples, subset = !subs) #list of samples that were removed due to low/high counts

#Save filtered, trimmed phyloseq object
saveRDS(ps16.trim, "ps16.trim.rds")
sum(otu_table(ps16.trim))

#Filter samples for cross contamination
#Additional package needed for trimming and filtering
library(metagMisc)

# read phyloseq object
ps16.trim <- readRDS("ps16.trim.rds")
#taxa_names(ps16.trim) <- paste0("ASV", sep = "_", seq(ntaxa(ps16.trim)))#replace sequences with ASV names

# Pre and post filtered tables by identifying read counts <50 and replace with 0
pre.filt2.data <- otu_table(ps16.trim)
post.filt2.data <- pre.filt2.data
row.means <- data.frame(matrix(NA, nrow(pre.filt2.data),2))
colnames(row.means) <- c("asv.id", "mean")
row.means$asv.id <- pre.filt2.data[,1]
for(i in 1:nrow(pre.filt2.data)){
  for(j in 2:ncol(pre.filt2.data)){
    if(pre.filt2.data[i,j] <= 50){
      post.filt2.data[i,j] <- 0
    }
  }
  row.means$mean[i] <- mean(as.numeric(post.filt2.data[i,2:ncol(post.filt2.data)]))
  print(paste(post.filt2.data[i,1], sprintf("has a mean of %f", mean(as.numeric(post.filt2.data[i,2:ncol(post.filt2.data)]))), sep = " "))
}


#Create a new phyloseq object using the filtered/transformed data
physampledata <- sample_data(ps16.trim)

random_tree = rtree(ntaxa(ps16.trim), rooted=TRUE, tip.label=taxa_names(ps16.trim))
plot(random_tree)

phytaxtable <- tax_table(ps16.trim)

physeq1 = merge_phyloseq(ps16.trim, physampledata, phytaxtable, random_tree)
physeq1

ps16.trim.r <- phyloseq(otu_table(post.filt2.data, taxa_are_rows=FALSE),
                        sample_data(physampledata),
                        tax_table(phytaxtable),
                        phy_tree(physeq1))

#Remove taxa with samples containing reads less than 50 reads threshold (all zeros in OTU table)
ps16.trim.r <- filter_taxa(ps16.trim.r, function(x) sum(x) > 0, TRUE)
ps16.trim
ps16.trim.r #filtered down to 50 total ASVs
#taxa_sums(ps16.trim.r) #use only if changed taxa names to ASV

#Save revised, filtered, trimmed phyloseq object
saveRDS(ps16.trim.r, "ps16.trim.r.rds")
sum(otu_table(ps16.trim))#started with sequences 3443874
sum(otu_table(ps16.trim.r))#remaining sequences 3441573
ps16.trim.r #filtered out no taxa

#Save OTU table of raw counts for revised filtered, trimmed dataset
taxa_names(ps16.trim.r) <- paste0("ASV", sep = "_", seq(ntaxa(ps16.trim.r)))#rename sequences to ASVs
write.csv(t(otu_table(ps16.trim.r)), file = "raw_counts_mw.r.csv")#csv file of OTU table

#Filtering Buchnera cross-contamination
# load filtered, trimmed phyloseq object
ps16.trim <- readRDS("ps16.trim.rds")

#Filter for Buchnera ASVs using 1% read cutoff
pre.Buch.filt.data <- otu_table(ps16.trim)
post.Buch.filt.data <- pre.Buch.filt.data

#Align full set of sequences with subset of Buchnera only sequences
#Search for columns (ASV) which are Buchnera
sample.sum <- rep(NA, nrow(post.Buch.filt.data))
find.Buchnera <- taxa_names(ps16.trim) %in% taxa_names(subset_taxa(ps16.trim, Genus=="Buchnera"))
for(i in 1:nrow(pre.Buch.filt.data)){
  sample.sum[i] <- sum(pre.Buch.filt.data[i,])
  for(j in which(find.Buchnera == TRUE)){
    print(sprintf("Sample %i, Buchnera OTU %i: threshold %f", i,j, (sample.sum[i]*0.01)))
    if(pre.Buch.filt.data[i,j] < (sample.sum[i]*0.01)){
      post.Buch.filt.data[i,j] <- 0
    }
  }
}

#Create a new phyloseq object using  the filtered/transformed data 
ps16.trim.Buchnofilt <- phyloseq(otu_table(pre.Buch.filt.data, taxa_are_rows=FALSE),
                                 sample_data(ps16.trim),
                                 tax_table(ps16.trim),
                                 phy_tree(physeq1))
ps16.trim.Buchfilt <- phyloseq(otu_table(post.Buch.filt.data, taxa_are_rows=FALSE),
                               sample_data(ps16.trim),
                               tax_table(ps16.trim),
                               phy_tree(physeq1))

#Save ASV table of raw counts for revised filtered, trimmed dataset
taxa_names(ps16.trim.Buchnofilt) <- paste0("ASV", sep = "_", seq(ntaxa(ps16.trim.Buchnofilt)))#rename sequences to ASVs
taxa_names(ps16.trim.Buchfilt) <- paste0("ASV", sep = "_", seq(ntaxa(ps16.trim.Buchfilt)))#rename sequences to ASVs
write.csv(t(otu_table(ps16.trim.Buchnofilt)), file = "raw_counts_preBuchfilt.csv")#csv file of OTU table
write.csv(t(otu_table(ps16.trim.Buchfilt)), file = "raw_counts_postBuchfilt.csv")#csv file of OTU table

#Discard taxa with samples containing no reads after filtering (all zeros)
ps16.trim.Buchfilt <- filter_taxa(ps16.trim.Buchfilt, function(x) sum(x) > 0, TRUE)
ps16.trim
ps16.trim.Buchfilt 

#Save revised Buchnera, filtered, trimmed phyloseq object
saveRDS(ps16.trim.Buchfilt, "ps16.trim.Buchfilt.rds")
sum(otu_table(ps16.trim))#started with sequences 3443874
sum(otu_table(ps16.trim.Buchfilt))#remaining sequences 3394287
ps16.trim.Buchfilt #filtered out Buchnera

#Rarefaction to normalize the sequencing depth for all samples
ps_rarefied <- rarefy_even_depth(ps16,
                                 rngseed = 1,
                                 sample.size = min(sample_sums(ps16)),
                                 replace = FALSE)

# Have a look at the result
otu_table(ps_rarefied)[1:5, 1:5]

#Compare sample reads before and after rarefaction
head(sample_sums(ps_rarefied))
head(sample_sums(ps16))

#Show ranks in the dataset
rank_names(ps_rarefied)

# Create table, number of features for each phyla
table(tax_table(ps_rarefied)[, "Phylum"], exclude = NULL)

# Filter out low abundance phyla
ps_filtered <- subset_taxa(ps_rarefied, !is.na(Phylum) & !Phylum %in% c("Abditibacteriota", "Bdellovibrionota", "Entotheonellaeota", "Hydrogenedentes", "Nitrospirota", "Calditrichota", "Deferrisomatota", "Fibrobacterota", "Latescibacterota", "Campylobacterota", "Deinococcota", "Methylomirabilota", "Armatimonadota", "Dependentiae", "Fusobacteriota", "Modulibacteria", "Spirochaetota", "Cloacimonadota", "Sumerlaeota"))

#Prevalence filtering
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_filtered),
               MARGIN = ifelse(taxa_are_rows(ps_filtered), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_filtered),
                    tax_table(ps_filtered))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps_filtered , "Phylum"))

# Define prevalence threshold as 5% of total samples (Remove ASVs present in <5% of samples)
prevalenceThreshold = 0.05 * nsamples(ps_filtered)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps_filtered)
saveRDS(ps1, "ps1.rds")

#Visualization of genus relative abundance
# How many genera would be present after filtering? 
length(get_taxa_unique(ps1, taxonomic.rank = "Genus"))

ps1_genus = tax_glom(ps1, "Genus", NArm = TRUE)

# Get top 10 genera
top10_genera <- names(sort(taxa_sums(ps1_genus), decreasing=TRUE))[1:10]

# Transform Taxa counts to relative abundance
ps1_genus_relabun <- transform_sample_counts(ps1_genus, function(OTU) OTU/sum(OTU) * 100)

# Extract the top 10 taxa and Regular Treatment Samples
ps1_genus_top10 <- prune_taxa(top10_genera, ps1_genus_relabun)

# Convert into dataframe
taxa_abundance_table_genus <- psmelt(ps1_genus_top10)

StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~ Time, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

StackedBarPlot_genus

#Taxonomy compositional analysis
#Filter out the genus Buchnera of the family Morganellaceae in search for other less abundant genera
nobuch=subset_taxa(ps16.trim.Buchfilt, Family!="Morganellaceae")
nobuch

#Visualization/Diversity
classtaxa <- get_taxadf(obj=nobuch, taxlevel=6)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass <- ggbartax(obj=classtaxa, facetNames="Time", topn = 10) + 
  xlab(NULL) +
  ylab("Relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
pclass

#Plot venn diagram
vennlist <- get_vennlist(obj=ps1, factorNames="Treatment")
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "#FD9347"),
                      cat.col=c("#00AED7", "#FD9347"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)

#Alpha diversity
#Plot richness boxplot
divIdx = c("Observed", "Chao1", "Shannon", "Simpson")
ad = estimate_richness(ps1, measures = divIdx)
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")
ad = ad %>% select(SampleID, Treatment, all_of(divIdx)) %>% 
  gather(key = "alpha", value = "measure", -c(SampleID, Treatment))

ggplot(ad, aes(Treatment, measure, color = Treatment)) +
  geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
  geom_quasirandom(size = 0.8, color = "black") + theme_bw() + 
  facet_wrap(~ alpha, scales = "free_y", nrow = 1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Treatment", y = "Alpha Diversity Measure")

#Beta diversity
#PERMANOVA
distme <- get_dist(ps1, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(ps1), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Group <- factor(sampleda$Treatment)
set.seed(1024)
adores <- adonis2(distme ~ Treatment, data=sampleda, permutation=9999)

#ANOSIM
dist = phyloseq::distance(ps1, method="bray", binary = TRUE) 
metadata <- data.frame(sample_data(ps1))
anosim(dist, metadata$Treatment)

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps1, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
myplot <- plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", shape = "Time")+stat_ellipse()
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Plot using plot_ordination function of phyloseq
dist = phyloseq::distance(ps1, method="bray")
ordination = ordinate(ps1, method="PCoA", distance=dist)
plot_ordination(ps1, ordination, color="Treatment", shape = "Time") + 
  theme_classic() +
  theme(strip.background = element_blank()) +
  stat_ellipse()






#Visualization of relative abundance
#Agglomerate to the genus level 
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)

# Merge sequence object into the phyloseq object
ps_raw <- merge_phyloseq(ps1, dna)
ps_raw

taxa_names(ps_raw) <- paste("ASV", 1:ntaxa(ps_raw), sep = "_")
taxa_names(ps_raw)[1:3]

# Transform Taxa counts to relative abundance
ps1_genus = tax_glom(ps_raw, "Genus", NArm = TRUE)

# Get top 10 genera
top10_genera <- names(sort(taxa_sums(ps1_genus), decreasing=TRUE))[1:10]

ps_genus <- phyloseq::tax_glom(ps1_genus, "Genus")
ps1_genus_relabun <- transform_sample_counts(ps1_genus, function(OTU) OTU/sum(OTU) * 100)
ps1_genus_top10 <- prune_taxa(top10_genera, ps1_genus_relabun)
phyloseq::taxa_names(ps1_genus_top10) <- phyloseq::tax_table(ps1_genus_top10)[, "Genus"]
phyloseq::otu_table(ps1_genus_top10)[1:5, 1:5]

# Convert into dataframe
taxa_abundance_table_genus <- psmelt(ps1_genus_top10)

#Melt and plot boxplot at genus level
phyloseq::psmelt(ps1_genus_top10) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Relative abundance\n") +
  facet_wrap(~ OTU, scales = "free")
