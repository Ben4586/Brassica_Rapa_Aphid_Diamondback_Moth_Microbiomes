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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 220), #trunclen will remove seqs shorter than provided value
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
mapfile <- Plant_metadata.csv" #make a csv file with sample metadata included in data folder. 
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
ps16 #filtered out 991 ASVs
sum(otu_table(ps16)) #number of seqs remaining
saveRDS(ps16, "ps16.rds")

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
ps_filtered <- subset_taxa(ps_rarefied, !is.na(Phylum) & !Phylum %in% c("Nitrospirota", "Cyanobacteria", "Fibrobacterota"))

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

# Define prevalence threshold as 5% of total samples #remove ASVs present in <5% of samples
prevalenceThreshold = 0.05 * nsamples(ps_filtered)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps_filtered)
saveRDS(ps1, "ps1.rds")

#Visualization/Diversity
classtaxa <- get_taxadf(obj=ps1, taxlevel=2) #taxlevel can be adjusted as needed (genus or phylum)
# The 10 most abundant taxonomy will be visualized by default (parameter `topn=10`). 
pclass <- ggbartax(obj=classtaxa, facetNames="Treatment", topn = 10) +
  xlab(NULL) +
  ylab("Relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
pclass

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
ps1_genus_relabun <- transform_sample_counts(ps_genus, function(OTU) OTU/sum(OTU) * 100)
ps1_genus_top10 <- prune_taxa(top10_genera, ps1_genus_relabun)
phyloseq::taxa_names(ps1_genus_top10) <- phyloseq::tax_table(ps1_genus_top10)[, "Genus"]
phyloseq::otu_table(ps1_genus_top10)[1:5, 1:5]

#Melt and plot boxplot at genus level
phyloseq::psmelt(ps1_genus_top10) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Relative abundance (%)\n") +
  facet_wrap(~ OTU, scales = "free")

#Agglomerate to the phylum level 
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)

# Merge sequence object into the phyloseq object
ps_raw <- merge_phyloseq(ps1, dna)
ps_raw

taxa_names(ps_raw) <- paste("ASV", 1:ntaxa(ps_raw), sep = "_")
taxa_names(ps_raw)[1:3]

#Convert to relative abundance
ps_phylum <- phyloseq::tax_glom(ps_raw, "Phylum")
ps1_phylum_relabun <- transform_sample_counts(ps_phylum, function(OTU) OTU/sum(OTU) * 100)

#Melt and plot boxplot at phylum level
phyloseq::psmelt(ps1_phylum_relabun) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Relative abundance (%)\n") +
  facet_wrap(~ OTU, scales = "free")

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
plot_ordination(ps.prop, ord.nmds.bray, color="Treatment") +
  theme_classic() +
  theme(strip.background = element_blank())+
  stat_ellipse()

# Plot using plot_ordination function of phyloseq
bray.pcoa <- ordinate(ps.prop, method = "PCoA", distance = "bray")
plot_ordination(ps.prop, bray.pcoa, color = "Treatment", axes = c(1,2)) +
  geom_point(size = 2) +
  labs(title = "PCoA of Bray Curtis Distances", color = "Treatment") +
  theme_classic() +
  theme(strip.background = element_blank())+
  stat_ellipse()



