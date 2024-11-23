library(dada2)
library(here)
library(tidyverse)

path <- here("data/raw_fastq")
processed_path <- here("data/processed/dada2")

fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(processed_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(processed_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=20)


errF <- learnErrors(filtFs, multithread=20)
errR <- learnErrors(filtRs, multithread=20)

dadaFs <- dada(filtFs, err=errF, multithread=20)
dadaRs <- dada(filtRs, err=errR, multithread=20)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Construct Sequence Table

seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
seqlength <- table(nchar(getSequences(seqtab)))
seqplot <- ggplot(as.data.frame(seqlength), aes(x=Var1, y = Freq)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  labs(x = "sequence length", y = "frequency")
ggsave(seqplot, filename = here("data/processed/dada2/seqplot.png"))


### Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=20, verbose=TRUE)

## Track Reads


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write_csv(track, here("data/processed/dada2/seq_tracking.csv"))

## Assign Taxonomy

taxa <- assignTaxonomy(seqtab.nochim, 
                       here("resources/databases/silva_nr99_v138.2_toGenus_trainset.fa.gz"), 
                       multithread=20)


### Add species

taxa <- addSpecies(taxa, "resources/databases/silva_v138.2_assignSpecies.fa.gz")

## Export data tables

write_csv(seqtab.nochim, "data/processed/dada2/seqtab.csv")
write_csv(taxa, "data/processed/dada2/taxa.csv")



