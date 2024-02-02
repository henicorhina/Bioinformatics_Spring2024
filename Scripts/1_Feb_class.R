library(msa)
library(tidyr)
library(dplyr)
library(genepop)
library(PopGenome)
library(Biostrings)
library(seqinr)
library(phangorn)

setwd("/Users/ojohnson/Documents/GitHub/Bioinformatics_Spring2024")

# read in albatross CytB sequences
seq_1 <- readDNAStringSet("Diomedea_exulans/sequence_1.fasta")
seq_2 <- readDNAStringSet("Diomedea_exulans/sequence_2.fasta")
seq_3 <- readDNAStringSet("Diomedea_exulans/sequence_3.fasta")
seq_4 <- readDNAStringSet("Diomedea_exulans/sequence_4.fasta")
seq_5 <- readDNAStringSet("Diomedea_exulans/sequence_5.fasta")
seq_6 <- readDNAStringSet("Diomedea_exulans/sequence_6.fasta")
seq_7 <- readDNAStringSet("Diomedea_exulans/sequence_7.fasta")
seq_8 <- readDNAStringSet("Diomedea_exulans/sequence_8.fasta")

# combine samples
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6, seq_7, seq_8)

# rename samples
names(seqs) <- c("sanfordi_U48946.1", "chionoptera_AF076048.1", "epomophora_AF076049.1",
                 "gibsoni_AF076050.1", "antipodensis_MH330008.1", "antipodensis_MH330010.1",
                 "exulans_U48947.1", "amsterdamensis_U48948.1")
# old names
# names(seqs) <- c("exulans_1B-111", "exulans_SP613", "exulans_SP622", "exulans_SP491", "exulans_SP616", "exulans_SP641", "exulans_SP623")

# MSA!
albatrossAln <- msa(seqs)

# alignment length
nchar(albatrossAln)
print(albatrossAln, show="complete")

# GC content
alFreq <- alphabetFrequency(albatrossAln)
alFreq
GC <- sum(alFreq[,"C"]) + sum(alFreq[,"G"]) 
AT <- sum(alFreq[,"A"]) + sum(alFreq[,"T"]) 
GC / (GC + AT )

# GC content with seqinr package
seq_1.1 <- read.fasta("Diomedea_exulans/sequence_1.fasta")
seqinr::GC(seq_1.1$U48946.1) # need to select the sequence from the vector


# identity matrix
albatrossAln2 <- msaConvert(albatrossAln, type="seqinr::alignment")
d <- dist.alignment(albatrossAln2, "identity")
as.matrix(d)
round(as.matrix(d)[1:8, "epomophora_AF076049.1", drop=FALSE], digits = 2) * 100

# translate
seq_1_AA <- Biostrings::translate(seq_1, no.init.codon=TRUE)
print(seq_1_AA)


# write alignment (harder than I expected)
albatrossAln_phyDat <- msaConvert(albatrossAln, type="phangorn::phyDat")
write.phyDat(albatrossAln_phyDat, "Diomedea_exulans/alignment.fasta", format = "fasta")


# write.fasta(albatrossAln2$seq, albatrossAln2$nam, "Diomedea_exulans/alignment.fasta")
# write.fasta(albatrossAln2, "Diomedea_exulans/alignment.fasta")
# 
# DNAStr <- as(albatrossAln, "AlbatrossDNAStringSet")
# writeXStringSet(DNAStr, file="myFile.fa")


###############################################################################
########################## example sequences ##################################
###############################################################################

# read in example sequence file
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences


# align using defaults
myFirstAlignment <- msa(mySequences)
myFirstAlignment

# convert to seqinr and compute distances
myFirstAln2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
d <- dist.alignment(myFirstAln2, "identity")
as.matrix(d)
round(as.matrix(d)[1:9, "PH4H_Homo_sapiens", drop=FALSE], digits = 2) * 100



###############################################################################
########################## Buffalo sequences ##################################
###############################################################################

buff_seq <- readDNAStringSet("egfr_flank.fasta")
buff_seqAln <- msa(buff_seq)
print(buff_seqAln, show="complete")

buff_seq.1 <- read.fasta("egfr_flank.fasta")
seqinr::GC(buff_seq.1[[1]]) # need to select the sequence from the vector

buff_seqAln2 <- msaConvert(buff_seqAln, type="seqinr::alignment")
buff_seqAln2.d <- dist.alignment(buff_seqAln2, "identity")
as.matrix(buff_seqAln2.d)
round(as.matrix(buff_seqAln2.d)[1:5, "ENSMUSG00000020122|ENSMUST00000125984", drop=FALSE], digits = 2) * 100

buff_seq_AA <- Biostrings::translate(buff_seq[1], no.init.codon=TRUE)
print(buff_seq_AA)

vmatchPattern("TAG", buff_seq)


