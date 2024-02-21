# load in all the libraries that you might need
# this should always be at the start of your script
library(msa)
library(tidyr)
library(dplyr)
library(genepop)
library(PopGenome)
library(Biostrings)
library(seqinr)
library(phangorn)

# set the working directory to the folder containing all of your scripts and data
# filepaths and files should always be in quotes. Variables in R should not.
setwd("/Users/ojohnson/Documents/GitHub/Bioinformatics_Spring2024")
setwd("/Users/ojohnson/Documents/Classes/2024S_Bioinformatics/Exams/Midterm_1_data")

seq_1 <- readDNAStringSet("sequences.fasta")
seqs <- seq_1
albatrossAln <- msa(seqs)
print(albatrossAln, show="complete")
albatrossAln2 <- msaConvert(albatrossAln, type="seqinr::alignment")
d <- dist.alignment(albatrossAln2, "similarity")
d
100 - (round(as.matrix(d)[, "Homo_sapiens_6", drop=FALSE], digits = 2) * 100)

seq_6 <- seqs$Homo_sapiens_6
seq_6_AA <- Biostrings::translate(seq_6)
print(seq_6_AA, show="complete")

output_file <- "amino_acid_sequence.fasta"
writeXStringSet(seq_6_AA, file = output_file,
                format = "fasta")
write.fasta(as.character(seq_6_AA), "seq_6_AA", output_file)

blast_test <- blast_nucleotide_to_nucleotide(
  query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
  subject = system.file('seqs/sbj_nn.fa', package = 'metablastr'),
  output.path = tempdir(),
  db.import  = FALSE)


AccessionNumbersGO <- GetProteinGOInfo(c("A0A0J9YWK4"))
GP <- GetPathology_Biotech(c("A0A0J9YWK4"))
GP <- GetPathology_Biotech(c("P68871"))



#### Read this file into R using the appropriate function ####
accession_numbers<- read.table("/Users/ojohnson/Downloads/AccNumbers.txt")

write.fasta()



# read in albatross Cytochrome B sequences 
# assign each one to a variable
# Note that these fasta files are contained in a folder called 'Diomedea_exulans'
seq_1 <- readDNAStringSet("sequence1.fasta")
seq_2 <- readDNAStringSet("sequence2.fasta")
seq_3 <- readDNAStringSet("sequence3.fasta")
seq_4 <- readDNAStringSet("sequence4.fasta")
seq_5 <- readDNAStringSet("sequence5.fasta")
seq_6 <- readDNAStringSet("sequence6.fasta")
seq_7 <- readDNAStringSet("sequence7.fasta")
seq_8 <- readDNAStringSet("sequence8.fasta")

# combine samples into a single variable using the combine ('c') function
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6, seq_7, seq_8)
seqs <- seq_1
# The default GenBank names are very very long, so lets
# rename the samples to something shorter and more legible
# we do this by assigning a list of characters (using the same 'c' function)
# to the 'names' of the seqs variable
names(seqs) <- c("sanfordi_U48946.1", "chionoptera_AF076048.1", "epomophora_AF076049.1",
                 "gibsoni_AF076050.1", "antipodensis_MH330008.1", "antipodensis_MH330010.1",
                 "exulans_U48947.1", "amsterdamensis_U48948.1")

# run the MSA! Assign it to a new variable
albatrossAln <- msa(seqs)

# check the alignment length, two different ways
nchar(albatrossAln)
print(albatrossAln, show="complete")

# Calculate the GC content
# first, calculate the alphabet frequency
alFreq <- alphabetFrequency(albatrossAln)
alFreq
# now pull out the total number of G's and C's
# here, the 'sum' function takes the sum, as you might expect
# the square brackets are for accessing rows and columns of a matrix
# values before the comma access rows, those after the comma access columns
# we want the columns
GC <- sum(alFreq[,"C"]) + sum(alFreq[,"G"]) 
AT <- sum(alFreq[,"A"]) + sum(alFreq[,"T"]) 
# and calculate the percentage that are G or C (out of the total)
GC / (GC + AT )

# calculate the GC content a different way, using the 'GC' function in the seqinr package
# we can only run this on one sample at a time, so let's read in the first sample using the 
# read.fasta() function in the seqinr package
# note that because there are multiple 'read.fasta()' functions, we need to specify that
# want to use the one in the seqinr package using the double colon
seq_1.seqinr <- seqinr::read.fasta("sequence1.fasta")
# then select the sequence data from the variable. Use the '$' to access it 
seqinr::GC(seq_1.seqinr$Homo_sapiens_1)

# calculate the identity matrix
# first, convert the alignment to the seqinr format using msaConvert
# because the dist.alignment() function is part of the seqinr package
albatrossAln2 <- msaConvert(albatrossAln, type="seqinr::alignment")
d <- dist.alignment(albatrossAln2, "similarity")
d
# this is a fancy way to compare only my 'epomophora' sample to the other samples in the matrix
# and convert the numbers to a percentage 
100 - (round(as.matrix(d)[, "epomophora_AF076049.1", drop=FALSE], digits = 2) * 100)

# translate one sample to amino acids
# we again need to specify which package to use because the 'translate' function 
# exists in both Biostrings and seqinr
seq_1_AA <- Biostrings::translate(seq_1)
print(seq_1_AA)


# write the alignment to a fasta file (harder than I expected to figure this out) 
albatrossAln_phyDat <- msaConvert(albatrossAln, type="phangorn::phyDat")
write.phyDat(albatrossAln_phyDat, "Diomedea_exulans/alignment.fasta", format = "fasta")




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


