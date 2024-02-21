# install.packages("protti", dependencies = TRUE)
# install.packages("r3dmol")
# install.packages("UniprotR")
# BiocManager::install("GenomicAlignments")
library(bio3d)
library(protti)
library(UniprotR)
library(tidyverse)
library(r3dmol)
library(phangorn)

# select one sample from last HW and translate to amino acids
seq_1 <- readDNAStringSet("Diomedea_exulans/sequence_1.fasta")
seq_1_AA <- Biostrings::translate(seq_1)
print(seq_1_AA, show="complete")
# write to fasta file
writeXStringSet(seq_1_AA, "Diomedea_exulans/sequence_1_AA.fasta")

# read in results from UniProt BLAST:
Accession.res <- read.csv("Diomedea_exulans/accessions.csv")
# format to character string by accessing the column of the dataframe containing the accession numbers
Accessions <- Accession.res$Accession

# example accession numbers to use if yours don't work:
# Accessions <- c("P0A799", "P08839")

# Or get some accession numbers from an online repository
# Accessions <- GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv") 

# get GO info for all 5 accession numbers
GeneOntologyObj <- GetProteinGOInfo(Accessions) 

# plot, save to file in new folder
PlotGoInfo(GeneOntologyObj, directorypath = "Diomedea_exulans_GO")
# nicer plot:
PlotGOAll(GOObj = GeneOntologyObj, Top = 10, directorypath = "Diomedea_exulans_GO", width = 8, height = 5)

# Interesting GO terms:
# all related to energy production and mitochondrial function
# e.g.: metal ion binding, respiratory chain complex, and inner membrane of the mitochondrion

# function to pull out the names of the taxon with the best match to our accession numbers
organism.names <- GetNamesTaxa(Accessions)

# get pathologies and diseases:
GPB <- GetPathology_Biotech(Accessions) 
diseases <- Get.diseases(GPB) # note that this needs the results of the GPB, so if you got no results then there will be nothing to run

# get UniProt metadata for your Accessions
uniprot_info <- fetch_uniprot(Accessions)

# this is where the 3D model ID numbers will be if they exist
xref_pdb <- uniprot_info$xref_pdb
xref_pdb_split <- str_split(xref_pdb, ";") # divide string by semicolons

# run just the first sample to get pdb file
ptsi_pgk_pdb_info <- fetch_pdb(xref_pdb_split[[1]][1]) 


# check these pdb ID numbers against the alphafold website: https://alphafold.ebi.ac.uk/
# alternatively, fetch alphafold predictions from the UniProt accession numbers
af <- fetch_alphafold_prediction(Accessions)
pdb.1 <- as.character(ptsi_pgk_pdb_info[1,1])

# plugging this in to a browser will give the 3D visualization
# the pdb ID is the from your pdb results (from fetch_pdb)
# see this website: https://3dmol.csb.pitt.edu/doc/tutorial-url.html
# https://3Dmol.org/viewer.html?pdb=1EZA&select=chain:A&style=cartoon;stick:radius~0.1&surface=opacity:0.8;colorscheme:whiteCarbon&select=chain:B&style=cartoon;line&select=resi:19,23,26;chain:B&style=stick&labelres=backgroundOpacity:0.8;fontSize:14


# for some reason the r3dmol code below stopped working
# maybe one of you can get it to work!
# xyz <- af$P0A799[,7:9]
# write.pdb(pdb = as.matrix(xyz), file = "P0A799.pdb")

r3dmol(                         # Set up the initial viewer
  id = "1EZA",
  viewer_spec = m_viewer_spec(
    cartoonQuality = 10,
    lowerZoomLimit = 50,
    upperZoomLimit = 350
  )
) %>%
  m_add_model(                  # Add model to scene
    data = pdb.1,
    format = "pdb"
  ) %>%
  m_zoom_to() %>%               # Zoom to encompass the whole scene
  m_set_style(                  # Set style of structures
    style = m_style_cartoon(
      color = "#00cc96"
    )
  ) %>%
  m_set_style(                  # Set style of specific selection
    sel = m_sel(ss = "s"),      # (selecting by secondary)
    style = m_style_cartoon(
      color = "#636efa",
      arrows = TRUE
    )
  ) %>%
  m_set_style(                  # Style the alpha helix
    sel = m_sel(ss = "h"),      # (selecting by alpha helix)
    style = m_style_cartoon(
      color = "#ff7f0e"
    )
  ) %>%
  m_rotate(                     # Rotate the scene by given angle on given axis
    angle = 90,
    axis = "y"
  ) %>%
  m_spin()                      # Animate the scene by spinning it
