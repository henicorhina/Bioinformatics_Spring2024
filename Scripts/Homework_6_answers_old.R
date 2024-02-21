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

# get some accession numbers from an online repository
Accessions <- GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv") 
TaxaObj <- GetNamesTaxa(Accessions) 
PlotGenesNetwork(TaxaObj)
GeneOntologyObj <- GetProteinGOInfo(Accessions) 
PlotGOBiological(GeneOntologyObj, Top = 10) 

alb_aa_blast <- read_csv("Diomedea_exulans/Book1.csv")
GetNamesTaxa(alb_aa_blast$Accession)
GeneOntologyObj <- GetProteinGOInfo(alb_aa_blast$Accession) 

my_pdb_ids <- c("P68871")

ptsi_pgk_pdb_information <- fetch_pdb(my_pdb_ids)

utils::data("ptsi_pgk")
uniprot_ids <- unique(ptsi_pgk$pg_protein_accessions)
uniprot_information <- fetch_uniprot(uniprot_ids)
xref_pdb <- uniprot_information$xref_pdb
af <- fetch_alphafold_prediction(uniprot_ids)

IDs <- str_split(xref_pdb, pattern = ";")
IDs <- c(IDs[[1]], IDs[[2]])
IDs <- IDs[IDs != ""]

pdb_information <- fetch_pdb(IDs)

AF_preds <- fetch_alphafold_prediction(uniprot_ids) 


r3dmol(                         # Set up the initial viewer
  viewer_spec = m_viewer_spec(
    cartoonQuality = 10,
    lowerZoomLimit = 50,
    upperZoomLimit = 350
  )
) %>%
  m_add_model(                  # Add model to scene
    data = pdb_6zsl,
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
