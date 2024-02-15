library(tidyr)
library(dplyr)
library(UniprotR)
library(protti)
library(stringr)

#amino acid sequence from last project written to new fasta file
writeXStringSet(rnaSeq1, "rufusAA.fasta", format="fasta")

#protein ascession numbers from uniprot read in
ProtSeq <- read.csv("protein_seq.fasta")
#needed as char string
AscNum <- ProtSeq[,1]

#gene ontology information
GO_Info <- GetProteinGOInfo(AscNum)
PlotGoInfo(GO_Info)

#from github, makes the GO info image and saves to my github folder
PlotGOAll(GOObj = GO_Info, Top = 10, directorypath = getwd(), width = 8, height = 5)

#gets the pathology info from uniprot 
PathInfo <- GetPathology_Biotech(AscNum)
Get.diseases(PathInfo, directorypath = getwd())

#gets info from uniprot
UnProtData <- fetch_uniprot(AscNum)

#no pbdID numbers returned, used ones given in hw
samppdbIDs <- c("1ZMR","2HWG")
str_info <- fetch_pdb(samppdbIDs)
threeDstruc <- fetch_alphafold_prediction(uniprot_ids = samppdbIDs, return_data_frame = TRUE)
# went to office hours on 2/15, process is right but when running code no disease information is returned, no protein database IDs either
# tried to run the code with the ones given in the homework, returns error 404 not found
#could not access 3D structure, so couldn't save that to the github folder