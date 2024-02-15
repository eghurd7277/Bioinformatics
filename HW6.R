writeXStringSet(rnaSeq1, "rufusAA.fasta", format="fasta")
ProtSeq <- read.csv("protein_seq.fasta")
AscNum <- ProtSeq[,1]
GO_Info <- GetProteinGOInfo(AscNum)
PlotGoInfo(GO_Info)
PlotGOAll(GOObj = GO_Info, Top = 10, directorypath = getwd(), width = 8, height = 5)
PathInfo <- GetPathology_Biotech(AscNum)
Get.diseases(PathInfo, directorypath = getwd())
UnProtData <- fetch_uniprot(AscNum)
uniprotIDs <- 
samppdbIDs <- c("1ZMR","2HWG")
str_info <- fetch_pdb(samppdbIDs)
threeDstruc <- fetch_alphafold_prediction(AscNum)