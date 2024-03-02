#MIDTERM 1 
#wd = /Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics

#loading necessary libraries
library(msa)
library(tidyr)
library(UniprotR)
library(protti)
library(stringr)
library(dplyr)
library(seqinr)
library(stringi)

#reading in the fasta file and assigning to variable using msa functions
MidSequences <- readDNAStringSet("sequences.fasta", format="fasta")

#creating an alignment and assigning to variable
MidAlign <- msa(MidSequences)
print(MidAlign, show = "complete")
#making alignment a matrix so can view only specific sequences together at a time
AlMat <- as.matrix(MidAlign)

#converting to seqinr type, setting up distance matrix to look at which sequences differ from others
MidAlign2 <- msaConvert(MidAlign, type="seqinr::alignment")
dis <- dist.alignment(MidAlign2, "identity")
distmatrix <- as.matrix(dis)
#Homo_sapiens_4, Homo_sapiens_6, and Homo_sapiens_10 differ from the others
#looking at alignment to see what is observable

#made smaller matrix of each differing sequence plus one as a reference, alignments not subsettable
compAlign <- AlMat[c(1:3,20),]

#Homo_sapiens_6 has a gap at position one (likely a deletion mutation), a mismatch at position 3 (A for C), 
#a mismatch at position 47 (G for A), a mismatch at position 134 (C for T),a mismatch at position 145 (T for A), 
#a mismatch at position 152 (G for C), and mismatch at position 586 (C for G)

#Homo_sapiens_4 has a mismatch at position 39 (A for C)

#Homo_sapiens_10 has a mismatch at position 39 (G for C), a mismatch at position 45 (T for A)

Consensus <- consensusString(MidAlign)
print(Consensus)

#consensus sequence copy-pasted to BLAST, result is the hemoglobin subunit beta (HBB) gene, accession num LC121775.1

#homosapiens 6 is most dissimilar (from looking at dist. matrix), made ch string for just that sequence and translated
hs6_dna <- as.character(MidSequences[6,])
hs6protseq <- translate(s2c(hs6_dna))

#wrote to fasta file, saves to github folder
write.fasta(hs6protseq,"Homo_sapiens_6_AA","hs6AA.fasta")
#Accession number for best match entry through UniProt is A0A0J9YWK4

#running the code similar to that of HW6 returned my same issue of getting a NULL result
# I chose to then obtain the information from the database myself, looking at closest match 
# this gene is associated with beta-thalessemia and heinz body anemia
# it does not look like the person has the disease, looking at which amino acid changes would result in pathogenic effects, none of them appear to line up with the sample sequence

#3d structure screenshot from UniProt saved in GitHub folder as Midterm1_ProtStruc
