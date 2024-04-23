#loading necessary packages
library(seqinr)
library (Biostrings)
library(BiocManager)
library(msa)
library(phangorn)
library(ape)

#reading in fasta files
catus <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/catus.fasta", format = "fasta")
silvestris <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/silvestris.fasta", format = "fasta")
nigripes <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/nigripes.fasta", format = "fasta")
margarita <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/margarita.fasta", format = "fasta")
chaus <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/chaus.fasta", format = "fasta")

#combining fasta files into one, rewriting seqence names
sequences <- c(catus,silvestris,nigripes,margarita,chaus)
names(sequences) <- c("catus_GU300839.1", "catus_GU300838.1", "silvestris_GU300940.1", "silvestris_GU300937.1", "nigripes_GU300900.1", "nigripes_GU300897.1", "margarita_GU300871.1", "margarita_GU300869.1", "chaus_GU300852.1", "chaus_GU300851.1")


#creating a MSA using muscle
alignment <- msa(sequences, method="Muscle")

#converting to seqinr format
aln2 <- msaConvert(alignment, type="seqinr::alignment")

#computing distance matrix
dmat <- dist.alignment(aln2,matrix="identity")
distmat <- as.matrix(dmat)
mydistmat <- as.dist(distmat)

#building a neighbor joining tree
myNJtree <- ape::nj(mydistmat)

#plotting the nj tree unrooted
plot(myNJtree, "unrooted")
#plotting the nj tree rooted
plot(myNJtree)

#Writing alignment to new file to do RaxML tree
write.phylip(alignment, filepath = "/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Final_Presentation/Alignment")