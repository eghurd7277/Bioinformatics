#filepath = "/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics"




mySeq <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Feb1_samp_seq/Feb1_Samp.fasta", format="fasta")

myFirstALign <- msa(mySeq)

myFirstALign

print(myFirstALign, show="complete")

LetterCont <- alphabetFrequency(myFirstALign)
ConsSeq <- consensusString(myFirstALign)
gaps <- str_count(ConsSeq,"-")


G <- sum(LetterCont[,3])
C <- sum(LetterCont[,2])
A <- sum(LetterCont[,1])
T <- sum(LetterCont[,4])
totalbp <- A+T+C+G
totalbp
GC_Cont <- (G+C)/totalbp


alignLength <- nchar(myFirstALign)
SecAlign <- msaConvert(myFirstALign,type="seqinr::alignment")

d <- dist.alignment(SecAlign, matrix = "identity")
dmat <- as.matrix(d)

Alignment_phyDat <- msaConvert(myFirstALign, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")



