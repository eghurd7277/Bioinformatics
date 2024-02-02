#filepath = "/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics"




mySeq <- readDNAStringSet("/Users/emmagrace3/Desktop/bioinformatics/GitHub/Bioinformatics/Feb1_samp_seq/Feb1_Samp.fasta", format="fasta")

myFirstALign <- msa(mySeq)

myFirstALign

print(myFirstALign, show="complete")

