##Aim: Generate input for picrust2
###Otu table in biomformat 

require("biomformat")

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(readr)
library(readxl)
library(Biostrings)


seqtab.nochim.ratshiperammon<- readRDS("/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/outputs/seqtab.nochim.ratshiperammon.Rds")

##Biom file asv table
asvmat<- as.matrix(t(seqtab.nochim.ratshiperammon))##Rows should be the ASVs and columns the samples

tmp<- as.data.frame(rownames(asvmat))
tmp[,2]<- paste0("ASV", 1:nrow(asvmat))
colnames(tmp)<- c("Sequence", "ASV")
rownames(asvmat) <- paste0("ASV", 1:nrow(asvmat))
biom.tmp<- make_biom(asvmat, matrix_element_type = "int")
write_biom(biom.tmp, biom_file = "/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/input_ratshiperammon/asv_ratshiperammon.biom")


##Sequences with corrected names
library(DECIPHER)

dna <- readDNAStringSet("/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/outputs/ASV.fasta")
names(dna)<- paste0("ASV", 1:length(dna)) 
writeXStringSet(dna, filepath = "/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/input_ratshiperammon/dna_ratshiperammon.fasta") 