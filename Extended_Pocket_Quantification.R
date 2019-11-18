#Import Required Libraries
Sys.setenv(R_INSTALL_STAGED=FALSE) #For Ubuntu on windows setup

if(!require(bios2mds)){
  install.packages("bios2mds", repos = "http://cran.us.r-project.org")
  library(bios2mds)
}

if(!require(devtools)){
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  library(devtools)
}

if(!require(phylotools)){
  install.packages("phylotools", repos = "http://cran.us.r-project.org")
  library(phylotools)
}

if(!require(seqinr)){
  install.packages("seqinr", repos = "http://cran.us.r-project.org")
  library(seqinr)
}

if(!require(dplyr)){
  install.packages("dplyr", repos = "http://cran.us.r-project.org")
  library(dplyr)
}

if(!require(openxlsx)){
  install.packages("openxlsx", repos = "http://cran.us.r-project.org")
  library(openxlsx)
}

install_github("mhahsler/rMSA") #MSA installation, requires devtools
library(rMSA)


if(!require(data.table)){
  install.packages("data.table", repos = "http://cran.us.r-project.org")
  library(data.table)
}

if(!require(ggplot2)){
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  library(ggplot2)
}

#install.packages("str_sub")
if(!require(stringr)){
  install.packages("stringr", repos = "http://cran.us.r-project.org")
  library(stringr)
}

#############################################################

#Read User Inputs into Variables
args <- commandArgs(trailingOnly = TRUE)
qname <- args[1]
PocketVar <- args[2]

if (PocketVar == "B") {
  Pocket <- c(7,9,24,34,45,63,66,67,70,99) #B-Pocket
}
if (PocketVar == "F") {
  Pocket <- c(77,80,81,84,94,115,122,142,145,146) #F-Pocket
}

#Temporary User Variables
qseq <- unlist(read.fasta(Sys.glob(file.path(paste("./InputFiles/",qname,sep = ""),"*.fasta")),seqonly = TRUE))

refseq <- "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV" #HLA-A2:01:01 sequence

#Align Query Sequence to alignment files

unlink("./Alignments_Copy", recursive = TRUE)
dir.create("Alignments_Copy")
file.copy("Alignments/.", "./Alignments_Copy/.", recursive = TRUE, overwrite = TRUE)

Alignnames <- list.files("./Alignments")
Alignnames <- substr(Alignnames, 1, nchar(Alignnames)-6)

for (i in 1:length(Alignnames)) {
clean.fasta.name(infile = paste0("Alignments_Copy/",Alignnames[i],".fasta") , outfile = paste0("Alignments_Copy/",Alignnames[i],"_clean.fasta")) #Removes spaces
}

for (i in 1:length(Alignnames)) {
  write.fasta(refseq, "HLA-A2:01:01:01", paste0("./Alignments_Copy/",Alignnames[i],"_clean.fasta"), open = "a", nbchar = 1000, as.string = TRUE) #Adds HLA Ref
}

for (i in 1:length(Alignnames)) {
  write.fasta(qseq, qname, paste0("./Alignments_Copy/",Alignnames[i],"_clean.fasta"), open = "a", nbchar = 1000, as.string = TRUE) #Adds query sequence
}

Addon <- data.frame(c("HLA-A2:01:01:01", qname), c(refseq, qseq)) #Sequence addon for later (adding names)
names(Addon) <- c("seq.name", "seq.text")

#Multiple Sequence Alignment

Alignment <- as.list(NULL)
for (i in 1:length(Alignnames)) {
  Alignment[[i]] <- readAAStringSet(paste0("./Alignments_Copy/",Alignnames[i],"_clean.fasta"))
}

Alignment_2 <- NULL
temp2 <- NULL
for (i in 1:length(Alignnames)) {
  tempname <- c(names(readAAStringSet(paste0("Alignments_Copy/",Alignnames[i],".fasta"))), "HLA-A2:01:01:01", qname) #Imports names
  temp <- mafft(Alignment[[i]], param = "--auto --clustalout --inputorder --amino")
    temp2 <- data.frame(temp@unmasked, tempname)
      Alignment_2[[i]] <- setDT(temp2, keep.rownames = TRUE)
}
temp <- NULL

#possible additional argument --lop -3 --lep -1

Alignment_2[[1]][!grepl("x", Alignment_2[[1]]$temp.unmasked), ]

saveRDS(Alignment_2, file = "Alignment_Complete.RDS")
Alignment_2 <- readRDS("Alignment_Complete.RDS")

#Alignment Position Alignment

aa <- c("G","A","R","H","K","D","E","S","T","N","Q","C","P","V","I","L","M","F","Y","W")

#Generates position for all amino acids in ref and query sequences for all alignments
poslist_qer_all <- NULL
poslist_ref_all <- NULL
poslist_qer <- NULL
poslist_ref <- NULL
for (j in 1:length(Alignment_2)) {
  for (i in 1:length(aa)) {
    poslist_qer[i] <- gregexpr(aa[i], Alignment_2[[j]][,2][(nrow(Alignment_2[[j]]) - 0)]) 
    poslist_ref[i] <- gregexpr(aa[i], Alignment_2[[j]][,2][(nrow(Alignment_2[[j]]) - 1)]) 
  }
  names(poslist_qer) <- aa
  poslist_qer_all[[j]] <- poslist_qer
  names(poslist_ref) <- aa
  poslist_ref_all[[j]] <- poslist_ref
  poslist_qer <- NULL
  poslist_ref <- NULL
}
names(poslist_qer_all) <- Alignnames
names(poslist_ref_all) <- Alignnames

#Pocket Labelling
#Pocket <- c(7,9,24,34,45,63,66,67,70,99) #B_Pocket - HLA
Pocket <- Pocket + 24 #Adjustment to include signal sequence for reference HLA
Pocket <- c(Pocket, (Pocket+1), (Pocket-1))

Alignment_refsequence_adjustment <- NULL
Alignment_refsequence_adjustment_OR <- NULL
Alignment_qersequence_adjustment <- NULL
Alignment_qersequence_adjustment_OR <- NULL
for (i in 1:length(poslist_ref_all)) {
Alignment_refsequence_adjustment[[i]] <- do.call(rbind, Map(data.frame, A = poslist_ref_all[[i]], B = aa))
Alignment_refsequence_adjustment_OR[[i]] <- Alignment_refsequence_adjustment[[i]][order(Alignment_refsequence_adjustment[[i]]$A),]
Alignment_qersequence_adjustment[[i]] <- do.call(rbind, Map(data.frame, A = poslist_qer_all[[i]], B = aa))
Alignment_qersequence_adjustment_OR[[i]] <- Alignment_qersequence_adjustment[[i]][order(Alignment_qersequence_adjustment[[i]]$A),]
}
Alignment_refsequence_adjustment <- NULL
Alignment_qersequence_adjustment <- NULL
names(Alignment_refsequence_adjustment_OR) <- Alignnames
names(Alignment_qersequence_adjustment_OR) <- Alignnames

#Calculate Binding Probabilities
#install.packages("str_sub")

BP <- read.table(paste0("InputFiles/",qname,"/prob.txt"), sep = "\t")
Probabilities <- as.numeric(str_sub(BP[,1], -5, -1))
Position <- c(1:nrow(BP))
AminoAcids <- str_sub(BP[,1], 3, 5)
BP_Results <- data.frame(Position, AminoAcids, Probabilities)

#Scoring Alignment
PMBEC <- read.xlsx("PMBEC_lower_-0.2gap.xlsx", colNames = TRUE, rowNames = TRUE) #PMBEC Matrix

length2NA <- function(x) { #Removes 0 length entries due to containing unspecified characters e.g. "x"
  if (length(x)==0) {
    return(NA)
  }
  else {
    return(x)
  }
}

length2Zero <- function(x) { #Converts 0 length entries to 0 probability due to containing unspecified sequence positions in BP_Results
  if (length(x)==0) {
    return(0)
  }
  else {
    return(x)
  }
}

Refscores <- as.list(NULL)
Final_Data <- as.list(NULL)
for (k in 1:length(Alignment_2)) {
  tempdf <- data.frame(Alignment_qersequence_adjustment_OR[[k]], BP_Results)
  Inputlist <- as.vector(0)
  for (j in 1:nrow(Alignment_2[[k]])) {  #Loop along MHCs
    score <- 0
    for (i in 1:length(Alignment_qersequence_adjustment_OR[[k]][,1])) { #Loop along sequence
      if(length2Zero(BP_Results[which(tempdf$A == Alignment_refsequence_adjustment_OR[[k]][i,]$A),3]) < 0.05) {
        score <- score + 0
      }
      else if(BP_Results[i,1] %in% Pocket) {
        score <- score + ((PMBEC[substring(Alignment_2[[k]][nrow(Alignment_2[[k]]),2],
                                          Alignment_refsequence_adjustment_OR[[k]][i,1],
                                          Alignment_refsequence_adjustment_OR[[k]][i,1]),
                                substring(Alignment_2[[k]][j,2],Alignment_refsequence_adjustment_OR[[k]][i,1],
                                          Alignment_refsequence_adjustment_OR[[k]][i,1])])*
          BP_Results[which(tempdf$A == Alignment_refsequence_adjustment_OR[[k]][i,]$A),3])
      }
      else {
        score <- score + 0
      }
      Inputlist[j] <- length2NA(score)
    }
  } 
  Final_Data[[k]] <- Inputlist
}



temp2 <- NULL
for(k in 1:length(Alignment_2)) {
  temp <- data.frame(Alignment_2[[k]]$tempname, Final_Data[[k]], Alignment_2[[k]]$temp.unmasked, 1:nrow(Alignment_2[[k]]), Alignnames[[k]])
  temp2 <- rbind(temp, temp2)
}

temp3 <- temp2[-(which(temp2$Alignment_2..k...tempname %in% qname)),] #Removes rows with query name match

temp4 <- data.frame(temp3, recode(temp3[,5], "A_prot_2" = "HLA_A",
                                  "A_prot" = "HLA_A",
                                  "B_prot" = "HLA_B",
                                  "B_prot_2" = "HLA_B",
                                  "C_prot" = "HLA_C",
                                  "C_prot_2" = "HLA_C",
                                  "Sus_2" = "Sus",
                                  "Sus_3" = "Sus",
                                  "Sus_1" = "Sus",
                                  "Patr_C" = "Patr",
                                  "Patr_B" = "Patr",
                                  "Patr_A" = "Patr"))

                  
        

###
#GRAPHS#
boxplot_1 <- ggplot(temp4, aes(x=temp4[,6], y=temp4[,2])) + geom_boxplot() +
  xlab("MHC Supertype") +
  ylab("Pocket Alignment Score") +
  ggtitle(paste(qname," ",PocketVar,"_Pocket", sep = "")) +
  geom_hline(yintercept = max(temp2[,2], na.rm = TRUE)) +
  geom_hline(yintercept = (max(temp2[,2], na.rm = TRUE)*0.9), linetype="dashed")

png(filename = paste("OutputFiles/",qname," ",PocketVar,"_Pocket.png", sep = ""), width = 735, height = 556)
boxplot_1
dev.off()

save.image(file = paste("OutputFiles/",qname," ",PocketVar,"_Pocket.RData", sep = ""))
