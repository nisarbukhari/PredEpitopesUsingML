###############################################################
#                                                             #
#       Feature Extraction of SARS-CoV-2 peptide sequences    #
#           (T-Cell Epitopes and Non T-Cell Epitopes)         #
#                                                             #
###############################################################
###############################################################
#                                                             #
# Feature extracted based on the physicochemical properties   #
# of protein sequences                                        #                                                                           #
###############################################################

library(Peptides)
library(peptider)
library(seqinr)
library(Biostrings)
sequence<-read.table("peptides.csv")
i=1
for (i in 1:length(sequence)){
  
  cat("\n Completed: ",i,"/", length(sequence))
  
  peptide         <- paste(unlist(sequence[i]), collapse='')
  
  ############################################
  #     Part I: Single Properties
  ############################################
  
  # F1: aliphaticIndex
  F1_aliphaticIndex  <- aIndex(peptide)
  
  # F2: bomanIndex
  F2_bomanIndex      <- boman(peptide)

  # F3: instaIndex
  F3_instaIndex      <- instaIndex(peptide)
  
  # F4: probabilityDetectionPeptide
  F4_probabilityDetectionPeptide <- ppeptide(peptide,libscheme = "NNK", N=10^8)
  
  # F5: numberNeighbors
  #F5_numberNeighbors <- getNofNeighbors(peptide,blosum = 1,method = "peptide", libscheme = "NNK")

  # Merging of Part 1 results
  resultPart1=data.frame(F1_aliphaticIndex,F2_bomanIndex,F3_instaIndex,F4_probabilityDetectionPeptide)#,F5_numberNeighbors)
  
  ############################################
  #     Part II: Double Properties
  ############################################
  
  # F5: homentIndex
  F5_homentIndex1     <- hmoment(seq = peptide, angle = 100, window = 11)
  F5_homentIndex2     <- hmoment(seq = peptide, angle = 160, window = 11)

  # F7: molecularWeight
  F6_molecularWeight1 <-mw(seq = peptide,monoisotopic = FALSE)
  F6_molecularWeight2 <-mw(seq = peptide,monoisotopic = TRUE)

  # Merging of Part 2 results
  resultPart2=data.frame(F5_homentIndex1, F5_homentIndex2,F6_molecularWeight1,F6_molecularWeight2)
  
  
  ############################################
  #     Part III: Multiple Properties
  ############################################
  
  # F7: peptideCharge
  pKscale=c("Bjellqvist", "Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon", "Stryer")
  F7_peptideCharge=c()
  for (j in 1:length(pKscale)){
    x=charge(seq= peptide,pH= seq(from = 5,to = 9,by = 1), pKscale= pKscale[j])
    F7_peptideCharge = c(F7_peptideCharge,x)
  }
  names(F7_peptideCharge)<-paste("F7_pCharge",c(1:length(F7_peptideCharge)),sep='')
  

  # F8: Hydrophobibity for 44 scales
  scale=c("Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari", "Chothia", "Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman", "Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic", "Kidera", "Kuhn", "KyteDoolittle", "Levitt", "Manavalan", "Miyazawa", "Parker", "Ponnuswamy", "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", "Tanford", "Welling", "Wilson", "Wolfenden", "Zimmerman", "interfaceScale_pH8", "interfaceScale_pH2", "octanolScale_pH8", "octanolScale_pH2", "oiScale_pH8","oiScale_pH2")
  F9_hydro=c()
  for (j in 1:length(scale)){
    x=  hydrophobicity(seq = peptide,scale = scale[j])
    F8_hydro=c(F8_hydro,x)
  }
  names(F8_hydro)<-paste("F8_hydro",c(1:length(F8_hydro)),sep='')

  
  # F9: isoElectricPoint at 9 pKscale
  pKscale=c("Bjellqvist", "EMBOSS", "Murray", "Sillero", "Solomon", "Stryer", "Lehninger", "Dawson","Rodwell")
  F10_isoElectricPoint=c()
  for (j in 1:length(pKscale)){
    x=  pI(peptide, pKscale = pKscale[j])
    F9_isoElectricPoint=c(F9_isoElectricPoint,x)
  }
  names(F9_isoElectricPoint)<-paste("F9_isoEP",c(1:length(F9_isoElectricPoint)),sep='')
  

  # F10: kideraFactors
  F10_kFactors       <- as.numeric(unlist(kideraFactors(seq = peptide)))
  names(F10_kFactors)<-paste("F10_kFactors",c(1:length(F10_kFactors)),sep='')


  # F11: aaComp
  F11_aaComp          <- as.numeric(unlist(aaComp(peptide)))
  names(F11_aaComp)   <-paste("F11_aaComp",c(1:length(F11_aaComp)),sep='')


  # Merging of Part 2 results
  resultPart3=c(F7_peptideCharge,F8_hydro,F9_isoElectricPoint,F10_kFactors,
                F11_aaComp)

  finalResult = c(resultPart1,resultPart2,resultPart3)
  
  if(i==1){
    write.table(finalResult, "dataSet.csv", sep = ",", row.names=F,col.names = T)
  }
  else{
    write.table(finalResult, "dataSet.csv", sep = ",", row.names=F, col.names = F, append = T)
  }

}

cat("\n Done")