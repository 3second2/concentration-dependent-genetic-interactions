# generating pymol script with colour
# code adapted based on P.Julien's script
#### 
#### Remember to change the input file and output file locations 
# line 33, OUTPUT_FILE path may need to be changes

load("path/positional_z_score.RData") # positional median, min, max z-score matrix to examine the sensitivity to mutational effects
pymol_scores<- list(l=c(), h=c()) # low, high
pymol_scores[[1]]<- positional_matrice[3, ] # low, positional median
pymol_scores[[2]]<- positional_matrice[4, ] # high, positional median
library(seqinr)

file_names= list(
  "ci_L_z_score_strucutral_median" , 
  "ci_H_z_score_strucutral_median"
)

# WT <- "CTTAAAGCAATTTATGAAAAAAAGAAAAATGAACTTGGCTTATCCCAGGAATCTGTCGCAGACAAGATGGGGATGGGGCAGTCAGGCGTTGGTGCTTTATTTAATGGCATCAATGCATTAAATGCTTATAACGCCGCATTGCTTGCAAAAATTCTCAAAGTTAGCGTTGAAGAATTT"
PDB_ID <- "3bdn"

colscale <- colorpanel(n=100, low="blue", hi="violet", mid="gray80") # gray80, yellow, and blues as in the scatter plot
defaut_color <- "[0.25,0.25,0.25]" # Color for protein residues not part of the domain of interest (as fraction of R G and B)
DNA_col <- "[1,1,1]" #   Color of DNA molecule (as fraction of R G and B) , 
WT_AA <- "LKAIYEKKKNELGLSQESVADKMGMGQSGVGALFNGINALNAYNAALLAKILKVSVEEF"
WT_AA_vec <- s2c(WT_AA)

possible_AAs<- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AA1to3 <- c(A="ALA", C="CYS", D="ASP", E="GLU", F="PHE", G="GLY", H="HIS", I="ILE", K="LYS", L="LEU", M="MET", N="ASN", P="PRO", Q="GLN", R="ARG", S="SER", T="THR", V="VAL", W="TRP", Y="TYR")

range_structure <- 18:(length(WT_AA_vec)+17) # Overlap between structure file and
for (i in 1: 2) {
  names(pymol_scores[[i]])<- c(1:59)
  OUTPUT_FILE= paste ("pymol_script_ci_",file_names[i], sep="")
  
  scores_offset <- as.numeric(names(pymol_scores[[i]])) + 17
  names(scores_offset) <-names(pymol_scores[[i]])
  
  pymol_scores2 <- pymol_scores[[i]][scores_offset[names(pymol_scores[[i]])] %in% range_structure] 
  names(pymol_scores2) <- as.character(range_structure)
  
  ### This is just in case part of the AA sequence is not in the pdb (otherwise this just returns the normal variant sequence)
  AAs_Struc <- WT_AA_vec[scores_offset[names(pymol_scores[[i]])] %in% range_structure]
  
  #### Assigning colors to AA based on their median score
  colStructure <- colscale[cut(c(pymol_scores2, max(abs(range(pymol_scores2))) * c(-1,1) ), 100)] # The max stuff is to have the scale centered on 0.
  
  
  #### This matrix associates all indices together
  IndexCorr <- matrix(c(1:length(WT_AA_vec), scores_offset, scores_offset), ncol=3, byrow=F, dimnames=list(WT_AA_vec, c("VariantIndex", "PDBIndex", "ProteinIndex")))
  
  # Redirection to script file
  
  sink(OUTPUT_FILE)
  # General commands
  #cat("load /Volumes/Nestor/Projects/Epistasis/Data/General/1K9Q.pdb\n")
  cat("fetch ", PDB_ID, "\n", sep="")
  cat("show_as cartoon\n")
  
  # Colouring the backbone
  cat("set_color default_col, ", defaut_color ,"\n")
  cat("color default_col\n")
  
  # Colouring the DNA
  cat("set_color DNA_col, ", DNA_col, "\n")
  cat("select DNA1, chain C\n")
  cat("color DNA_col, DNA1\n")
  cat("select DNA2, chain D\n")
  cat("color DNA_col, DNA2\n")
  
  # A and B chains as sphere. Comment lines below to disable
  cat("select chainA, chain A\n")
  cat("show_as spheres, chainA\n")
  
  cat("select chainB, chain B\n")
  cat("show_as spheres, chainB\n")
  
  # Now colouring each individual residue
  for (j in 1:length(range_structure)) {
    
    n <- paste("sel",j,sep="")  
    cat("select ", n, ", /", PDB_ID, "//A/", AA1to3[AAs_Struc[j]], "`", names(pymol_scores2)[j], "\n", sep="")
    cat("set_color newcol" ,j, ", [", sep="") ;cat(round(col2rgb(colStructure[j])[,1] / 255, 3), sep=","); cat("]\n", sep="")
    cat("color newcol", j,", ", n ,"\n", sep="")
    #	cat("delete sele\n", sep="")
    
    # Duplicating to also colour 
    n <- paste("sel",j,sep="")	
    cat("select ", n, ", /", PDB_ID, "//B/", AA1to3[AAs_Struc[j]], "`", names(pymol_scores2)[j], "\n", sep="")
    cat("set_color newcol" ,j, ", [", sep="") ;cat(round(col2rgb(colStructure[j])[,1] / 255, 3), sep=","); cat("]\n", sep="")
    cat("color newcol", j,", ", n ,"\n", sep="")
    #	cat("delete sele\n", sep="")
    
  }
  sink()
  
}
