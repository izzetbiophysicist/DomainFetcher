library(bio3d)

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
pdb_vec <- read.csv(args[1], header=FALSE)
pdb_vec <- pdb_vec[,1]

pdb_pfam_map <- read.table("pdb_pfam_mapping.txt", skip=1, sep="\t", header = 1)

pdb_pfam_map$PDB <- toupper(pdb_pfam_map$PDB)

pdb_not_found <- c()
removed <- c()
system("mkdir domains")

  
for(current_pdb in pdb_vec){
    
    subset_pfam <- pdb_pfam_map[which(pdb_pfam_map$PDB == current_pdb),]
    
    removed <- c(removed,subset_pfam$PDB[which(subset_pfam$AUTH_PDBRES_START == "None")])
    removed <- c(removed,subset_pfam$PDB[which(subset_pfam$AUTH_PDBRES_END == "None")])
    
    w_nostart <- which(subset_pfam$AUTH_PDBRES_START == "None")
    if(length(w_nostart) != 0){
      subset_pfam <- subset_pfam[-which(subset_pfam$AUTH_PDBRES_START == "None"),]
    }
    w_noend <- which(subset_pfam$AUTH_PDBRES_END == "None")
    if(length(w_noend) != 0){
      subset_pfam <- subset_pfam[-which(subset_pfam$AUTH_PDBRES_END == "None"),]
    }
   
    
    if(nrow(subset_pfam) != 0){
      pdb <- tryCatch(read.pdb(current_pdb), error=function(e) NULL)
      if(is.null(pdb) != TRUE){
      
       for(j in 1:nrow(subset_pfam)){
         
         res_start <- subset_pfam$AUTH_PDBRES_START[j]
         res_stop <- subset_pfam$AUTH_PDBRES_END[j]
         
         domain <- trim.pdb(pdb, atom.select(pdb, chain=subset_pfam$CHAIN[j], resno=c(res_start:res_stop)))
         system(paste("mkdir ./domains/", subset_pfam$PFAM_ACCESSION[j],sep=""))
         
         write.pdb(domain, file=paste("./domains/",subset_pfam$PFAM_ACCESSION[j],"/",current_pdb,"_",subset_pfam$CHAIN[j], ".pdb", sep=""))
         
        }
      } else {
        pdb_not_found <- c(pdb_not_found, current_pdb)
      }
    }
  }

write.csv(pdb_not_found, "not_found.csv")
write.csv(removed, "removed.csv")