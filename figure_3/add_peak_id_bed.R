#clear workspace
rm(list = ls())
graphics.off()

#script to call in bash files to add peak ids to bed files, get into correct format before running Homer

args = commandArgs(trailingOnly=T)

for(i in 1:length(args)){ 

  #read in peak file, originally narrowPeak
  trim <- read.table(args[i], sep="\t", header=F) #a tsv file with header, in bed format otherwise
  
  #create peak ids
  paste_fxn <- function(i){ paste0("peak_",i) }

  #add them to peak file
  peak_ids <- sapply( 1:nrow(trim), paste_fxn)

  #select columns
  trim_2 <- cbind(trim[,c(1:3)], peak_ids) #keep only chr, start, end columns
  
  #process file name
  file_split <- strsplit(args[i],"/")
  
  count <- length(file_split[[1]])
  
  file_name <- file_split[[1]][count]
  
  file_path <- strsplit(args[i], file_name)[[1]]
  
  full_path <- paste(file_path, "homer_",file_name,sep = "")
  
  #save new bed file with peak ids, adding prefix homer_
  write.table(trim_2, file=full_path, sep="\t", row.names = F, col.names = F, quote = F) 
  
}