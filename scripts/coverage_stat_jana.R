library(data.table)

# get arguments
args = commandArgs(trailingOnly=TRUE)

# load file
file <- read.table(args[1], sep="\t", skip = 1)
colnames(file) <- c("Chrom", "Start", "End", "Exon", "Base_number", "Coverage")

file_name <- basename(args[1])
path <- dirname(args[1])
file$sample <- gsub(".PBcov.cons.txt", "", file_name)

# create table with coverage per exon statistics
file_ID <- file$sample[1]
Exon_stat <- list()
#Exon_stat <- file[, list(Coverage_max = max(Coverage), Coverage_min = min(Coverage), Coverage_mean = round(mean(Coverage), digit = 2),
#                         Coverage_median = median(as.double(Coverage))), by=.(Chrom, Exon, Start, End)]

  for(j in levels(file$Exon)){
    
    # only particular exon
    tmp <- file[grep(paste0("^",j,"$"),file$Exon),]
    
    # count statistics
    tmp_final <- data.frame(Chrom = tmp$Chrom[1], Exon = tmp$Exon[1], Start = tmp$Start[1], End = tmp$End[1],
                            Coverage_max = max(tmp$Coverage), Coverage_min = min(tmp$Coverage),
                            Coverage_mean = round(mean(tmp$Coverage), digit = 2), Coverage_median = median(tmp$Coverage))
    Exon_stat[[j]] <- tmp_final
  }
  
#save file
Exons_sample_stat <- rbindlist(Exon_stat)
#Exons_sample_stat <- Exon_stat
write.table(Exons_sample_stat, file = paste0(path,"/",file_ID,".perexon_stat.txt"), quote = F, row.names = F, sep = "\t")




