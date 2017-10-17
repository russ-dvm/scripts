library(tidyr)
library(plyr)
library(dplyr)

##Read in the various disease/pathogen subgroups
missing.file.list <- list.files(path="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/", pattern="missing.annotated")
path="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/"
all_missing <- paste(path, missing.file.list, sep="")

##The data is read in as a list of dataframes, one for each condition. Labels -must- be applied using the python script "assign_info_to_plink_output.py" before hand. The filenames are not used here to impart information to the data frames.
all_missing_values <- lapply(all_missing, read.table, h=T)

#Do a bit of sanity checking:
#Does the number of diseased animals (phenotyped and non-phenotyped) always add up? (total should be 421)
#The missed calls should be the same across each disease subtype for each rsID in controls
check_diseased <- ldply(all_missing_values, function(x){
    non_pheno <- x[1,]$N_GENO
    dis <- x[3,]$N_GENO
    total <- non_pheno+dis
    condish <- x[2,]$condition
    total_condish <- data.frame(condish, total)
    return(total_condish)
    
})
# View(check_diseased)
check_normal <- ldply(all_missing_values, function(x){
  
  normals <- filter(x, CLST == 1)
  big_frame <- data.frame(normals$SNP, normals$N_MISS, normals$condition)
  return(big_frame)
  
})
check_normal <- spread(check_normal, normals.condition, normals.N_MISS)
# View(check_normal)
normal_values <- check_normal[c(1,2)]
colnames(normal_values) <- c("SNP", "Healthy")

##Make a df of all the missing values for each group (column) for each rsid (row)
diseased_only <- ldply(all_missing_values, filter, CLST == 2)

#Remove columns that still contain variable values, or the spread step won't accomplish what we want
diseased_only <- diseased_only[c(-5,-6,-7)]
diseased_only <- spread(diseased_only, condition, N_MISS)

#Add in the healthy values
complete_missing_df <- merge(normal_values, diseased_only)

write.table(complete_missing_df, file="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/missing_values.txt", quote=F, sep="\t", row.names = F)


