library(plyr)
library(tidyr)
library(dplyr)

###Function that takes a dataframe and false discovery rate as inputs, and calculates the Benjamini-Hochberg critical value and the significance of all rows in the data frame. Requires the original P value to be in a column named "P"
benj_hoch <- function(x, fdr){
  
  #order data by p value, from lowest to highest
  sorted.by.p <- x[order(x$P, decreasing=F),]
  #determine number of rows
  rownum <- nrow(sorted.by.p)
  #assign ranks to each row, lowest p value = 1
  sorted.by.p$rank <- c(1:rownum)
  #calculate the benjamini hochberg critical value - the integer = False Discovery Rate
  benj.hoch <- (sorted.by.p$rank/rownum)*fdr
  #add benjamini hochbergs to the data frame
  sorted.by.p$bh <- benj.hoch
  #determine if the p value less than the BH critical value
  sig <- sorted.by.p$P < sorted.by.p$bh
  #add the significance to the data frame
  sorted.by.p$significant <- sig
  #Take a subset composed only of values that are significant
  all.significant <- subset(sorted.by.p, significant == TRUE)
  #Take the last significant one (i.e. lowest rank)
  lowest.significant <- tail(all.significant, 1)
  lowest.rank <- lowest.significant$rank
  #"The largest P value that is greater than the BH critical value is significant, and ALL of the P values smaller than it are also significant, even the ones that aren't less than their BH critical value" http://www.biostathandbook.com/multiplecomparisons.html
  if (length(lowest.rank) > 0) {
    for (y in 1:lowest.rank){
      sorted.by.p$significant[y] <- TRUE
    }
  }  
  
  return(sorted.by.p)
  
}

##Global disease vs normal
global <- read.table("~/porcine/sequenom_results/plink_analysis/pig_sequenom.assoc.fisher.annotated", h=T, stringsAsFactors = F)

##Read in the various disease/pathogen subgroups
file.list <- list.files(path="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/", pattern="fisher.annotated")
path="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/"
all_files <- paste(path, file.list, sep="")

##The data is read in as a list of dataframes, one for each condition. Labels -must- be applied using the python script "assign_genes_to_plink_output.py" before hand. The filenames are not used here to impart information to the data frames.
all_conditions <- lapply(all_files, read.table, h=T)

##Perform multiple testing correction on all the dataframes using ldply (from the plyr package), which applies the function over each list item and returns a single dataframe composed of all list items.
all_conditions_df <- ldply(all_conditions, benj_hoch, 0.05)

#use tidyr to combine the frequency of alt_allele, whether it is BH-significant, and the original P value into a singal column
all_conditions_df_united <- unite(all_conditions_df, sig_p, c(significant, P))
all_conditions_df_united <- unite(all_conditions_df_united, freq_sig_p, c(F_A, sig_p), sep="<")

#Remove extraneous variables that eff up the spread step
all_conditions_df_united_without_var <- all_conditions_df_united[c(-8, -11, -12)]

#Rename columns
all_conditions_df_united_without_var <- unite(all_conditions_df_united_without_var, Healthy, F_U)
all_conditions_df_united_without_var <- unite(all_conditions_df_united_without_var, Ref, A2)
all_conditions_df_united_without_var <- unite(all_conditions_df_united_without_var, Alt, A1)

#Spread (un-melt?) the data, so that you have a table similar to the Keirstead paper
final <- spread(all_conditions_df_united_without_var, condition, freq_sig_p)

write.table(final, file="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/all_data_tabulated_fdr-0.05.txt", quote=F, sep="\t", row.names = F)
