library(tidyverse)

#Run "missing_values.R and benjamini_for_sequenom.R" first.

#Goal of this script is to combine values from the two dataframes so that all the information is in a single df. The number of missing values should be in brackets (n). Default seperator of "_" is fine.

#Place brackets around all entries in the missing_values dataframe.
#remove data no longer used at this stage
complete_missing_df_new <- complete_missing_df[c(-3,-4,-5)]

#Add brackets to entries
for (i in 2:ncol(complete_missing_df_new)){
  
  complete_missing_df_new[,i] <- sub("^", "(", complete_missing_df_new[,i])
  complete_missing_df_new[,i] <- sub("$", ")", complete_missing_df_new[,i])

}

#merge the two dfs.
mega_table <- merge(final, complete_missing_df_new, by.x = "SNP", by.y = "SNP")
unified_table <- mega_table
x_columns_indices <- grep("\\.x", colnames(mega_table))

for (j in x_columns_indices){
  x_col_name <- colnames(mega_table)[j]
  y_match <- sub("\\.x", "\\.y", x_col_name)
  unified_name <- sub("\\.x", "", x_col_name)
  unified_table <- unite_(unified_table, paste(unified_name), c(x_col_name, y_match), sep=">")
}

write.table(unified_table, file="~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/unified_table_fdr-0.05.txt", row.names = F, quote = F, sep="\t")
