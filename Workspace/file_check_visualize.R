library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyverse)

run <- "A01861_0063"
df_file_check <- read.delim(paste(run, "_output/File_Check_", run, sep = ''), sep = '\t')
#df_file_check <- df_file_check[,c('X0','X1','X2', 'FASTQ', 'fastp', "bwa_decontaminate", "preiVar_sars_alignmt","postiVar_sars_alignmt","freyja_variant_outputs","freya_final_output","freya_output_meets_mincriteria")]
df_file_check <- df_file_check[,3:13]
colnames(df_file_check) <- c('Sample','Run','Run_Name', 'FASTQ', 'fastp', "bwa_decontaminate", "preiVar_sars_alignmt","postiVar_sars_alignmt","freyja_variant_outputs","freya_final_output","freya_output_meets_mincriteria")
df_file_check <- as.data.frame(df_file_check)
df_file_check$count <- apply(df_file_check, 1, function(x) length(which(x == 'False')))
df_file_check <- df_file_check[order(df_file_check$count), ]
df_file_check <- df_file_check[,1:11]
df_file_check_long <- gather(df_file_check, stage, status, FASTQ:fastp:bwa_decontaminate:preiVar_sars_alignmt:postiVar_sars_alignmt:freyja_variant_outputs:freya_final_output:freya_output_meets_mincriteria, factor_key=TRUE)
df_file_check_long$status <- gsub('False', 'Fail', df_file_check_long$status)
df_file_check_long$status <- gsub('True', 'Pass', df_file_check_long$status)
df_file_check_long$status <- as.factor(df_file_check_long$status)
simple_variant_colours <- c('Pass' = "#31a354",'Fail' = "#e34a33")
for (run in unique(sort(df_file_check_long$Run))){
    file_status <- df_file_check_long %>% 
    filter(Run == run) %>%
    ggplot(aes(stage, Sample, fill = status)) + 
    geom_tile(color = "black", lwd = 0.5, linetype = 1) + 
    theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)) +
    scale_fill_manual(values = simple_variant_colours)
  
    ggsave(paste(run, 'file_status.png',sep = '_'), plot = file_status, device = 'png', path = paste("/Users/Sutcliffe/Work/06_CentrEau/02_Analysis/06_Quality_Control/", run, "_output", sep = ''), width = 15 , height =17, units = "cm", limitsize = FALSE)
  
}


