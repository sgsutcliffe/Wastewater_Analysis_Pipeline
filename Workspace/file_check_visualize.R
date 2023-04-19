# Read data and select columns
#Version
#file_check_visualize_V1
### Arguments ###
args = commandArgs(trailingOnly=TRUE)
df_file_check <- read.delim(args[1], sep = '\t')
df_file_check <- df_file_check[,2:10]
colnames(df_file_check) <- c('Sample','FASTQ', 'fastp', "bwa_decontaminate", "preiVar_sars_alignmt","postiVar_sars_alignmt","freyja_variant_outputs","freya_final_output","freya_output_meets_mincriteria")
df_file_check <- as.data.frame(df_file_check)

# replace True with 1 and False with 0
df_file_check[, 2:9] <- ifelse(df_file_check[, 2:9] == "True", 1, 0)

# set the first column as row names
rownames(df_file_check) <- df_file_check[,1]
df_file_check <- df_file_check[, -1] # remove the first column

# Open a new graphics device
cairo_pdf(filename = paste(args[2],"/file_check_status-test.pdf", sep = ""), width = 10, height = 13, pointsize = 6)
heatmap(as.matrix(df_file_check),
                       scale="none", 
                       col = c("#e34a33","#31a354"), 
                       Rowv = NA,
                       Colv = NA,
                       xlab = "Stage", ylab = "Sample",
                       main = "File Status of Samples Analyzed in Pipeline",
                       cex.axis = 0.1, cexRow = 0.5,
                       cex.lab = 1,
                       margins = c(5,15))

# Save plot
#plot(file_status)
dev.off()
