SB03_annotations <- read.csv("~/Desktop/SB03_raw_data_and_annotations/SB03_annotations.csv")
hash <- read.csv("~/Desktop/SB03_raw_data_and_annotations/Hash_count_table.csv")
H5Fopen("~/Desktop/SB03_raw_data_and_annotations/raw_gene_bc_matrices_h5.h5")
group <- paste(SB03_annotations$hash_ID, SB03_annotations$`Cell Type`, sep = "_")  
group <- factor(group) 
table(group) 
Counts <- read.delim("~/Desktop/SB03_raw_data_and_annotations/SB03_annotations.csv",
                     sep=',')

                     row.names = SB03_annotations$index, 
                     col.names = paste(SB03_annotations$hash_ID, SB03_annotations$`Cell Type`, sep = "_")) 
