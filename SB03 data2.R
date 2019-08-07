#################
## libraries ####
#################
# You might need to install these libraries and their dependencies if you dont have them yet
library(scater)
library(Seurat)
library(org.Mm.eg.db)
library(scran)
library(TxDb.Mmusculus.UCSC.mm10.ensGene) 
############
## data ####
############
direct <- "~/Desktop/"
# 
hash <- read.csv(file.path(direct, "SB03_raw_data_and_annotations",
                           "/Hash_count_table.csv"), stringsAsFactors=FALSE) 
metadata <-   read.csv(file.path(direct, "SB03_raw_data_and_annotations",
                             "/SB03_annotations.csv"), stringsAsFactors=FALSE) 
h5closeAll()
rhdf5::h5ls(file.path(direct, "SB03_raw_data_and_annotations",
                      "raw_gene_bc_matrices_h5.h5"))

counts <- Seurat::Read10X_h5(file.path(direct, "SB03_raw_data_and_annotations",
                                       "raw_gene_bc_matrices_h5.h5"))


sc_raw <- SingleCellExperiment(list(counts=counts), rowData=counts@Dimnames[[1]])

####

################
## analysis ####
################
## process data ####
#etadata <- read.delim("~/Desktop/SB03_raw_data_and_annotations/SB03_annotations.csv",
#                       check.names = FALSE, header = TRUE, sep=",") 
#raw_index <- gsub("-1", "", SB03_annotations$index) 
#SB03_annotations$index <- raw_index
##Matching for Metadata##
colnames(sc_raw) <- gsub("-1", "", colnames(sc_raw)) 
metadata$index <- gsub("-1", "", metadata$index) 

sc <- sc_raw[,colnames(sc_raw) %in% metadata$index]
m <- match(colnames(sc), metadata$index)       
metadata_sorted <- metadata[m,]    
head(metadata_sorted)  
metadata(sc) <- metadata_sorted 

##Finding ENSEMBL IDs of the genes##
sym <- mapIds(org.Mm.eg.db, keys = rownames(sc), keytype = "SYMBOL",
              column = "ENSEMBL") 
rowData(sc)$SYMBOL <- rownames(sc) 
rowData(sc)$ENSEMBL <- sym
head(rowData(sc)) 
rownames(sc) <- uniquifyFeatureNames(rowData(sc)$SYMBOL, rowData(sc)$ENSEMBL)  
head(rownames(sc))

##Mapping the IDs##
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sc)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID") 
rowData(sc)$CHR <- location 
summary(location=="ChrM") 
mito <- which(rowData(sc)$CHR=="chrM") 
sc <- calculateQCMetrics(sc, feature_controls = list(Mt=mito)) 
head(colnames(colData(sc)), 10) 

## Quality Control ####

## Violin plot##
plotColData(sc, y="n_counts", x="Cell.Type")  
##Scatterplot Graph##
par(mfrow=c(1,3)) 
plot(sc$total_features_by_counts, sc$total_counts, xlab="Number of expressed genes",
     ylab="Library size") 
fontsize <- theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12)) 
##Highest Expression Graph##
plotHighestExprs(sc, n=50) + fontsize 
ave.counts <- calcAverage(sc, use_size_factors=FALSE) 
hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count")) 
demo.keep <- ave.counts >= -4 
demo.keep <-  1 >= demo.keep   
filtered.sc <- sc[demo.keep,] 
summary(demo.keep) 
num.cells <- nexprs(sc, byrow=TRUE) 
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count")) 
to.keep <- num.cells > 0 
sc <- sc[to.keep,] 
summary(to.keep) 
sc <- computeSumFactors(sc) 
summary(sizeFactors(sc)) 
plot(sc$total_counts, sizeFactors(sc), log="xy", xlab="Library size", ylab="Size factor", 
     col=c("red","black")[sc$is_cell_control], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex = 1.2, 
       legend = levels(sc$is_cell_control))
sc <- normalize(sc) 
var.fit <- trendVar(sc, parametric=TRUE, block=sc$t, loess.args=list(span=0.3))  
