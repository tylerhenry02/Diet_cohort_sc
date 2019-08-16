#################
## libraries ####
#################
# You might need to install these libraries and their dependencies if you dont have them yet
library(scater)
library(Seurat)
library(org.Mm.eg.db)
library(scran)
library(TxDb.Mmusculus.UCSC.mm10.ensGene) 
library(limma) 
############
## data ####
############
direct <- "~/Desktop/"
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

## Mitochondrial and Ribosomal gene list ####
Mito_ribo <- c('mt-Nd1', 'mt-Nd2', 'mt-Co1', 'mt-Co2', 'mt-Atp8', 'mt-Atp6','mt-Co3','mt-Nd3',
               'mt-Nd4l', 'mt-Nd4','mt-Nd5','mt-Nd6','mt-Cytb','Rpl7','Rpl31','Rpl37a',
              'Rps6kc1','Rpl7a','Rpl12','Rpl35','Rps21','Rpl39','Rpl10','Rps4x',
                'Rps6ka6','Rpl36a','Rps6ka3','Rpl22l1','Rps3a1','Rps27','Rpl34',
                'Rps20','Rps6','Rps8','Rps6ka1','Rpl11','Rpl22','Rpl9','Rpl5','Rplp0',
                'Rpl6','Rpl21','Rpl32','Rps9','Rpl28','Rps5','Rps19','Rps16',
                'Rps11','Rpl13a','Rpl18','Rps17','Rps3','Rpl27a','Rps13',
                'Rps15a','Rplp2','Rps12','Rps15','Rpl6l','Rpl41','Rps26','Rpl18a',
                'Rps18-ps3','Rps26-ps1','Rpl13','Rpl21-ps4','Rpl15','Rps24','Rpl23a-ps3',
                'Rpl13-ps3','Rps2-ps6','Rps25','Rpl10-ps3','Rplp1','Rpl4','Rps27l','Rpl29',
                'Rps27rt','Rpsa','Rpl14','Rps27a','Rpl26','Rpl23a','Rpl9-ps1','Rps6kb1',
                'Rpl23','Rpl19','Rpl27','Rpl38','Rps23','Rpl36-ps3','Rps7','Rpl10l',
                'Rps29','Rpl36al','Rps6kl1','Rps6ka5','Rpl37','Rpl30','Rpl7a-ps3',
                'Rpl8','Rpl3','Rps19bp1','Rpl39l','Rpl35a','Rpl24','Rps6ka2',
                'Rps2','Rpl3l','Rps10','Rpl10a','Rps28','Rps18','Rpl7l1',
                'Rpl36','Rpl7a-ps5','Rpl36-ps4','Rpl27-ps3','Rps14','Rpl17','Rps6kb2',
                'Rps6ka4','Rpl9-ps6','Rpl13a-ps1','Rps12-ps3','Xist')
####
################
## analysis ####
################
## process data ####
#etadata <- read.delim("~/Desktop/SB03_raw_data_and_annotations/SB03_annotations.csv",
#                       check.names = FALSE, header = TRUE, sep=",") 
#raw_index <- gsub("-1", "", SB03_annotations$index) 
#SB03_annotations$index <- raw_index
##Matching for Metadata####
colnames(sc_raw) <- gsub("-1", "", colnames(sc_raw)) 
metadata$index <- gsub("-1", "", metadata$index) 

sc <- sc_raw[,colnames(sc_raw) %in% metadata$index]
m <- match(colnames(sc), metadata$index)       
metadata_sorted <- metadata[m,]    
head(metadata_sorted)  
metadata(sc) <- metadata_sorted 

##Finding ENSEMBL IDs of the genes####
sym <- mapIds(org.Mm.eg.db, keys = rownames(sc), keytype = "SYMBOL",
              column = "ENSEMBL") 
rowData(sc)$SYMBOL <- rownames(sc) 

## Trying to match Mito_ribo genes to genes in dataset ####
length(match(sym, Mito_ribo)) 
matched_mito <- match(Mito_ribo, sym)  
rowData(sc)$ENSEMBL <- sym
head(rowData(sc)) 
rownames(sc) <- uniquifyFeatureNames(rowData(sc)$SYMBOL, rowData(sc)$ENSEMBL)  
head(rownames(sc))

##Mapping the IDs####
ENSEMBL_noNA <- rowData(sc)$ENSEMBL[!is.na(rowData(sc)$ENSEMBL)]
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=ENSEMBL_noNA, 
                   column="CDSCHROM", keytype="GENEID", asNA=TRUE) 
# Checking that location IDs are matched to input IDs
#sum(names(location) == ENSEMBL_noNA)
rowData(sc)$CHR <- NA
rowData(sc)$CHR[!is.na(rowData(sc)$ENSEMBL)] <- location

# manually assigning mitochondrial positions
mito <- which(grepl("mt-", rowData(sc)$SYMBOL))
Ribo <- which(grepl("rpl|rps", rowData(sc)$SYMBOL)) 
rowData(sc)$CHR[grepl("mt-", rowData(sc)$SYMBOL)] <- "chrM"
sc <- calculateQCMetrics(sc, feature_controls = list(Mt=mito, Rb=Ribo)) 
head(colnames(colData(sc)), 10) 

 

##multiplot ####
colData(sc)$Diet <- metadata$hash_ID
multiplot(
  plotColData(sc, y="total_counts", x="Diet", colour_by="Diet"),
  plotColData(sc, y="total_features_by_counts", x="Diet", colour_by = "Diet"),
  plotColData(sc, y="pct_counts_feature_control", x="Diet", colour_by = "Diet"),
  plotColData(sc, y="pct_counts_Mt", x="Diet", colour_by = "Diet"), 
cols = 2)
plotColData(sc, y="pct_counts_Rb", x="Diet", colour_by = "Diet") 
##Multiscatter ####
par(mfrow=c(1,3))
plot(sc$total_features_by_counts, sc$total_counts/1e6, xlab="Number of expressed genes",
     ylab="Library size") 
plot(sc$total_features_by_counts, sc$pct_counts_endogenous, xlab="Number of expressed genes",
     ylab="Endogenous proportion (%)") 
plot(sc$total_features_by_counts, sc$pct_counts_Mt, xlab="Number of expressed genes",
     ylab="Mitochondrial proportion (%)")

## Finding outliers ####
libsize.drop <- isOutlier(sc$total_counts, nmads = 3, type="lower",
                          log = TRUE, batch = sc$DietCell) 
feature.drop <- isOutlier(sc$total_features_by_counts, nmads = 3, type = "lower", 
                          log = TRUE, batch = sc$DietCell)
libsize.drop <- isOutlier(sc$total_counts, nmads = 3, type="higher",
                          log = TRUE, batch = sc$DietCell) 
feature.drop <- isOutlier(sc$total_features_by_counts, nmads = 3, type = "higher", 
                          log = TRUE, batch = sc$DietCell)
keep <- !(libsize.drop | feature.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=sum(keep)) 
sc$PassQC <- keep 
sc <- sc[,keep]
dim(sc) 
attr(libsize.drop, "thresholds") 
attr(feature.drop, "thresholds") 
## assign experimental setup to colData ####
colData(sc)$Diet <- metadata$hash_ID
colData(sc)$Cell_type <- metadata$Cell.Type
colData(sc)$index <- metadata$index
colData(sc)$n_counts <- metadata$n_counts
sc$pct
## Quality Control ####
## Violin plot####
plotColData(sc, y="total_counts", x="Diet")   
plotColData(sc[,sc$total_counts > 8000 & sc$total_counts < 22000],
            y="total_counts", x="Diet")   

##Scatterplot Graph####
par(mfrow=c(1,3)) 
plot(sc$total_features_by_counts, sc$total_counts, xlab="Number of expressed genes",
     ylab="Library size") 
fontsize <- theme(axis.text = element_text(size = 8), axis.title = element_text(size = 12)) 
##Highest Expression Graph##
plotHighestExprs(sc, n=50) + fontsize 
ave.counts <- calcAverage(sc, use_size_factors=FALSE) 

## Histogram ####
hist(sc[,sc$total_counts > 2500 & sc$total_counts < 22000], breaks=100, main="",
     col="grey80",
     xlab=expression("total_counts")) 
filtered.sc <- sc[demo.keep,] 
summary(demo.keep) 
num.cells <- nexprs(sc, byrow=TRUE) 

## Scatter Graph ####
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count")) 
to.keep <- num.cells > 0 
sc <- sc[to.keep,] 
summary(to.keep) 
sc <- computeSumFactors(sc) 
summary(sizeFactors(sc)) 

## Normalization ####
sc <- normalize(sc)

##Variance Graph ####
var.fit <- trendVar(sc, parametric=TRUE, block=sc$t, loess.args=list(span=0.3), 
                    use.spikes=FALSE)  
var.out <- decomposeVar(sc, var.fit) 
head(var.out) 
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression") 
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE) 
cur.spike <- isSpike(sc) 
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16) 
chosen.genes <- order(var.out$bio, decreasing = TRUE)[1:10] 
plotExpression(sc, features = rownames(var.out)[chosen.genes]) + fontsize 

## Removing Batch Effects ####
assay(sc, "corrected") <- removeBatchEffect(logcounts(sc), 
                                            design = model.matrix(~metadata$hash_ID),
                                            batch = metadata$Cell.Type) 
assayNames(sc) 
sc <- denoisePCA(sc, technical=var.out, assay.type="corrected") 
dim(reducedDim(sc, "PCA")) 
plotReducedDim(sc, use_dimred="PCA", ncomponents=3, colour_by="Diet") + fontsize

##TSNE plots #### 
set.seed(100) 
out5 <- plotTSNE(sc, run_args = list(use_dimred="PCA", perplexity=5),
                 colour_by="Cell_type") + fontsize + ggtitle("Perplexity = 5") 
set.seed(100) 
out10 <- plotTSNE(sc, run_args = list(use_dimred="PCA", perplexity=10), 
                  colour_by="Cell_type") + fontsize + ggtitle("Perplexity = 10") 
set.seed(100) 
out20 <- plotTSNE(sc, run_args = list(use_dimred="PCA", perplexity=20), 
                  colour_by="Cell_type") + fontsize + ggtitle("Perplexity = 20") 
multiplot(out5, out10, out20, cols=3) 
set.seed(100) 
sc <- runTSNE(sc, use_dimred = "PCA", perplexity = 20) 
reducedDimNames(sc) 
pcs <- reducedDim(sc, "PCA")

##CLusters ####
my.dist <- dist(pcs) 
my.tree <- hclust(my.dist, method = "ward.D2") 
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
                                    minClusterSize = 10, verbose=0)) 
table(my.clusters, sc$Cell_type) 
table(my.clusters, sc$Diet)
sc$cluster <- factor(my.clusters) 
plotTSNE(sc, colour_by="cluster") + fontsize
clust.col <- scater:::.get_palette("tableau10medium") 
sil <- silhouette(my.clusters, dist = my.dist) 
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])] 
sil.cols <- sil.cols[order(-sil[,1], sil[,3])] 
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), border=sil.cols, 
     col=sil.cols, do.col.sort=FALSE) 
markers <- findMarkers(sc, my.clusters, block=sc$index)  
marker.set <- markers[["1"]] 
head(marker.set, 10) 
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sc, features=top.markers, columns=order(sc$cluster), 
            colour_columns_by=c("cluster", "Cell_type", "Diet"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5,5))