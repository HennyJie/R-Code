

### load in library and data
library(monocle)
library(HSMMSingleCell)
library(reshape2)


HSMM_expr_matrix = exp(datDCG)
# HSMM_sample_sheet = data.frame(Group = c(rep("ES",48),rep("MFE",44)) )
HSMM_sample_sheet = data.frame(Group = c(rep("Cancer",430),rep("Cell Line",102)) )
HSMM_sample_sheet = data.frame(Group = c(rep("MGH26",118),rep("MGH28",94),rep("MGH29",75),rep("MGH30",73),rep("MGH31",70)
                                         ,rep("Cell Line 1",44),rep("Cell Line 2",58)) )


rownames(HSMM_sample_sheet) = colnames(datDCG)
HSMM_gene_annotation = data.frame(gene_short_name = rownames(datDCG))
rownames(HSMM_gene_annotation) = rownames(datDCG)

## Construct CellDataSet object
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)

HSMM <- estimateSizeFactors(HSMM)
# HSMM <- estimateDispersions(HSMM)

# ## Detect genes - this is a filtering step to select genes for downstream clustering
# HSMM <- detectGenes(HSMM, min_expr = 0.1)
# head(fData(HSMM))
# 
# expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))
# length(expressed_genes)
# head(pData(HSMM)) ## The last column can be used to filter cells

#############################
## cell clustering
#############################
HSMM <- clusterCells(HSMM, num_clusters=7)

pdf(file="H:MFEES.pdf",width = 5, height =5 )
plot_cell_trajectory(HSMM, 1, 2, color="Group",cell_size=3)
dev.off()
