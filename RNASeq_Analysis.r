library("cummeRbund")
cuff_data <- readCufflinks("C:/Users/HP/OneDrive/Desktop/diff_out/diff_out")
cuff_data
dens <- csDensity(genes(cuff_data));
dens
dend<-csDendro(genes(cuff_data));
dend
m <- MAplot(genes(cuff_data),"H1hesc","Cd20");
m
V <- csVolcanoMatrix(genes(cuff_data));
V

c<-csScatter(genes(cuff_data), 'H1hesc', 'Cd20')
c
gene_diff_data <- diffData(genes(cuff_data))
gene_diff_data
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
sig_gene_data
nrow(sig_gene_data)
up_gene_data <- subset(sig_gene_data, (log2_fold_change > 1))
up_gene_data 
down_gene_data <- subset(sig_gene_data, (log2_fold_change < -1))
down_gene_data 
geneids <- c(up_gene_data$gene_id)
geneids
myGenes <- getGenes(cuff_data, geneids)
myGenes
csHeatmap(myGenes, cluster="both")
mygene1 <- getGene(cuff_data,'MADCAM1')
mygene1
expressionBarplot(mygene1)
expressionBarplot(isoforms(mygene1))
csVolcano(genes(cuff_data), "H1hesc", "Cd20");


