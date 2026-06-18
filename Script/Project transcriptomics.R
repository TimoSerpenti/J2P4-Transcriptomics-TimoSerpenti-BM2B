getwd()
install.packages("BiocManager")
BiocManager::install('Rsubread')
library(Rsubread)
buildindex(
  basename = 'Reuma',
  reference = 'Reuma/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 8000,
  indexSplit = TRUE)
align.reuma1 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785979_2_subset40k.fastq", output_file = "reuma1.BAM")
align.reuma2 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "reuma2.BAM")
align.reuma3 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785986_2_subset40k.fastq", output_file = "reuma3.BAM")
align.reuma3 <- align(index = "Reuma", readfile1 ="Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785988_2_subset40k.fastq", output_file = "reuma4.BAM")
align.rctrl1 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "rctrl1.BAM")
align.rctrl2 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785820_2_subset40k.fastq", output_file = "rctrl2.BAM")
align.rctrl3 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "rctrl3.BAM")
align.rctrl4 <- align(index = "Reuma", readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2= "Data_RA_raw/SRR4785831_2_subset40k.fastq", output_file = "rctrl4.BAM")
BiocManager::install('Rsamtools')
library(Rsamtools)
reumasamples <- c('reuma1', 'reuma2', 'reuma3', 'reuma4', 'rctrl1', 'rctrl2', 'rctrl3', 'rctrl4')
lapply(reumasamples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
lapply(reumasamples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})
allreumasamples <- c('reuma1.BAM', 'reuma2.BAM', 'reuma3.BAM', 'reuma4.BAM', 'rctrl1.BAM', 'rctrl2.BAM', 'rctrl3.BAM', 'rctrl4.BAM')
count_matrixreuma <- featureCounts( files = allreumasamples,
                               annot.ext = "Reuma/ncbi_dataset/data/GCF_000001405.40/genomic.gtf",
                               isPairedEnd = TRUE,
                               isGTFAnnotationFile = TRUE,
                               useMetaFeatures = TRUE)
str(count_matrixreuma)
countsreuma <- count_matrixreuma$counts
head(countsreuma)
colnames(countsreuma) <- c('reuma1', 'reuma2', 'reuma3', 'reuma4', 'rctrl1', 'rctrl2', 'rctrl3', 'rctrl4')
head(countsreuma)
write.csv(countsreuma, "reuma_countmatrix.csv")

countsreuma <- read.csv("reuma_countmatrix.csv", row.names = 1)
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")
library("DESeq2")
library("KEGGREST")
library(EnhancedVolcano)
library(pathview)
samples <- c("reuma", "reuma", "reuma","reuma", "control", "control", "control", "control")
samples_table <- data.frame(samples)
rownames(samples_table) <- c('reuma1', 'reuma2', 'reuma3', 'reuma4', 'rctrl1', 'rctrl2', 'rctrl3', 'rctrl4')
head(samples_table)

reuma_table<- read.table("count_matrix_RA.txt", row.names= 1)
completesamples <- c("control", "control", "control", "control","reuma", "reuma", "reuma","reuma")
complete_table <- data.frame(completesamples)
rownames(samples_table) <- c('rctrl1', 'rctrl2', 'rctrl3', 'rctrl4','reuma1', 'reuma2', 'reuma3', 'reuma4')
colnames(reuma_table) <- c( 'rctrl1', 'rctrl2', 'rctrl3', 'rctrl4','reuma1', 'reuma2', 'reuma3', 'reuma4')
head(complete_table)
ddsreuma <- DESeqDataSetFromMatrix(countData = reuma_table,
                              colData = complete_table,
                              design = ~ completesamples)




ddsreuma <- DESeq(ddsreuma)
reumaresultaten <- results(ddsreuma)
write.table(reumaresultaten, file = 'ResultatenReuma.csv', row.names = TRUE, col.names = TRUE)
sum(reumaresultaten$padj < 0.05 & reumaresultaten$log2FoldChange > 1, na.rm = TRUE)
sum(reumaresultaten$padj < 0.05 & reumaresultaten$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- reumaresultaten[order(reumaresultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- reumaresultaten[order(reumaresultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- reumaresultaten[order(reumaresultaten$padj, decreasing = FALSE), ]
EnhancedVolcano(reumaresultaten,
                lab = rownames(reumaresultaten),
                x = 'log2FoldChange',
                y = 'padj')
EnhancedVolcano(
  reumaresultaten,
  lab = rownames(reumaresultaten),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1,
  labSize = 3
)
dev.copy(png, 'VolcanoplotREUMA2.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

  
  
library(dplyr)

rownames.all = rownames(reumaresultaten)

reumaresultaten = as.data.frame(reumaresultaten)

deg = reumaresultaten %>% filter(padj < 0.05, log2FoldChange > 1)

rownames_deg = rownames(deg)

gene.vector = as.integer(rownames.all %in% rownames_deg)
names(gene.vector) = rownames.all
View(gene.vector)

BiocManager::install("goseq")

library(goseq)
library(geneLenDataBase)
library(org.Hs.eg.db)
library(AnnotationDbi)

PWF = nullp(DEgenes = gene.vector, genome = "hg19", id = "geneSymbol")

GO.wall = goseq(PWF, "hg19", "geneSymbol")

class(GO.wall)
head(GO.wall)
nrow(GO.wall)

enriched.GO = GO.wall$category[GO.wall$over_represented_pvalue<.05]

class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

library(GO.db)
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")}, file="SigGo.txt")
GO.sig <- GO.wall %>% 
  filter(over_represented_pvalue < 0.05)
topGO <- GO.sig %>%
  arrange(over_represented_pvalue) %>%
  head(20)
topGO
library(ggplot2)
ggplot(topGO,
       aes(x = reorder(term, -log10(over_represented_pvalue)),
           y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("GO term") +
  ylab("-log10(p-value)") +
  ggtitle("Top enriched GO terms in Rheumatoid Arthritis")

genes_go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0006954",
  keytype = "GOALL",
  columns = "SYMBOL"
)
go_deg <- intersect(genes_go$SYMBOL, rownames_deg)

go_deg
reumaresultaten[1] <- NULL
reumaresultaten[2:5] <- NULL
ReumaPathway <- reumaresultaten$log2FoldChange
names(ReumaPathway) <- rownames(reumaresultaten)

library(pathview)
pathview(
  gene.data = reumaresultaten,
  pathway.id = "hsa04668",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 5))

    