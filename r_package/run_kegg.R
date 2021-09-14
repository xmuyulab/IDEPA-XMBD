library(org.Hs.eg.db)
library(clusterProfiler)

options(warn =-1)
args <- commandArgs()

genes_path = args[6]

save_path = args[7]

genes <- read.delim(genes_path, header = F, stringsAsFactors = FALSE)

kegg_output = function(genes, save_path)
{
  colnames(genes)="SYMBOL"
  genes=bitr(genes$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

  kegg <- enrichKEGG(
    gene = genes$ENTREZID,
    keyType = 'kegg',
    organism = 'hsa',
    pAdjustMethod = 'fdr',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
  )
  kegg = setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


  write.csv(kegg@result,
            save_path, 
            sep=',', col.names = T, row.names = F)
  return(kegg@result)
}

kegg_result = kegg_output(genes = genes,
                          save_path = save_path)
