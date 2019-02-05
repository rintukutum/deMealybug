rm(list=ls())
library(ballgown)
sampleFolders <- list.files(
  path = '../28jan19/',
  pattern = 'MBVB',
  full.names=TRUE
)
bg_dr<- ballgown(samples = sampleFolders)
pheno.mb <- read.csv('../28jan19/mealybug_phenodata.csv')
pData(bg_dr) <- pheno.mb
p.df <- data.frame(sex = pData(bg_dr)[,'sex'])
rownames(p.df) <- rownames(pData(bg_dr))
design <- model.matrix(~.,p.df)
## extract only the FPKM from the ballgown object
fpkm.expr <- bg_dr@expr$trans[,grep('FPKM',colnames(bg_dr@expr$trans))]
colnames(fpkm.expr) <- gsub('FPKM.','',colnames(fpkm.expr))
## limma
model.y <- limma::lmFit(fpkm.expr,design)
#  Moderated t-statistic
fitBayes_model.y <- limma::eBayes(model.y)
top_model.y <- limma::topTable(
  fitBayes_model.y,
  coef=2,
  number = nrow(fpkm.expr),
  adjust.method = 'BH'
)
## add the comparison
top_model.y$comp <- paste0(rev(levels(p.df$sex)),collapse = '/')
## Get gene ID information
mstrg.geneIdx <- geneIDs(bg_dr)
top_model.y$geneIDs <- mstrg.geneIdx[rownames(top_model.y)]
##
## compute q-value 
top_model.y$qvalue <- qvalue::qvalue(p=top_model.y$P.Value)$qvalues
dir.create('./data',showWarnings = FALSE)
write.csv(
  top_model.y,
  './data/deMealybug-limma-output.csv',
  row.names = FALSE
)