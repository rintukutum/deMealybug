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
## extract gene level info from the ballgown object
gene.expr <- gexpr(bg_dr)
## limma
model.y <- limma::lmFit(gene.expr,design)
#  Moderated t-statistic
fitBayes_model.y <- limma::eBayes(model.y)
top_model.y <- limma::topTable(
  fitBayes_model.y,
  coef=2,
  number = nrow(gene.expr),
  adjust.method = 'BH'
)
## add the comparison
top_model.y$comp <- paste0(rev(levels(p.df$sex)),collapse = '/')
top_model.y$geneIDs <- rownames(top_model.y)
##
## compute q-value 
top_model.y$qvalue <- qvalue::qvalue(p=top_model.y$P.Value)$qvalues
dir.create('./data',showWarnings = FALSE)
write.csv(
  top_model.y,
  './data/deMealybug-limma-output-gene-level.csv',
  row.names = FALSE
)