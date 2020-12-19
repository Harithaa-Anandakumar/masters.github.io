##making pseudoBulk
#figure 4.13
dat <- read.table("data/combined_d15-d30.txt", 
                  header = TRUE,
                  check.names = FALSE)
expDat <- read.table("data/new.txt", 
                  header = TRUE,
                  check.names = FALSE)
library(SCDC)
library(scater)
dat <- as.data.frame(dat)
dim(dat)
genes <-dat$GeneSymbol
#dat <- dat[,2:6393]

rownames(expDat) <- genes
head(dat)
fdata <- rownames(dat)
pdata <- colnames(dat)
eset <- getESET(dat, fdata = fdata, pdata = pdata)
generateBulk_norep(eset$pdata,
                   ct.varname = "cluster", sample = "sample", ct.sub = c("d15:S1","d15:S2","d30:S1","d30:S2"), nbulk = 3)
?generateBulk_norep
qc.seger <- readRDS("data/qc_segerstolpe.rds")
qc.seger@phenoData@data
pseudo.seger.rand <- generateBulk_norep(qc.seger$sc.eset.qc, ct.varname = "cluster", sample = "sample", ct.sub = c("alpha","beta","delta","gamma"), nbulk = 3)
round(pseudo.seger.rand$true_p,3)
######################################################
#' Construct Pseudo bulk samples
#' @description Construct Pseudo bulk samples by actual number of cells per subject
#' @name generateBulk_allcells
#' @param eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @return pseudo bulk samples ExpressionSet, and actual cell-type proportions
#' @export
#' Then create the pseudo bulk samples as follows.
#'  ct.varname specifies the variable name of your single cell clustering result variable. 
#'  sample specifies the variable indicating subjects or individuals 
#'  information.
#'  ct.sub specifies the subset of cell types you want to use to construct pseudo bulk samples.
generateBulk_allcells <- function(eset, ct.varname, sample, disease = NULL, ct.sub = NULL){
  if (is.null(ct.sub)){
    ct.sub <- unique(eset@phenoData@data[,ct.varname])
  }
  eset <- eset[, eset@phenoData@data[,ct.varname] %in% ct.sub]
  cluster.id <- eset@phenoData@data[,ct.varname]
  sample.id <- eset@phenoData@data[,sample]
  condition.id <- eset@phenoData@data[,disease]
  
  ## expression
  pseudo.exprs <- sapply(unique(sample.id), function(sid){
    y <- exprs(eset)[, sample.id %in% sid]
    rowSums(y, na.rm = T)
  })
  colnames(pseudo.exprs) <- unique(sample.id)
  ## true proportion: sample by cell types
  ncount <- table(sample.id, cluster.id)
  true.prop <- ncount / rowSums(ncount, na.rm = T)
  true.prop <- true.prop[complete.cases(true.prop),]
  ## eset for pseudo bulk sample
  if (is.null(disease)){
    pseudo.disease <- NA
  } else {
    pseudo.disease <- sapply(unique(sample.id), function(sid){
      condition.id[sample.id == sid][1]
    })
  }
  pseudo.pdata <- data.frame(sample = colnames(pseudo.exprs),
                             disease = pseudo.disease)
  pseudo.fdata <- data.frame(genes = rownames(pseudo.exprs))
  rownames(pseudo.fdata) <- rownames(pseudo.exprs)
  pseudo_eset <- getESET(exprs = pseudo.exprs,
                         fdata = pseudo.fdata,
                         pdata = pseudo.pdata)
  return(list(truep = true.prop, pseudo_eset = pseudo_eset))
}

## expression
pseudo.exprs <-  sapply(unique(colnames(expDat)), function(sid){
  y <- exprs(eset)[, colnames(expDat) %in% sid ]
  rowSums(y, na.rm = T)
})


trialPseudoBulk <- pseudo.exprs %>% 
  as_tibble(rownames = 'Genes') %>% 
 mutate(equalSplit=rowSums(.[2:5])/4,
        uneqSpOne = (`d30:S2`*10/100) + (`d30:S1`*70/100) + 
          (`d15:S1`*5/100) + (`d15:S2`*15/100),
        uneqSpTwo = (`d30:S2`*80/100)  + (`d15:S2`*20/100))


trialPseudoBulk <- trialPseudoBulk[,c(1,6,7,8)]
(trialPseudoBulk)
write.csv(trialPseudoBulk, "trialPb2.csv", quote = F, row.names = F)
t