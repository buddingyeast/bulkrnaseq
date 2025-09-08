#MUSIC analysis

#installation may require R>4.3 and some bioconductor packages e.g. SingleCellExperiment
#to be installed separately
devtools::install_version("RefFreeEWAS",version="2.2")
BiocManager::install("TOAST")
BiocManager::install("SingleCellExperiment")
BiocManager::install("Biobase")
devtools::install_github("renozao/xbioc")

#install
BiocManager::install("Biobase")
devtools::install_github('xuranw/MuSiC')
devtools::install_github("Jiaxin-Fan/MuSiC2")

#load stuff
library(MuSiC)
library(MuSiC2)
library(dplyr)
library(Seurat)
library(Biobase)
library(SingleCellExperiment)

#change working directory
setwd('/gpfs/scratch/finer03/')
load("./my_work_space.RData")

#load transcript table
cts<-read.csv("./Ming_Clean_sorted.csv",row.names=1, header=TRUE)

#load metadata
coldat<-read.csv("./Table.S2.sorted.csv")
#coldat <- coldat %>% mutate(genotype=paste0(Ret, "_", Ednrb))

#subset and filter
#coldat_fil<-subset(coldat, rownames(coldat) != c("E145_51","E145_18"))
#cts_fil<-select(cts, -c("E145_51","E145_18"))

#split time
coldat_14<- subset(coldat, age=="E14.5")
cts_14<-cts[,1:28]
coldat_0<- subset(coldat_fil, age=="P0")
cts_0<-cts_fil[,29:58]
cts_fil<-subset(cts_14, rownames(cts_14) %in% genelist$E14.5_Multiplicative_Up)

#data wrangle for ExpressionSet format
ctsmat<-as.matrix(cts_fil)
metadata <- data.frame(labelDescription= c("id","genotype", "sex","age"),row.names=c("id","genotype", "sex","age"))
pheno<-new("AnnotatedDataFrame",data=coldat_fil,varMetadata=metadata)
colnames(ctsmat)<-rownames(coldat_14)
bulkset<-ExpressionSet(assayData = ctsmat, phenoData=pheno)

#split time
ctsmat14<-as.matrix(cts_14)
ctsmat0<-as.matrix(cts_0)
metadata <- data.frame(labelDescription= c("id","genotype", "sex","age"),row.names=c("id","genotype", "sex","age"))
pheno14<-new("AnnotatedDataFrame",data=coldat_14,varMetadata=metadata)
#may be error with row and column names if data not loaded exactly as above, to fix use #colnames(ctsmat)<-rownames(coldat_14)
bulkset14<-ExpressionSet(assayData = ctsmat14, phenoData=pheno14)

#P0
pheno0<-new("AnnotatedDataFrame",data=coldat_0,varMetadata=metadata)
bulkset0<-ExpressionSet(assayData = ctsmat0, phenoData=pheno0)

#check to be true for format to work
all(rownames(coldat_fil)==colnames(ctsmat))

#load scdata
scdata<-readRDS("./E145_filtered_allcells_WTonly.rds")
DefaultAssay(object = scdata)<-"RNA"
scdata[["RNA"]] <- as(scdata[["RNA"]], Class="Assay")

#convert to singlecellexperiment
sce <- as.SingleCellExperiment(scdata)
rm(scdata)

bulk.control.mtx = exprs(bulkset)[, bulkset$genotype == 'WT']
bulk.case.mtx = exprs(bulkset)[, bulkset$genotype == 'Ret_Het_Piebald_Hom']

#run music
#I had to modify my scRNAseq data inputs
#clusters changed to ident which is in metadata
#sample_ID was changed to sample in metadata
#for both variables check the singlecellexperiment object and look under colData listData to find variable names
#also had to modify the object to change log counts to counts with following code:

logcounts<-assay(sce,"logcounts")
counts(sce)<-2^logcounts
assayNames(sce) 
#should see counts as an assay name
set.seed(123)
est = music2_prop_t_statistics(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = sce, clusters = 'ident', samples = 'seq_folder', select.ct = c('Epithelial','IESC','Immune Cells','Neuronal','Enteroendocrine','Progenitor/Glia','Myofibroblasts','Smooth Muscle','Early Differentiating Neurons','Endothelial','Mesenchyme','Enterocytes','Erythrocytes'), n_resample=10, sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)

#faster? Toast analysis
est = my_music2_prop_toast(bulk.control.mtx = bulk.control.mtx, bulk.case.mtx = bulk.case.mtx, sc.sce = sce, clusters = 'cluster', samples = 'sample', select.ct = c('Neurons/Glia','Myofibroblasts','Endothelial','Mesothelial','Smooth Muscle','Immune','Mesenchymal Stem Cells'), prop_r=0.1, cutoff_c=10^(-3), cutoff_r=10^(-3), cap=0.3)

#plots

est.prop = est$Est.prop
prop_all = cbind('proportion'=c(est.prop),'sampleID'=rep(rownames(est.prop),times=ncol(est.prop)),'celltype'=rep(colnames(est.prop), each=nrow(est.prop)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
prop_all$group = ifelse(prop_all$sampleID == c('SRR5973222','SRR5973223','SRR5973224','SRR5973225','SRR5973226','SRR5973227','WTP0_1','WTP0_10','WTP0_9','WTP0_2','WTP0_5','WTP0_6'), 'WT', 'Ret_Het_Piebald_Hom')

write.table(prop_all, "propall.txt", sep=',')

ggplot(prop_all, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#bug in toast, use this function
my_music2_prop_toast <- function (bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, 
                                  samples, select.ct, expr_low = 20, prop_r = 0.1, eps_c = 0.05, 
                                  eps_r = 0.01, cutoff_c = 10^(-3), cutoff_r = 10^(-3), cap = 0.3, 
                                  maxiter = 200, markers = NULL, ct.cov = FALSE, cell_size = NULL,
                                  centered = FALSE, normalize = FALSE) {
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if (length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx), 
                                    nrow(bulk.case.mtx))) {
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, 
  ])
  Pheno = data.frame(condition = factor(c(rep("control", ncol(bulk.control.mtx)), 
                                          rep("case", ncol(bulk.case.mtx))), levels = c("control", 
                                                                                        "case")))
  rownames(Pheno) = colnames(bulk.mtx)
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if (length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))) {
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr >= expr_low])
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)]
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce, 
                            clusters = clusters, samples = samples, select.ct = select.ct, 
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                            iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                            normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case, sc.sce = sc.sce, 
                             clusters = clusters, samples = samples, select.ct = select.ct, 
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                             iter.max = 1000, nu = 1e-04, eps = 0.01, centered = centered, 
                             normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  iter = 1
  ncell = length(select.ct)
  id_conv = NULL
  while (iter <= maxiter) {
    Y_raw = log1p(bulk.mtx)
    design = Pheno
    Prop <- prop_all[rownames(Pheno), ]
    Design_out <- makeDesign(design, Prop)
    fitted_model <- fitModel(Design_out, Y_raw)
    res_table <- csTest(fitted_model, coef = "condition", verbose = F)
    mex = apply(prop_all, 2, mean)
    lr = NULL
    for (celltype in select.ct) {
      m = mex[celltype]
      DE = res_table[[celltype]]
      pval = DE$fdr
      names(pval) = rownames(DE)
      pval = pval[names(pval) %in% exp_genel]
      if (m >= prop_r) {
        lr = c(lr, names(pval[pval <= cutoff_c & pval <= 
                                quantile(pval, prob = cap)]))
      }
      else {
        lr = c(lr, names(pval[pval <= cutoff_r & pval <= 
                                quantile(pval, prob = cap)]))
      }
    }
    lr = unique(lr)
    l = setdiff(gene_all, lr)
    sc.iter.sce = sc.sce[l, ]
    if (length(id_conv) > 0) {
      case_sample = bulk.case[, !colnames(bulk.case) %in% 
                                id_conv]
    }
    else {
      case_sample = bulk.case
    }
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           verbose = F)$Est.prop.weighted
    prop_CASE = rbind(prop_case, prop_case_fix)
    if (length(id_conv) == 1) {
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    prop_all = rbind(prop_control, prop_CASE)
    prop_case = prop_case[rownames(prop_case_ini), ]
    pc = abs(prop_case - prop_case_ini)
    conv = pc
    conv[, ] = 1
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= 
                                                prop_r] < eps_r, 0, 1)
    pc[prop_case_ini > prop_r] = pc[prop_case_ini > prop_r]/prop_case_ini[prop_case_ini > 
                                                                            prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > 
                                               prop_r] < eps_c, 0, 1)
    convf = apply(conv, 1, function(x) {
      all(x == 0)
    })
    all_converge = FALSE
    id_conv = c(id_conv, names(convf[convf == TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% 
                                id_conv, ]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv, 
    ]
    if (is.vector(prop_case_ini)) {
      all_converge = TRUE
      break
    }
    else if (nrow(prop_case_ini) == 0) {
      all_converge = TRUE
      break
    }
    iter = iter + 1
  }
  if (all_converge) {
    return(list(Est.prop = prop_all, convergence = TRUE, 
                n.iter = iter, DE.genes = lr))
  }
  else {
    return(list(Est.prop = prop_all, convergence = FALSE, 
                id.not.converge = rownames(prop_case_ini)))
  }
}

#music analysis
est = music_prop(bulk.mtx = bulk.mtx, sc.sce = sce, clusters = 'ident', samples = 'seq_folder', select.ct = c('Epithelial','IESC','Immune Cells','Neuronal','Enteroendocrine','Progenitor/Glia','Myofibroblasts','Smooth Muscle','Early Differentiating Neurons','Endothelial','Mesenchyme','Enterocytes','Erythrocytes'), verbose=F)

est.prop = est$Est.prop.weighted
prop_all = cbind('proportion'=c(est.prop),'sampleID'=rep(rownames(est.prop),times=ncol(est.prop)),'celltype'=rep(colnames(est.prop), each=nrow(est.prop)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
#rep below is replicate funciton, i.e. repeat x, 13 times to make prop_all$group vector
prop_all$group = rep(coldat_0$genotype, 13)

#Plot
gplot(prop_all, aes(x=celltype, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.5,color='gray36')+
  facet_grid(.~group)+
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12,angle = 45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
write.table(prop_all, "propall.txt", sep=',')
