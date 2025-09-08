### get genes for specific chromosomes and plot heatmaps
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(dplyr)
# BiocManager::install("biomaRt")

chr_list <- paste0("Chr", c(1:22, "X", "Y")) # list chromosomes 1-22, x and y
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # get ensembl database

# searchAttributes(mart = ensembl, pattern = "synonym")


genes <- getBM(attributes = c('external_gene_name' , 'hgnc_symbol', "external_synonym", 'ensembl_transcript_id', 'ensembl_gene_id', 
                              'chromosome_name', 'start_position', 'end_position', 'band'), mart = ensembl) # get genes and locations
genes$chromosome_name <- paste0("Chr", genes$chromosome_name) # add "Chr" to chromosome numbers
genes$Arm <- substr(genes$band, 1, 1) # create "Arm" column (p or q) from cytobands
genes$hgnc_symbol <- ifelse(genes$hgnc_symbol == "", genes$external_gene_name, genes$hgnc_symbol) # fill in some missing values
genesdf <- data.frame(Gene = genes$hgnc_symbol, 
                      Gene_Synonym = genes$external_synonym,
                      Chromosome = genes$chromosome_name, 
                      Arm = genes$Arm, 
                      Cytoband = paste0(gsub("Chr", "", genes$chromosome_name), genes$band),
                      Start = genes$start_position) # get columns we want to use
genesdf <- unique(genesdf[ genesdf$Chromosome %in% chr_list,]) # get rid of duplicated rows
genesdf <- genesdf[-which(genesdf$Gene == ""), ] # get rid of missing genes
genesdf$Chromosome <- factor(genesdf$Chromosome, levels = chr_list) # set chr to factors

# grep("-", genesdf$Gene, value = T)
genesdf <- genesdf[order(genesdf$Chromosome, genesdf$Start),] # order rows by chr and start

#set exp_table and metadata, make sure to size normalize
exp_table<-read.csv("~/Desktop/karan_batch2_counts_mean.csv",row.names = 1)
meta<-read.csv("~/Desktop/Karanmetadata_batch2_mean.csv",row.names=1)

### obtain genes from specific chromosome
# get_chr_specific_genes(chr = "Chr9")
get_chr_specific_genes <- function(chr) {
  chrdf <- genesdf[genesdf$Chromosome %in% chr, ] # filter genes df for specific chromosome
  chrdf <- unique(rbind( chrdf[,-2], dplyr::rename(chrdf[,-1], Gene=Gene_Synonym) )) # add synonym names 
  chrdf <- chrdf[!chrdf$Gene == "", ] # remove rows with missing name
  
  chrdf <- chrdf[order(chrdf$Start, decreasing = F),] # order by start site
  
  chrdf$Arm <- factor(chrdf$Arm, levels = c("p", "q"))
  chrdf$Cytoband_simplified <- gsub("\\..*", "", chrdf$Cytoband) # remove sub bands 
  chrdf$Cytoband <- factor(chrdf$Cytoband, levels = unique(chrdf$Cytoband))
  chrdf$Cytoband_simplified <- factor(chrdf$Cytoband_simplified, levels = unique(chrdf$Cytoband_simplified))
  
  rownames(chrdf) <- NULL
  
  return(chrdf)
}

#run the above function
chr_specific_table<-get_chr_specific_genes("Chr21")

####### get locations of specific genes
gene_locations <- function(genes) {
  
  df <- genesdf[genesdf$Gene %in% genes | genesdf$Gene_Synonym %in% genes, ] # filter genes df for genes of interest
  df <- unique(rbind( df[,-2], dplyr::rename(df[,-1], Gene=Gene_Synonym) )) # add synonym names 
  df <- df[!df$Gene == "", ] # remove rows with missing name
  df <- df[df$Gene %in% genes, ] # make sure only names from original gene list match
  df <- df[order(df$Chromosome, df$Start, decreasing = F),] # order by start site
  df$Arm <- factor(df$Arm, levels = c("p", "q"))
  df$Cytoband_simplified <- gsub("\\..*", "", df$Cytoband) # remove sub bands 
  df$Cytoband <- factor(df$Cytoband, levels = unique(df$Cytoband))
  df$Cytoband_simplified <- factor(df$Cytoband_simplified, levels = unique(df$Cytoband_simplified))
  rownames(df) <- NULL
  
  ### check if any genes not found
  if ( any(!genes %in% df$Gene) ) {
    print( "Genes not found: " )
    print(genes[!genes %in% df$Gene])
  }
  
  return(df)
}

#hscr genes human nomenclature, for above function
hscr_genes<-gene_locations(c('RET','GDNF','GFRA1','NRTN','SOX10','EDNRB','EDN3','ECE1','ZEB2','PHOX2B','KIFBP','TCF4','L1CAM','ELP1','SEMA3C','SEMA3D','NRG1','ACSS2','ADAMTS17','ENO3','PRXL2A','SH3PXD2A','SLC27A4','UBR4'))


#to set colors for plot below
col_fun = colorRamp2(c(0, 2, 8), c("blue", "green", "red"))

#how to run function below
plot_chr_specific_heatmap(exp_table= exp_table, datatype = "RNA Z-score", scale="z-score",meta = meta, chr_specific_table = chr_specific_table, chromosome = "21",group_colors = list(foo = col_fun),path="~/Desktop/")

### plot heatmap for different omics for genes on specific chr or specific genes
plot_chr_specific_heatmap <- function(exp_table, datatype, scale="z-score", meta, chr_specific_table, 
                                      chromosome, group_colors, path, h=8, w=10, 
                                      cluster_rows=F, clustcols=F, 
                                      show_column_names=F, show_row_names=F) {
  ### get genes
  samegenes <- unique(intersect(rownames(exp_table), chr_specific_table$Gene))
  
  chr_specific_table_filt <- chr_specific_table[chr_specific_table$Gene %in% samegenes,] # filter genes for only ones in table and gene list
  chr_specific_table_filt <- chr_specific_table_filt[!duplicated(chr_specific_table_filt$Gene),] # get rid of dup
  chr_specific_table_filt <- chr_specific_table_filt[order(chr_specific_table_filt$Start, decreasing = F),] # make sure genes ordered
  
  exp_table_filt <- exp_table[chr_specific_table_filt$Gene,] # filter counts table
  
  ### input tables
  if (scale=="z-score") {
    mat <- t(scale(t(exp_table_filt))) # z-score input matrix
  } else { 
    mat <- as.matrix(exp_table_filt) }
  
  sideannotdf <- data.frame(row.names = chr_specific_table_filt$Gene, Arm = chr_specific_table_filt$Arm) # side annot for hm
  topannotdf <- meta # top annot for hm
  
  if(all(rownames(mat) == rownames(sideannotdf)) == F) { stop("Check genes") }
  if(all(colnames(mat) == rownames(topannotdf)) == F) { stop("Check samples") }
  
  ### set annotations
  splitcyto = chr_specific_table_filt$Cytoband_simplified # split rows by cytobands
  cytotext <- as.list(unique(splitcyto)) # add cytoband names to annotations
  names(cytotext) <- unique(splitcyto)
  
  # side_annot = rowAnnotation(df = sideannotdf, col = list(Arm = c("p" = "#C3C3C3", "q" = "#898989"))) # set side colors
  side_annot = rowAnnotation(textbox = anno_textbox(splitcyto, cytotext, side = "left", by = "anno_block", 
                                                    gp = gpar(col = "black", fontsize=9)))
  
  top_annot = HeatmapAnnotation(df = topannotdf, col = group_colors) # set top colors
  
  ### plot hm
  if (ncol(mat) > 10) { hmw=3 } else { hmw=10 }
  hm <- Heatmap(mat, name = datatype, 
                cluster_rows = cluster_rows, cluster_columns = clustcols, 
                show_column_names = show_column_names, show_row_names = show_row_names, 
                column_names_side = "bottom", row_names_side = "right",
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 8),
                row_split = splitcyto, row_title = NULL,
                column_split = topannotdf$Group, column_title = NULL,
                col = colorRamp2(c(-2, 0, 2), c("#0c23b5", "white", "#9f0909")),
                top_annotation = top_annot, 
                left_annotation = side_annot,
                width = ncol(mat)*unit(hmw, "mm"))
  pdf(paste0(path, chromosome, "_Genes_Heatmap.pdf"), height = h, width = w)
  print(hm)
  junk <- dev.off() 
}


###for genes###
plot_gene_specific_heatmap <- function(exp_table, datatype, scale="z-score", meta, chr_specific_table, 
                                      genes, group_colors, path, h=8, w=10, 
                                      cluster_rows=T, clustcols=F, 
                                      show_column_names=F, show_row_names=F) {
  ### get genes
  samegenes <- unique(intersect(rownames(exp_table), chr_specific_table$Gene))
  
  chr_specific_table_filt <- chr_specific_table[chr_specific_table$Gene %in% samegenes,] # filter genes for only ones in table and gene list
  chr_specific_table_filt <- chr_specific_table_filt[!duplicated(chr_specific_table_filt$Gene),] # get rid of dup
  chr_specific_table_filt <- chr_specific_table_filt[order(chr_specific_table_filt$Start, decreasing = F),] # make sure genes ordered
  
  exp_table_filt <- exp_table[chr_specific_table_filt$Gene,] # filter counts table
  
  ### input tables
  if (scale=="z-score") {
    mat <- t(scale(t(exp_table_filt))) # z-score input matrix
  } else { 
    mat <- as.matrix(exp_table_filt) }
  
  sideannotdf <- data.frame(row.names = chr_specific_table_filt$Gene, gene_name = chr_specific_table_filt$Gene) # side annot for hm
  topannotdf <- meta # top annot for hm
  
  if(all(rownames(mat) == rownames(sideannotdf)) == F) { stop("Check genes") }
  if(all(colnames(mat) == rownames(topannotdf)) == F) { stop("Check samples") }
  
  ### set annotations
  splitcyto = chr_specific_table_filt$Gene 
  cytotext <- as.list(unique(splitcyto)) 
  names(cytotext) <- unique(splitcyto)
  
  side_annot = rowAnnotation(foo = anno_text(splitcyto, cytotext, location = grid::unit(0.5, "npc"), just = "center"))
  top_annot = HeatmapAnnotation(df = topannotdf, col = group_colors)
  
  ### plot hm
  if (ncol(mat) > 10) { hmw=3 } else { hmw=10 }
  hm <- Heatmap(mat, name = datatype, 
                cluster_rows = cluster_rows, cluster_columns = clustcols, 
                show_column_names = show_column_names, show_row_names = show_row_names, 
                column_names_side = "bottom", row_names_side = "right",
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 8),
                column_split = topannotdf$Group, column_title = NULL,
                col = colorRamp2(c(-2, 0, 2), c("#0c23b5", "white", "#9f0909")),
                top_annotation = top_annot, 
                left_annotation = side_annot,
                width = ncol(mat)*unit(hmw, "mm"))
  pdf(paste0(path,genes,"_Genes_Heatmap.pdf"), height = h, width = w)
  print(hm)
  junk <- dev.off() 
}

#to filter genes with 0 variance
gene_variances <- apply(exp_table_filt, 1, var)
zero_variance_genes <- names(gene_variances[gene_variances == 0])
exp_table_filt <- exp_table_filt[!rownames(exp_table_filt) %in% zero_variance_genes, ]

