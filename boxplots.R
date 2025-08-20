#vln and box plots samples

#gene plots
plotCounts(dds, gene="ENSMUSG00000022122.15", intgroup="genotype")
x<- plotCounts(dds, gene="ENSMUSG00000022122.15", intgroup="genotype", main = "Ednrb E14.5 Expression", returnData = TRUE)
#above is Ednrb Ensembl stable id

#dot plot
x$name<-coldat$sex
ggplot(x,aes(x=genotype, y=count, color="genotype")) + ggtitle("Ednrb E14.5 Expression")+ geom_point(position=position_jitter(w=.1,h=0))+ scale_y_log10(breaks=c(500,1000,2000,3000))+ geom_text_repel(aes(label=name)) + scale_x_discrete(guide = guide_axis(n.dodge=2))

#split into sexes and gene
library(ggpubr)
x<- plotCounts(dds, gene="Ret", intgroup="genotype", main = "Ret E14.5 Expression", returnData = TRUE)
x$sex<-coldat$sex
p<-ggplot(x, aes(x=sex, y=count)) + geom_point(position=position_jitter(w=.1,h=0))+ ggtitle("Ret E14.5 Expression")+ facet_wrap(~genotype) +theme(axis.text=element_text(size=16),axis.title=element_text(size=18), plot.title = element_text(size=20))
p+ stat_compare_means(method="t.test", label="p.format")

#save
ggsave("~/Desktop/test.pdf",last_plot(), device="pdf", dpi = 300, width = 10, height = 8, units = "cm",scale = 2)

#violin plot,requires x plot counts from above#
x<- plotCounts(dds, gene="Ednrb", intgroup="genotype", returnData = TRUE)
p<-ggplot(x, aes(x=genotype, y=count)) + geom_violin() + ggtitle("Ednrb E14.5 Expression")
p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")

#box plot with error bars
library(ggsignif)
library(ggpubr)

design<-as.data.frame(colData(dds))
myorder<-factor(design$genotype, levels = c('WT','Piebald_Het','Piebald_Hom','Ret_Het_Piebald_Het','Ret_Het_Piebald_Hom'),ordered=TRUE)

x<- plotCounts(dds, gene="Tmeff2", intgroup="genotype", returnData = TRUE)
ggplot(x, aes(x=myorder, y=count)) + geom_boxplot() +  ggtitle("Tmeff2")+ geom_signif(comparisons = list(c("Ret_Het_Piebald_Hom","WT"),c("Ret_Het_Piebald_Het","WT"),c("Piebald_Hom","WT"),c("Piebald_Het","WT")), map_signif_level=TRUE, y_position=c(4000,4400,4900,5300))+ scale_x_discrete(guide = guide_axis(n.dodge=2))+theme(axis.text=element_text(size=16),axis.title=element_text(size=18), plot.title = element_text(size=20))+ ylim(0,5800)

#box plot with anova test
ggplot(x, aes(x=genotype, y=count)) + geom_boxplot() +  ggtitle("Ednrb E14.5 Expression")+  stat_compare_means(method = "anova", label.y = 5)


###box plots of all hscr genes###
library(tidyr)
library(tibble)
library(ggsignif)
library(dplyr)

#load transcript table, note: gene column needs to be labelled "gene_names" to work
cts<-read.csv("~/Desktop/Ming_Clean.csv",row.names=1, header=TRUE)
coldat<-read.csv("~/Desktop/Table.S2.csv")
dds <- DESeqDataSetFromMatrix(countData=round(cts),colData=coldat,design=~genotype)

#normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
ncts<-as.data.frame(normalized_counts)
ncts <- tibble::rownames_to_column(ncts, "gene_names")

hscr_genes<-c('Ret','Gdnf','Gfra1','Nrtn','Sox10','Ednrb','Edn3','Ece1','Zeb2','Phox2b','Kifbp','Tcf4','L1cam','Elp1','Sema3c','Sema3d','Nrg1','Acss2','Adamts17','Eno3','Prxl2a','Sh3pxd2a','Slc27a4','Ubr4')
hscr_genes<-c('Gart','Sumo3','Hmgn1','Son')

hscr_gene_counts<-subset(ncts, gene_names %in% hscr_genes)

gather_hscr<-hscr_gene_counts %>% gather(colnames(hscr_gene_counts)[2:29], key = "id", value = "count")
gather_hscr<-inner_join(coldat,gather_hscr)

gather_hscr$myorder<-factor(gather_hscr$genotype, levels = c('WT','Piebald_Het','Piebald_Hom','Ret_Het_Piebald_Het','Ret_Het_Piebald_Hom'),ordered=TRUE)

ggplot(gather_hscr) +  ggtitle("E14.5 Expression") + geom_boxplot(aes(x = myorder, y = count))+facet_wrap(~gene_names, scales="free_y") +theme(axis.text.x =element_text(size=0),axis.title.y=element_text(size=18), plot.title = element_text(size=20))
#+ylim(0,3000)

#with error bars, note this does not do multiple test correction
library(ggpubr)
ggerrorplot(gather_hscr, x = "myorder", y = "count",facet.by = "gene_names",add = "boxplot",scales="free",desc_stat="median")+ theme(axis.text.x =element_text(size=5,angle=25,hjust=1),axis.title.y=element_text(size=14), plot.title = element_text(size=14)) + stat_compare_means(aes(label = ..p.signif..), method = "t.test", ref.group = "WT",label.x.npc = 0.5,label.y.npc = .9)

#cleaner version of above error bar plot
ggerrorplot(gather_hscr, x = "myorder", y = "count",facet.by = "merged_col",add = "boxplot",ylim = c(0, 5000),desc_stat = "median",ggtheme = theme_gray()) + stat_compare_means(aes(label = ..p.signif..), ref.group = "Ret_WT",label.x.npc = 0.5,label.y.npc = .95)+theme(strip.background = element_blank(),strip.text.x = element_blank())

#can use for multiple test correction, manually edit pvalue 
library(rstatix)
stat.test <- gather_hscr %>% group_by(gene_names)%>%t_test(count~genotype) %>%adjust_pvalue(method = "hochberg") %>%add_significance("p.adj")
stat.test.filter<-subset(stat.test, p.adj.signif != "ns")
ggerrorplot(gather_hscr, x = "myorder", y = "count",facet.by = "gene_names",add = "boxplot",scales="free",desc_stat="median")+ theme(axis.text.x =element_text(size=5,angle=25,hjust=1),axis.title.y=element_text(size=14), plot.title = element_text(size=14)) + stat_pvalue_manual(stat.test.filter, y.position = 700, step.increase = 0.1,label = "p.adj.signif")

