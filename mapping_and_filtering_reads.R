#Running Igor Dolgalev's Salmon pipeline on big purple

#download SNS'Seq N Slide' github code 
git clone --depth 1 https://github.com/igordot/sns 

#cd
cd ./sns/
  
  #specify fastq file location
  srun --time=2:00:00 gather-fastqs <fastq_dir>
  
  #specify reference genome
  srun --time=2:00:00 generate-settings <genome>
  
  #note <genome> can be mm10, hg38 etc. on big purple and this will 
  #default to Igor's reference genome library
  #may break if he stops maintaining code
  
  #run salmon
  srun --time=2:00:00 run rna-salmon

#converting ensembl stable transcript to gene names/symbols
#remove decimal
cut -f1 quant.salmon.tximport.counts.txt | sed -e 's/\.[0-9]*//g' -e 's/ *$//'>temp.txt
#paste new ids
cut -f2-26 quant.salmon.tximport.counts.txt >temp2.txt
paste temp.txt temp2.txt >temp3.txt

#acquired gene symbols for Ensembl stable transcript ids from biomart on ensembl.org

#note this made me realize there are ~30k lncRNAs and miRNAs, pseudogenes, mtRNA, Riken defined RNA and removal of lncRNA, miRNA, pseudogenes reduced list to ~20k genes instead of 50k stable transcripts#

#my reasoning is that the library construction probably requires special considerations such as size to actually examine these types of RNAs

#GTF file of mouse genome:https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz)
gtf <- rtracklayer::import('~/Desktop/gencode.vM28.annotation.gtf') 
gtf_df=as.data.frame(gtf)
gtf_df2 = subset(gtf_df, type == 'gene')
gtf_df3 = subset(gtf_df2, gene_type == "protein_coding")
write.csv(gtf_df3,"~/Desktop/protein_coding.csv", row.names = T)
cts_filtered<-cts[gtf_df3$gene_id,]
#export to excl and use filter to remove NA's
#leaves approx 20k genes
write.csv(cts_filtered,"~/Desktop/temp.csv")