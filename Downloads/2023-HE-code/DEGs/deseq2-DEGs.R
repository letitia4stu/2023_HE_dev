rm(list=ls())
library(GenomicFeatures)

library(DESeq2)

library(dplyr)
setwd("/Users/lichaoran/dataprocess/2016-nature-HE/HE_bowtie_new_version/")

###load and set output file

file_in <- "deseq_count.csv"

file_design <- "deseq_7group_design.csv"

file_compare <- "he201655_sample_compare.csv"

file_deg_num = paste("data/","DE_",file_in,sep="")  ##number of DEGs in all compare

file_final_csv  = paste("data/","DE",file_in, "Final_Out.csv",sep="")

file_final_genelist = paste("data/","DEG_geneid_allcomapre.txt",sep="")

################################################################

###########get DEGs in different compare counts file

data_in = read.csv(file_in, head=TRUE,row.names =1, check.names = FALSE)

mycompare=read.csv(file_compare,head=TRUE)

mydesign=read.csv(file_design, head = TRUE)[,2:4]

##filter counts

countData = as.data.frame(data_in)

dim(countData)

mycounts_filter <- countData[rowSums(countData) != 0,]

dim(mycounts_filter)

##check group

head(mycompare)

head(mydesign)

total_num=dim(mycompare)[1]  ##get the num of groups

tracking = 0

gene_num_out = c()

pvalue_cut = 0.01

condition_name = c()

for (index_num in c(1:total_num)){
  
  tracking = tracking + 1
  
  test_name = as.character(mycompare[index_num,1])
  
  group1 = as.character(mycompare[index_num, 2])
  
  group2 = as.character(mycompare[index_num, 3])
  
  sample1 = mydesign[mydesign$Group == group1,]
  
  sample2 = mydesign[mydesign$Group == group2,]
  
  allsample = rbind(sample1, sample2)
  
  counts_sample = as.character(allsample$counts_id)
  
  groupreal = as.character(allsample$Group)
  
  countData = mycounts_filter[, counts_sample]
  
  colData = as.data.frame(cbind(counts_sample, groupreal),stringsAsFactors = TRUE)
  
  names(colData) = c("sample", "condition")
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                
                                colData = colData,
                                
                                design = ~ condition)
  
  dds <- DESeq(dds)
  
  resSFtreatment <- results(dds, cooksCutoff =FALSE, contrast=c("condition",group2,group1))
  
  out_test = as.data.frame(resSFtreatment)
  
  final_each = cbind( out_test$log2FoldChange,out_test$pvalue,out_test$padj)
  
  rownames(final_each) = resSFtreatment@rownames  
  
  names = c('(logFC)','(pvalue)','(Qvalue)')
  
  final_name = paste(test_name,'_',names,sep="")
  
  colnames(final_each) = final_name
  
  if (tracking == 1){
    
    final_table = final_each
    
  }else{
    
    final_table = cbind(final_table, final_each)
    
  }
  
  gene_sel = out_test[((!is.na(out_test$pvalue))&(!is.na(out_test$log2FoldChange)))&out_test$pvalue < 0.05 & abs(out_test$log2FoldChange) >= 1 & out_test$padj<0.05, ]
  
  gene_sel<-na.omit(gene_sel)  ##delete rows contain NA
  
  gene_sel_up = gene_sel [gene_sel$log2FoldChange>0,]
  
  gene_sel_do = gene_sel [gene_sel$log2FoldChange<0,]
  
  file_out_up = paste("data/","UP_",test_name,".txt",sep="")
  
  file_out_do = paste("data/","DOWN_",test_name,".txt",sep="")
  
  gene_list_up = rownames(gene_sel_up)
  
  gene_list_do = rownames(gene_sel_do)
  
  all_ub_down = list(gene_list_up, gene_list_do)
  
  nameup <- paste(test_name,"_up",sep="")
  
  namedown <- paste(test_name,"_down",sep="")
  
  names(all_ub_down) = c(nameup, namedown)
  
  if (tracking == 1){
    
    final_genelist = all_ub_down
    
  }else{
    
    final_genelist = c(final_genelist, all_ub_down)
    
  }
  
  gene_num_up = length(gene_list_up)
  
  gene_num_do = length(gene_list_do)
  
  write.table(gene_sel_up, file= file_out_up , row.names = TRUE,col.names = TRUE)
  
  write.table(gene_sel_do, file= file_out_do , row.names = TRUE,col.names = TRUE)
  
  condition_name = c(condition_name, paste("UP",test_name,sep=""),paste("DO",test_name,sep=""))
  
  gene_num_out = c(gene_num_out,gene_num_up,gene_num_do)
  
}

final_DEGs_list<-do.call(cbind, lapply(lapply(final_genelist, unlist),`length<-`, max(lengths(final_genelist))))

final_DEGs_list

out_final2 = cbind(condition_name, gene_num_out)

colnames(out_final2) = c("Tests", "DEG number")

write.table(final_DEGs_list, file=file_final_genelist, row.names = FALSE, sep="\t", na = "") ###all DEGs list

write.table(out_final2, file=file_deg_num, row.names = FALSE, sep="\t")    ##DEG number in different compare

