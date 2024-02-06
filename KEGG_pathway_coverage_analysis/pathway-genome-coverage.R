####划分同源的重做的
setwd('/Users/lichaoran/fsdownload/tmdchongzuo/')
library('pathview')
library('clusterProfiler')
library('org.Dm.eg.db')
#gene.ensprot=read.csv('/Users/lichaoran/fsdownload/dme-he-id-2313.csv',header = T,row.names = 1)
#gene.ensprot=read.csv('/Users/lichaoran/fsdownload/win10-exp/dme-he-id-5883.csv',header = T,row.names = 1)
gene.ensprot=read.csv('/Users/lichaoran/fsdownload/win10-exp/dme-he-id-haoduo.csv',header = T,row.names = 1)
gene.ensprot=as.matrix(gene.ensprot)
#单条通路尝试
pv.out <- pathview(gene.data = gene.ensprot[,1], #基因数据
                   gene.idtype = "Uniprot",  #基因使用的I
                   pathway.id = '04330',#通路ID
                   species = "dme", 
                   same.layer = T,
                   out.suffix = 'Notch-he-uniprot') 
# 定义需要循环的pathwayid向量
paras <- read.table('/Users/lichaoran/Downloads/dme-kegg/dme-all-kegg/dme-all-de00514-00511.txt',sep = '\t')
#changge dme009880 to  009880
paras$V1 <- substring(paras$V1, 4)
#
paras <- separate(paras,col='V2',into=c('V2','V3'),sep='- Drosophila')
paras$V2 <- gsub(" ", "", paras$V2)
paras$V2 <- gsub("/", "", paras$V2)
##

df_list=list()
# for循环对不同的pathwayid进行循环操作，并写入不同的文件
for(i in 1:nrow(paras)) {
  # 从数据框params中获取当前循环的pathwayid和suffix参数
  current_path <- paras$V1[i]
  current_suffix <- paras$V2[i]
  # 进行pathview操作，并指定pathwayid参数为当前循环的i值
  pv.out <- pathview(gene.data = gene.ensprot,
                     pathway.id = current_path, species = "dme",
                     limit = list(gene = 1.5), 
                     low = 'blue',
                     mid = 'white',
                     high = 'red',
                     kegg.native = T, gene.idtype = 'UNIPROT',
                     out.suffix = current_suffix, same.layer = F)
  # 提取plot.data.gene数据并将其转为数据框
  pv <- data.frame(lapply(pv.out[["plot.data.gene"]], as.character), stringsAsFactors=FALSE)
  # 将当前pathwayid的数据写入文件
  df_list[[current_path]]=pv
  write.csv(pv, paste0("/Users/lichaoran/fsdownload/tmdchongzuo/pathway-concrete",current_path,current_suffix, "-kegg.csv"))
}

'''
df_list <- list()

for (i in 1:nrow(paras)) {
  current_path <- paras$V1[i]
  current_suffix <- paras$V2[i]
  
  # 使用tryCatch包装pathview操作
  tryCatch({
    pv.out <- pathview(gene.data = gene.ensprot,
                       pathway.id = current_path, species = "dme",
                       limit = list(gene = 1.5), 
                       low = 'blue',
                       mid = 'white',
                       high = 'red',
                       kegg.native = T, gene.idtype = 'UNIPROT',
                       out.suffix = current_suffix, same.layer = F)
    
    pv <- data.frame(lapply(pv.out[["plot.data.gene"]], as.character), stringsAsFactors = FALSE)
    
    # 将当前pathwayid的数据写入df_list
    df_list[[current_path]] <- pv
  }, error = function(e) {
    # 异常处理，如果无法提取plot.data.gene则跳过这个pathwayid
    cat(sprintf("Warning: Number of mappable nodes is below", current_path, e$message))
    next
  })
  
}
# 处理完所有pathwayid后，将df_list中的数据写入文件
for (i in seq_along(df_list)) {
  if (is.null(df_list[[i]])) {
    cat(sprintf("pathway %s skipped!\n", i))
    next
  }
  
  current_path <- paras$V1[i]
  current_suffix <- paras$V2[i]
  write.csv(df_list[[i]], sprintf("pathway_%s_%s.csv", current_path, current_suffix))
}
'''

#####计算节点覆盖度
# 假设我们有N个dataframe，存储在一个列表df_list中
N <- length(df_list)

# 创建一个空dataframe，用于存储比例数值
ratio_df <- data.frame(ratio = rep(NA, N),pid=names(df_list))

# 循环计算比例数值，并将结果写入ratio_df中
for (i in 1:N) {
  df <- df_list[[i]]
  allmap <- sum(df[["all.mapped"]]!= "") # 计算第三列不为空的数目
  keggm  <- sum(df[["kegg.names"]]!= "") # 计算第一列不为空的数目
  ratio_df[i, "ratio"] <- allmap / keggm # 计算比例数值并写入ratio_df中
}

# 查看计算结果
ratio_df
ratio_df2=merge(ratio_df,paras,by.x='pid',by.y='V1')

write.csv(ratio_df2,'nnngene-all-dme-kegg-cover-ratio.csv')

##
gene.ensprot2=as.matrix(read.csv('hebing-he2016-kegg-uniprot-dmeku-idtransfer-dduplicate.csv',row.names = 1))

df <- df_list[[10]]
sum(df[, 3]!= "")

######
gene.ensprot=read.csv('/Users/lichaoran/Downloads/figure2/pathplot/slide-dme-he-exp-10win.csv',row.names = 1)
gene.ensprot=as.matrix(gene.ensprot)
#合并多样本表达量数据
gene.data.scale <- as.data.frame(
  t(apply(gene.ensprot,1,function(x){scale(x)[,1]
  })))
#04310WNT
pv.out <-pathview(gene.data = gene.data.scale,
                  pathway.id = "04391", species = "dme",
                  limit = list(gene = 1.5), 
                  low = 'blue',
                  mid = 'white',
                  high = 'red',
                  kegg.native = T,gene.idtype = 'Uniprot',
                  out.suffix = "2", same.layer = F)
pv = data.frame(lapply(pv.out[["plot.data.gene"]], as.character), stringsAsFactors=FALSE)
# write file
write.csv(pv,'04310WNT-kegg.csv')