library('pathview')
library('tidyr')
rm(list=ls())
#文件路径
setwd("/Users/lichaoran/fsdownload/win10-exp/ortho-dev-all")

# 获取当前工作目录（即当前文件夹）中的所有csv文件
folder_path <- "/Users/lichaoran/fsdownload/win10-exp/ortho-dev-all"
csv_files <- list.files(path = folder_path, pattern = "*.csv")


# 循环遍历每个csv文件，并将其读取为矩阵
for (csv_file in csv_files) {
  
  # 读取csv文件为一个矩阵
  matrix_data <- as.matrix(read.csv(file.path(folder_path,csv_file), row.names = 1))
  
  # 获取csv文件名的前缀
  prefix <- tools::file_path_sans_ext(csv_file)
  
  # 新建一个同名的矩阵，并将读取出的数据存入其中
  assign(paste0(prefix, "_matrix"), matrix_data)
  
  # 打印新建的矩阵的名称和维度
  cat(sprintf("Matrix '%s' is created with dimension %d x %d.\n", paste0(prefix, "_matrix"), nrow(matrix_data), ncol(matrix_data)))
  
}
# 获取所有已经生成的矩阵名
matrix_names <- ls(pattern = "_matrix")
# 删除指定的矩阵名
#matrix_names <- matrix_names[!matrix_names %in% c("Elongation_2_matrix", "Segmentation_3_matrix")]

# 循环遍历每个矩阵，计算出新的数据框
#for (matrix_name in matrix_names) {
  
  # 从矩阵名中提取前缀
  #prefix <- gsub("_matrix", "", matrix_name)
  
  # 获取指定的矩阵，并进行运算
  #gene.ensprot <- get(matrix_name)
  #gexp <- log2(gene.ensprot + 1)
  #gene.data.scale <- as.data.frame(t(apply(gene.ensprot,1,function(x){scale(x)[,1]})))
  
  # 以矩阵名作为前缀，生成新的数据框
 # assign(paste0(prefix, "_matrix"), gene.data.scale)
}

# 定义多个示例 gene_data_scale
# 获取满足pattern条件的对象名称
matrix_names <- ls(pattern = "_matrix")


# 根据对象名称获取对象的值，并将它们放入list中
gene_data_scale_list <- list()
for (name in matrix_names) {
  gene_matrix <- get(name)
  gene_data_scale_list <- c(gene_data_scale_list, list(gene_matrix))
}
'''
ls(pattern = "_matrix")
gene_data_scale_list <- list( Cleavage_1_matrix,
                              Elongation_2_matrix,
                              Segmentation_3_matrix,
                              Limb_bud_and_pharynx_4_matrix,
                              Ganglia_and_midgut_5_matrix,
                              Morphologically_distinction_6_matrix,
                              Hatching_7_matrix)
'''

# 创建文件夹，用于存储多个 gene_data_scale 的结果文件夹
#dir.create("pathway_results")
# 循环创建文件夹
for (i in 8:10) {
 dir.create(paste0(i))
}
# 调用函数进行循环操作
# 定义需要循环的pathwayid向量
params <- read.table('/Users/lichaoran/Downloads/dme-kegg/dme-all-kegg/dme-all-de00514-00511.txt',sep = '\t')
#中间报错,接上
#params <- params[110:nrow(params), ]
#changge dme009880 to  009880
params$V1 <- substring(params$V1, 4)
#
#paras <- separate(paras,col='V2',into=c('V2','V3'),sep='- Drosophila')
#paras$V2 <- gsub(" ", "", paras$V2)
params$V2 <- gsub("/", "", params$V2)
colnames(params)=c('pathwayids','suffix')
##

#seq_along(gene_data_scale_list
pathview_loop <- function(gene_data_scale_list, params) {
  for (j in 8:10) {
    gene_data_scale <- gene_data_scale_list[[j]]
    df_list <- list()
    
    for (i in 1:nrow(params)) {
      
      current_path <- params$pathwayids[i]
      current_suffix <- params$suffix[i]
      
      pv.out <- pathview(gene.data = gene_data_scale,
                         pathway.id = current_path, species = "dme",
                         limit = list(gene = 1.5), 
                         low = 'blue',
                         mid = 'white',
                         high = 'red',
                         kegg.native = T, gene.idtype = 'UNIPROT',
                         out.suffix = current_suffix, same.layer = F)
      
      pv <- data.frame(lapply(pv.out[["plot.data.gene"]], as.character), stringsAsFactors=FALSE)
      
      df_list[[current_path]] <- pv
      
      filename <- paste0(j,current_suffix, "_", current_path, "-kegg.csv")
      filepath <- file.path(j, filename)
      write.csv(pv, filepath, row.names = FALSE)
      
    }
    
    N <- length(df_list)
    
    ratio_df <- data.frame(ratio = rep(NA, N), pid = params$pathwayids)
    
    for (i in 1:N) {
      df <- df_list[[i]]
      allmap <- sum(df$all.mapped!= "")
      keggm  <- sum(df$kegg.names!= "")
      ratio_df[i, "ratio"] <- allmap / keggm 
    }
    
    ratio_df2 <- merge(ratio_df, params, by.x = 'pid', by.y = 'pathwayids')
    
    filename <- paste0("dme-kegg-cover-ratio_", j, ".csv")
    filepath <- file.path('pathway_results',filename)
    write.csv(ratio_df2, filepath, row.names = FALSE)
    
  }
  
}

pathview_loop(gene_data_scale_list, params)

# 设置工作目录
setwd("/Users/lichaoran/fsdownload/win10-exp/ortho-dev-all/pathway_results") 
rm(list=ls())
# 读取文件名，存储到文件名列表中
file_list <- list.files(pattern="*.csv")

# 创建一个空 DataFrame，用于存储所有 CSV 的数据
combined_data <- data.frame()

# 循环遍历文件名列表
for (f in file_list) {
  # 读取每个 CSV 文件
  data <- read.csv(f, header=TRUE, sep=",", stringsAsFactors=FALSE)
  
  # 删除第一列
  data <- data[-1]
  #data <- data[-7,]
  
  # 获取第二列名称
  col_name <- data[,2]
  
  # 将第二列数据存储为列名为第三列内容的 DataFrame
  new_data <- t(data.frame(data[,1]))
  colnames(new_data) <- col_name
  
  # 合并所有数据
  combined_data <- rbind(combined_data, new_data)
}

# 显示结果
write.csv(t(combined_data),'hebing_dev_stages_kegg_3.csv')
