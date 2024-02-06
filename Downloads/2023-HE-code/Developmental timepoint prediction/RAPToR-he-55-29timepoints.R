rm(list=ls())
##read 55 sample read count expression txt
data_folder <- "/Users/lichaoran/dataprocess/2016-nature-HE/HE_bowtie_new_version/"
g_file_he201655 <- paste0(data_folder, "he201655.txt")
X_he201655 <- read.table(g_file_he201655, h=T,sep = '\t', stringsAsFactors = F, row.names = 1)
# To TPM
x=colSums(X_he201655)
countsum=as.data.frame(t(rbind(x,X_he201655)))
sct=(countsum/countsum[,1])*1000000
X_he20162=sct[,-1]
X_he201655=as.data.frame(t(X_he20162))
tail(X_he201655)

# P_he201655
P_he201655=read.csv("/Users/lichaoran/dataprocess/2016-nature-HE/HE_bowtie_new_version/P_he201655.csv")

colnames(P_he201655)[3] <- "strain"
P_he201655$strain <- factor(P_he201655$strain)
P_he201655$title <- make.names(P_he201655$title)
head(P_he201655, n = 5)

# get age 
P_he201655$age <- as.numeric(P_he201655$age)
P_he201655$age2 <- as.numeric(P_he201655$age2)
X_he201655 <- X_he201655[,P_he201655$title]
he201655 <- list(g = X_he201655, p = P_he201655)

data_folder2 <- "/Users/lichaoran/Downloads/test/extdata/"
save(he201655, file = paste0(data_folder2, "he201655.RData"), compress = "xz")


load("/Users/lichaoran/Downloads/test/extdata/he201655.RData")

tail(he201655$g)
he201655$g <- limma::normalizeBetweenArrays(he201655$g, method = "quantile")
he201655$g <- log1p(he201655$g)

he201655$g[1:5,1:4]
head(he201655$p, n = 5)
pca_he201655 <- stats::prcomp(t(he201655$g), rank = 25,
                              center = TRUE, scale = FALSE)
par(mfrow = c(2,4))
invisible(
  sapply(seq_len(8), function(i){
    plot(he201655$p$age, pca_he201655$x[,i], lwd = 2,col = he201655$p$strain,
         xlab = "age", ylab = "PC", main = paste0("PC", i))
    
    # connect the dots
    sapply(seq_along(levels(he201655$p$strain)), function(l){
      s <- which(he201655$p$strain == levels(he201655$p$strain)[l])
      points(he201655$p$age[s], pca_he201655$x[s,i], col = l, 
             type = 'l', lty = 2)
    })
    
    if(i == 1)
      legend("topleft", bty = 'n', legend = c("Hypsibius dujardini"),
             pch = c(rep(1)), lty = c(rep(NA, 1)), col = c(1), lwd = 3)
  })
)

######
#Correlation
cor_he201655 <- cor(he201655$g, method = "spearman")
#Plots
# Heatmap
ord <- order(he201655$p$age)
heatmap(cor_he201655[ord, ord], Colv = NA, Rowv = NA, scale = "none", keep.dendro = F, margins = c(1,1)
        ,labRow = "", labCol = "")
#heatmap(cor_he201655[ord, ord], Colv = NA, Rowv = NA, scale = "none", keep.dendro = F, margins = c(1,1),
#RowSideColors = as.numeric(dsaeschimann2017$p$organism_ch1[ord]), labRow = "", labCol = "")
par(xpd = T) # text may have to be tweaked to plot size
mtext(text = unique(he201655$p$age), side = 1, line = 4, at = seq(-.1,1.05, l = 11))

#Validation
#Predict
# setup newdat
n.inter <- 55 # nb of new timepoints
newdat <- data.frame(
  age = seq(min(he201655$p$age), max(he201655$p$age), l = n.inter),
  strain = rep("Hypsibius exemplaris", n.inter) # we want to predict as N2 
)
head(newdat)
#age              strain
#1  30.0000 Hypsibius dujardini
#2 126.8852 Hypsibius dujardini
#3 223.7705 Hypsibius dujardini
#4 320.6557 Hypsibius dujardini
#5 417.5410 Hypsibius dujardini
#6 514.4262 Hypsibius dujardini

# predict 
pred_m_he201655 <- predict(m_he201655, newdata = newdat)
pred_m_he201655_comp <- predict(m_he201655, newdata = newdat, as.c = TRUE)

#Build reference & stage samples
# make a 'reference object' 
r_he201655 <- list(interpGE = pred_m_he201655, time.series = newdat$age)
ae_he201655 <- ae(he201655$g, r_he201655$interpGE, r_he201655$time.series)
r_he201655 <- list(interpGE = pred_m_he201655, time.series = newdat$age)
ae_test_he201655 <- ae(he201655$g, r_he201655$interpGE, r_he201655$time.series)
par(mfrow = c(1,2))
rg <- range(c(ae_test_he201655$age.estimates[,1], he201655$p$age))
# Plot 1
plot(ae_test_he201655$age.estimates[,1]~he201655$p$age, 
     xlab = "Chronological age", ylab = "Estimated age (he201655)", 
     xlim = rg, ylim = rg,
     main = "Chron. vs Estimated ages for he201655\n(on he201655 reference)", lwd = 2, 
     col = factor(he201655$p$strain))
# connect the dots
invisible(sapply(levels(factor(he201655$p$strain)), function(l){
  s <- he201655$p$strain == l
  points(ae_test_he201655$age.estimates[s,1]~he201655$p$age[s], type = 'l', 
         lty = 2, col = which(l==levels(factor(he201655$p$strain))))
}))
abline(a = 0, b = 1, lty = 3, lwd = 2) # x = y

legend("bottomright", legend = c("Hypsibius exemplaris"), 
       lwd=3, col=c(1:4, 1), bty='n', pch = c(1,1,1,1,NA), lty = c(rep(NA, 4), 3))


plot(ae_he201655, groups = he201655$p$strain, show.boot_estimates = T)
write.csv(ae_test_he201655[["age.estimates"]],"he201655-age.estimates.csv")



