###############################################
############### Counts to RPKM ################
###############################################
# STEPS: 
# 1. Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
# 2. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
# 3. Divide the RPM values by the length of the gene, in kilo bases. This gives you RPKM.

##############
require (dplyr)
rm(list = ls())
#read count file
mat<-read.csv("Counts.csv")
row.names(mat)<-mat[,1]
gene<-as.data.frame(mat[,3]) ## read gene length
colnames(gene)<-"gene"

#calculation begins
total<- mat %>% summarise_all(sum)
per_million<-total[1,]/1000000

#divide each read count by per million factor
seq<-colnames(mat)
for (i in seq){
    print(i)
        mat[,i]<-(mat[,i])/(per_million[,i])
}

#divide by gene length
seq1<-c(1:nrows)   #nrows= total number of genes
for (i in seq1){
    print(i)
        mat[i,]<-(mat[i,])/(gene[i,])
}
final<-(mat)*1000
write.csv(final,file="rpkm_final.csv")
