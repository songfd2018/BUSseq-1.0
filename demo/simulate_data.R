rm(list=ls())
library(ggplot2)
library(Rtsne)

set.seed(12356)
################################
# Set the Synthetic Parameters #
################################
#The number of batches
B<-3

#The number of cells per batch
nb<-c(150,150,150)

#The total number of cells
N<-sum(nb)

#The number of genes
G<-1000

#The number of cell types
K<-4

# the project name
proj <- "demo"

#The first column of gamma.syn denotes the intercept of 
#the logistic regression for dropout events
#The second column of gamma.syn denotes the odds ratios 
#of the logistic regression for dropout events
gamma.syn<-matrix(0,B,2)
gamma.syn[1,]<-c(-1,-0.05)
gamma.syn[2,]<-c(-1,-0.05)
gamma.syn[3,]<-c(-0.5,-0.1)

#the log-scale baseline expression levels
alpha.syn<-rep(1,G)
alpha.syn[1:(G*0.1)]<-4

#the cell-type effects 
beta.syn<-matrix(0,G,K)

#the first cell type is regarded as the reference cell type
beta.syn[,1] <- 0

#the effects of the second cell type
beta.syn[1:(G * 0.1),2] <- -3
beta.syn[(G * 0.1 + 1):(G * 0.3),2] <- 3
beta.syn[(G * 0.3 + 1):G,2] <- 0

#the effects of the third cell type
beta.syn[1:(G * 0.1),3]<- -3
beta.syn[(G * 0.1 + 1):(G * 0.3),3] <- 0
beta.syn[(G * 0.3 + 1):(G * 0.45),3] <- 3
beta.syn[(G * 0.45 + 1):G,3] <- 0

#the effects of the forth cell type
beta.syn[1:(G * 0.1),4]<- -3
beta.syn[(G * 0.1 + 1):(G * 0.45),4] <- 0
beta.syn[(G * 0.45 + 1):(G * 0.6),4] <- 3
beta.syn[(G * 0.6 + 1):G,4] <- 0

#the batch effects
nu.syn<-matrix(NA,G,B)

#the first batch is regarded as the reference batch
nu.syn[,1] <- 0

#the effect of the second batch
nu.syn[1:(G * 0.4),2]<-2
nu.syn[(G * 0.4 + 1):(G * 0.8),2]<-1
nu.syn[(G * 0.8 + 1):G,2]<-2

#the effect of the third batch
nu.syn[1:(G * 0.3),3]<-3
nu.syn[(G * 0.3 + 1):(G * 0.6),3]<-2
nu.syn[(G * 0.6 + 1):G,3]<-1

#the cell-specific size factors
delta.syn <- list()
for(b in 1:B){
  delta.syn[[b]] <- rep(NA, nb[b])
}

#The frist cell in each bathc is regarded as the reference cell
delta.syn[[1]][1:(nb[1] * 0.5)]<-0
delta.syn[[1]][(nb[1] * 0.5 + 1):(nb[1] * 0.9)]<-1
delta.syn[[1]][(nb[1] * 0.9 + 1):nb[1]]<-2

#the second batch
delta.syn[[2]][1:(nb[2] * 0.5)]<-0
delta.syn[[2]][(nb[2] * 0.5 + 1):(nb[2] * 0.7)]<-2
delta.syn[[2]][(nb[2] * 0.7 + 1):nb[2]]<--1

#the third batch
delta.syn[[3]][1:(nb[3] * 0.3)]<-0
delta.syn[[3]][(nb[3] * 0.3 + 1):nb[3]]<-1

#the batch-specific and gene-specific overdispersion parameters
phi.syn<-matrix(10,B,G)

#the cell-type proportions in each batch
pi.syn <- matrix(NA,K,B)

#the frist batch
pi.syn[,1]<-c(0.4,0.3,0.2,0.1)

#the second batch
pi.syn[,2]<-c(0.2,0.24,0.2,0.36)

#the third batch
pi.syn[,3]<-c(0.2,0.18,0.32,0.3)

##############################################
# Simulate Latent Varibles and Observed data #
##############################################
#the cell-type indicators of each cell
w<-list()

#the first batch
nw<-nb[1] * pi.syn[,1]
w[[1]]<- rep(1:4,nw)

#the second batch
nw<-nb[2] * pi.syn[,2]
w[[2]]<- rep(1:4,nw)


#the third batch
nw<-nb[3] * pi.syn[,3]
w[[3]]<- rep(1:4,nw)

#the indicators for dropout events
z<-list()

#the underlying true expression levels
x<-list()

#the observed expression levels
y<-list()

#the logarithm of mean expreesion level of each gene in each cell
log.mu<-list()

for(b in 1:B){
  z[[b]] <- matrix(NA, G, nb[b])
  x[[b]] <- matrix(NA, G, nb[b])
  y[[b]] <- matrix(NA, G, nb[b])
  log.mu[[b]] <- matrix(NA, G, nb[b])
}

#generate the latent variable and observed data
for(b in 1:B){
  for(i in 1:nb[b]){
    log.mu[[b]][,i] <- alpha.syn + beta.syn[,w[[b]][i]] 
    log.mu[[b]][,i] <- log.mu[[b]][,i] + nu.syn[,b]
    log.mu[[b]][,i] <- log.mu[[b]][,i] + delta.syn[[b]][i]
    
    for(j in 1:G){
      x[[b]][j,i]<-rnbinom(1,phi.syn[b,j],
                           mu=exp(log.mu[[b]][j,i]))
      logit_pi <- gamma.syn[b,1] + gamma.syn[b,2] * x[[b]][j,i]
      
      z[[b]][j,i]<-rbinom(1,1,prob = 1/(1+exp(-logit_pi)))
      if(z[[b]][j,i]==0){
        y[[b]][j,i]<- x[[b]][j,i]
      }else{
        y[[b]][j,i]<- 0
      }
    }
  }
}

if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

save.image(paste0(proj,"_v1.RData"))

# Write out the read count matrix, dimsion information and metadata
readcount <- do.call(cbind,y)
write.table(readcount, file = paste0("RawCountData/count_data_",proj,"_v1.txt"),row.names = FALSE,col.names = FALSE)

dim <- c(N,G,B,nb)
write.table(dim, file = paste0("RawCountData/dim_",proj,"_v1.txt"),row.names = FALSE, col.names = FALSE)

metadata <- data.frame( batch = paste0("Batch_",rep(1:B,nb)),
                        celltype = paste0("Type_",unlist(w)))
write.table(metadata, file = paste0("RawCountData/metadata_",proj,"_v1.txt"), quote = FALSE, row.names = FALSE)

##############################################
# Plot PCA and t-SNE plots before correction #
##############################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#set cell type colorings
celltype_color<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
celltype.cols <- celltype_color[unlist(w)]
plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=celltype.cols[subset])
  dev.off()
}

batch_color<-c("#F88A7E", "#FFD87D", "#ABD978")
batch.cols<-batch_color[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=batch.cols[subset]) 
  dev.off()
}

# t-SNE
set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(readcount))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

unc_by_celltype<- "Image/tsne_simulation_uncorrected_by_celltype.jpeg"
plot_by_celltype(unc_by_celltype, tsne_uncorrected$Y)

unc_by_batch <- "Image/tsne_simulation_uncorrected_by_batch.jpeg"
plot_by_batch(unc_by_batch, tsne_uncorrected$Y)

# legend
pdf(file="Image/legend.pdf", width=10, height=8)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=names(table(metadata$batch)), pch=20, cex=2.5, col=batch_color, title = "Batch", bty="n")
legend(x=-0, y=4, legend=names(table(metadata$celltype)), pch=20, cex=2.5, col=celltype_color, title = "Cell Type", bty="n")
dev.off()