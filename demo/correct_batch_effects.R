rm(list=ls())
library(ggplot2)
library(Rtsne)

proj <- "demo"
K <- 4

###########################
# Correct read count data #
###########################
# Load dimension
dim <- read.table(paste0("RawCountData/dim_",proj,"_v1.txt"))
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v1.txt"),header = TRUE)

# load w_est
w.est <- read.table(paste0("Inference_K",K,"/w_est.txt"))
w_BUSseq <- unlist(w.est)

# load alpha_est
alpha.est <- read.table(paste0("Inference_K",K,"/alpha_est.txt"))
alpha.est <- unlist(alpha.est)
# load beta_est
beta.est <- read.table(paste0("Inference_K",K,"/beta_est.txt"))
beta.est <- matrix(unlist(beta.est),G,K)
logmu.est<-beta.est+alpha.est

# load nu_est
nu.est <- read.table(paste0("Inference_K",K,"/nu_est.txt"))
nu.est <- matrix(unlist(nu.est),G,B)

# load delta_est
delta.est <- read.table(paste0("Inference_K",K,"/delta_est.txt"))
delta.est <- unlist(delta.est)

# load gamma_est
gamma.est <- read.table(paste0("Inference_K",K,"/gamma_est.txt"))
gamma.est <- matrix(unlist(gamma.est),B,2)

# load phi_est
phi.est <- read.table(paste0("Inference_K",K,"/phi_est.txt"))
phi.est <- matrix(unlist(phi.est),G,B)

# load pi_est
pi.est <- read.table(paste0("Inference_K",K,"/pi_est.txt"))
pi.est <- matrix(unlist(pi.est),B,K)
order.est<-order(pi.est[1,],decreasing = T)

# load p_est
p.est <- read.table(paste0("Inference_K",K,"/p_est.txt"))

# load tau0_est
tau0.est <- read.table(paste0("Inference_K",K,"/tau0_est.txt"))

# load PPI_est
PPI.est <- read.table(paste0("Inference_K",K,"/PPI_est.txt"))
D.est <- unlist(read.table(paste0("Inference_K",K,"/IG_est.txt")))

# load X_imputed
x_imputed <- read.table(paste0("Inference_K",K,"/x_imputed.txt"))

############################
# Batch effects correction #
############################
adjusted_values = function(Truereads, .B, .nb, .N, .G, .K,
                           .alpha,.beta,.nu,.delta,.phi,.w){
  .w = .w + 1
  .b = unlist(sapply(seq_len(.B), function(b) rep(b, .nb[b])))
  #uniform random numbers are obtained outside of parallel code to handle seed issue.
  global_u = array(runif(.N*.G),c(.G,.N))
  get_corrected_read = function(i){
    b = .b[i]
    # The following are vectors of length .G
    p_x <- pnbinom(Truereads[,i], size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
    p_xminus1 <- pnbinom(Truereads[,i] - 1, size = .phi[,b], mu = exp(.alpha + .beta[,.w[i]] + .nu[,b] + .delta[i]))
    local_u = mapply(function(u) min(u, .999), global_u[,i]*(p_x - p_xminus1) + p_xminus1)
    return(qnbinom(local_u, size = .phi[,1], mu = exp(.alpha + .beta[,.w[i]])))
  }
  return(mapply(get_corrected_read, seq_len(.N)))
  
}

print("Calculate corrected read counts:")
x_corrected<-adjusted_values(x_imputed, B, nb, N, G, K,
                             alpha.est,beta.est,nu.est,delta.est,phi.est,w_BUSseq)
write.table(x_corrected, file = "x_corrected.txt", row.names = FALSE, col.names = FALSE)

#######
# ARI #
#######
library(mclust)
ARI_BUSseq <- adjustedRandIndex(metadata$celltype,w_BUSseq)

#############################################
# Plot PCA and t-SNE plots after correction #
#############################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#set cell type colorings
celltype_color<-c("#EB4334","#FBBD06","#35AA53","#4586F3")
celltype.cols <- celltype_color[metadata$celltype]
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
  legend("bottomleft", legend=names(table(metadata$celltype)), pch=20, cex=4, col=celltype_color, title = "Cell Type")
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
  legend("bottomleft", legend=names(table(metadata$batch)), pch=20, cex=4, col=batch_color, title = "Batch")
  dev.off()
}

set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log1p(x_corrected[D.est==1,]))))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

BUSseq_by_celltype<- paste0("Image/tsne_",proj,"_BUSseq_by_celltype.jpeg")
plot_by_celltype(BUSseq_by_celltype, tsne_BUSseq_dist$Y)

BUSseq_by_batch <- paste0("Image/tsne_",proj,"_BUSseq_by_batch.jpeg")
plot_by_batch(BUSseq_by_batch, tsne_BUSseq_dist$Y)
