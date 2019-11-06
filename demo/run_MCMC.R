rm(list=ls())
######################
# Run MCMC algorithm #
######################

set.seed(2350)
# the range of K to test
k_sel <- 2:5
# project name
proj <- "demo"
# number of the overall MCMC iterations
iter_num <- 1000
# number of iterations for each output
iter_out <- 500
# number of burn-in iterations
burnin <- 500
# number of cores running in parallel
cores <- 4 

for(k in k_sel){
  seed <- round(runif(1,0,10000))
  # run MCMC alogrithm
  system(paste0("../BUSseq -d./ -r./RawCountData/ -p ",proj," -v 1 -K ",k,
                " -i ",iter_num," -o ",iter_out, " -s ",seed," -c ",cores))
  
  # conduct posterior inference
  system(paste0("../BUSseq_inference -d./ -r./RawCountData/ -p ",proj," -v 1 -K ",k,
                " -i ",iter_num," -b ",burnin," -c ",cores))
}

###############################
# Select the optimal K by BIC #
###############################
BIC_k <- rep(Inf, length(k_sel))
for(i in 1:length(k_sel)){
  BIC_k[i] <- unlist(read.table(paste0("Inference_K", k_sel[i], "/BIC.txt")))
}


jpeg(paste("Image/",proj,"_v1_BIC.jpg",sep=""), width = 800, height = 600)
par(mar=c(5,5,2,2))
plot(k_sel, BIC_k, type="n", xlab = "K", ylab = "BIC", cex.axis = 3, cex.lab = 3)
points(k_sel, BIC_k, type="b", pch=19, cex=3)
dev.off()

k_opt <- k_sel[which.min(BIC_k)]
message(paste0("The optimal cell type number is ",k_opt,"."))