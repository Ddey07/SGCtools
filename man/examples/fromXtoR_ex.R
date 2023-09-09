require(SGCTools)
# generate data and estimate latent correlation matrix
n <- 100
p.seq <- 5
k <- 1
p <- p.seq[k]
# Generate random correlation matrix
S <- clusterGeneration::rcorrmatrix(p)
# Generate latent normal
L <- mvtnorm::rmvnorm(n,sigma = S)
# Generate latent non-paranormal (monotone transformations)
Z <- cbind(L[,1]^3, exp(L[,2]), L[,3], L[,4], L[,5])
# Generate observed X : C, B, O, O, T
X <- cbind(Z[,1],Z[,2] >= 1, as.numeric(cut(Z[,3],breaks=c(-Inf,-0.2,0.1,0.3,Inf)))-1 ,as.numeric(cut(Z[,4],breaks=c(-Inf,-0.3,0.05,0.2,Inf)))-1,ifelse(Z[,5]>0,Z[,5],0))
# Calculate Kendall's Tau of observed data
K = Kendall_mixed(X)
# Calculate bridged estimates of latent correlation
R=fromXtoRMixed(X,use.nearPD=TRUE)

# Calculate predictions of latent normal
Lhat <-getLatentPreds(X, lat.cov = R$hatR)
