using MultivariateStats
# suppose Xtr and Xte are training and testing data matrix,# with each observation in a column
# train a PCA modelM = fit(PCA, Xtr; maxoutdim=100)
# apply PCA model to testing setYte = transform(M, Xte)
# reconstruct testing observations (approximately)Xr = reconstruct(M, Yte)
