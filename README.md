SparseLogisticPCA
=================
This is an implementation of the sparse logistic PCA algorithm from "[Sparse logistic principal components analysis for binary data](http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.aoas/1287409387)" by Lee, Huang, and Hu (2010). It uses the uniform bound for the log likelihood. The function is in the file sparse_logistic_pca.R.

I attempted to recreate the SNP data that was used as an example in the paper. The SNP data comes from release 16 of [HapMap](http://snp.cshl.org/downloads/datafreeze.html.en) data. I used the [full, non-redundant data](http://snp.cshl.org/downloads/genotypes/2005-06_16c.1_phaseI/full/non-redundant/). The list of SNPs used is in table S1 of [this paper](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0001382#pone.0001382.s008) and in the file locations.csv. The data was manipulated using the file combineData.R. The final binary data is in SNPBinaryMatrix.csv

I was not able to perfectly recreate the dataset that Lee, et al. did, but the results are similar.