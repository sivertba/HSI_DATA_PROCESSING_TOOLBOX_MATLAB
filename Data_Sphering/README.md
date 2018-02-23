# Data Sphering

The purpose of data sphering is to allow users to analyze data structure characterized by high-order
statistics. Before doing so, the data samples characterized by the first two orders of statistics, that is,
mean and variances/covariances must be removed. The data sphering is designed to accomplish this
task. It first removes the data sample mean by setting data set centered at the origin and then
de-correlates data samples by zeroing all covariances via diagonalization of data sample covariance
matrix. Finally, it normalizes data sample variances to 1 by placing all de-correlated data samples on
the unit sphere. 

So, by means of matrix diagonalization and variance normalization, all data samples
characterized by high-order statistics are either inside the sphere characterized by sub-Gaussian or
outside the sphere characterized by super-Gaussian. Technically speaking, a whitening processing is
a part of data sphering that only de-correlates data samples without normalization. However, in statistical
signal processing and communications community as well as in many application whitening is
indeed data sphering. 

In the book "Hyperspectral Data Processing: Algorithm Design and Analysis", we particularly make a distinction between them. In other words,
whitening only de-correlates data samples by making co-variances zero but does not normalize variances
to 1. It is only a part of data sphering. Details of data sphering can be found in Chapter 6 of the "Hyperspectral Data Processing: Algorithm Design and Analysis" book.