function number=cchengNWHFC(HIM,t)
%% Code Comments:
% NWHFC gives the VD number estimated by given false alarm property using
% NWHFC method.
%
% There are two parameters, NWHFC(HIM,t) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% t is the false alarm probability.
%
% HFC uses the NWHFC algorithm developed by Dr. Chein-I Chang,
% see http://www.umbc.edu/rssipl/. The Matlab code was
% programmed by Jing Wang in Remote Sensing Signal and
% Image Processing Lab.
%

%% Authors Comments
%	Algorithm name: noise-whitened Harsanyi–Farrand–Chang (HFC) method.
%	Authors: Chein-I Chang
%	Category: preprocessing
%	Designed criteria: eigenvlaues of sample correlation/covariance matrix
%	Designed method: Neyman–Pearson detection theory
%	Typical use (LOI's addressed): estimation of number of spectral distinct signatures
%	Inputs: reflectance or radiance cube
%	Outputs: a positive integer and false alarm probability
%	Assumptions: no prior knowledge required
%	Sensitivity to LOI (target knowledge): moderate
%	Sensitivity to noise: low
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		The NWHFC method is a noise-whitened version of the HFC method and was developed in
%		Chang (2003a) and Chang and Du (2004). It requires a reliable method to estimate the noise
%		covariance matrix.


%% Code

[XX,YY,bnd] = size(HIM);
pxl_no = XX*YY;
r = (reshape(HIM,pxl_no,bnd))';
R = (r*r')/pxl_no;
u = mean(r,2);
K = R-u*u';
%======Noise estimation=====
K_Inverse=inv(K);
tuta=diag(K_Inverse);
K_noise=1./tuta;
K_noise=diag(K_noise);
%=====Noise whitening===
y=sqrtm(K_noise)\r;
y=reshape(y',XX,YY,bnd);
%=====Call HFC to estimate===
number=cchengHFC(y,t);