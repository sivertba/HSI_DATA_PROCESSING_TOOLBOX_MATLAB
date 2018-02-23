function [PCs]=cchengPCA(HIM,M)
%% Code Comments:
%PCA does what PCA does best
%Should consider to sue the built-in MATLAB PCA version
%
% There are two parameters,cchengPCA(HIM,M) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% M appears to be the number of components or what not ...
%
%
%% Authors Comments
%	Algorithm name: principal components analysis (PCA)
%	Category: component analysis-based transform
%	Designed criteria: eigenvalues
%	Designed method: sample covariance matrix
%	Typical use (LOI’s addressed): data representation by eigenvectors
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale images
%	Assumptions: no prior target knowledge required
%	Sensitivity to LOI (target knowledge): targets of high-order statistics
%	Sensitivity to Noise: low
%	Operating Bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		PCA (Gonzalez and Woods, 2002) is probably the most widely used component analysis transform
%		that allows users to data in terms of eigenvalues/eigenvectors via eigen-decomposition
%		where data are preserved and retained according to data variances in the descending order.


%% Code
bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);
x=reshape(HIM,xx*yy,bnd);
x=x';
L=size(x,1);
K=size(x,2);
u=mean(x,2); %dimension of u is 1*L
x_hat=x-u*ones(1,K);
m=mean(x,2);
C=(x*x')/size(x,2)-m*m';
%===========
[V,D]=eig(C);
d=(diag(D))';
[~,Index]=sort(d);
for m=1:L
    D_sort(1,m)=d(1,Index(1,Lm));
    V_sort(:,m)=V(:,Index(1,Lm));
end
D=diag(D_sort);
V=V_sort;
D=D(1:M,1:M);
V=V(:,1:M);
%====for the matrix with full column rank, so the
A=V';
x_whitened=A*(x_hat);
PCs=x_whitened;

