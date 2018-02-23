function sphered_data=cchengData_Sphering(HIM)
%% Code Comments:
% For information about data_sphering see README
%
% There is one parameter, data_sphering(HIM) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
%

%% Authors Comments
%	Algorithm name: data sphering
%	Category: preprocessing
%	Designed criteria: eigenvalues
%	Designed method: sample covariance matrix
%	Typical use (LOI’s addressed): preprocessing
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale images
%	Assumptions: no prior knowledge required
%	Sensitivity to LOI (target knowledge): high to high-order statistics
%	Sensitivity to noise: low
%	Operating Bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		Data sphering, also known as a whitening processing in statistical signal processing and
%		communication, is a commonly used method to remove data sample vectors characterized
%		by the first two-order statistics. It is a required preprocessing step prior to ICA (Hyvarinen
%		and Oja, 2001).
%


%% Code
% Initial variables
[row, column, band]=size(HIM);
% Center the data
HIM_row=reshape(HIM,row*column,band);
HIM_zeromean=HIM_row-repmat(mean(HIM_row),row*column,1);
cov=HIM_zeromean'*HIM_zeromean/(row*column);
% Eigen decomposition
[V, D]=eig(cov);
% Transform the data set
for i=1:band
    % TODO: Preallocate data
    sphered_data(i,:)=(V(:,i)'*HIM_zeromean')./(D(i,i)^.5*row);
end
% Transform the data back
sphered_data=reshape(sphered_data',row,column,band);