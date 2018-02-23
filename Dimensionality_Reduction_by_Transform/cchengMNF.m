function [HIM_MNF,Matrix_of_Vector]=cchengMNF(HIM)
%% Code Comments:
%There were none but I would assume that it safe to say that
% There is one parameter,MNF(HIM) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
%

%% Authors Comments
%	Algorithm name: maximum noise fraction (MNF)
%	Authors: A.A. Green, M. Berman, P. Switzer, and M.D. Craig
%	Category: component analysis-based transform
%	Designed criteria: signal-to-noise ratio
%	Designed method: sample covariance matrix
%	Typical use (LOI’s addressed): data representation by SNR
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale images
%	Assumptions: no prior target knowledge required
%	Sensitivity to LOI (target knowledge): high-order statistics
%	Sensitivity to Noise: high
%	Operating Bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		MNF was developed by Green et al. (1988) and further modified by Lee et al. (1990) and
%		referred to as Noise Adjusted Principal Component (NAPC) transform. It represents data in
%		terms of SNR rather than data variances as PCA does.


%% Code
[height,width,NumberOfSpectrum]=size(HIM);
HIM_MNF=zeros(height,width,NumberOfSpectrum);
% ————————————————— begin to compute Matrix_of_Vector ———————————————— %
meanSpect=zeros(NumberOfSpectrum,1);
for II=1:height
    for JJ=1:width
        meanSpect=meanSpect+squeeze(HIM(II,JJ,:))/height/width;
    end
end
TotalCovariance=zeros(NumberOfSpectrum,NumberOfSpectrum);
for II=1:height
    for JJ=1:width
        TotalCovariance=TotalCovariance+(squeeze(HIM(II,JJ,:))- ...
        meanSpect)*(squeeze(HIM(II,JJ,:))-meanSpect)'/height/width;
    end
end
Matrix_F=zeros(NumberOfSpectrum,NumberOfSpectrum);
Cov_inv=inv(TotalCovariance);
for II=1:NumberOfSpectrum
    Matrix_F(II,II)=sqrt(Cov_inv(II,II));
end
adjusted_Cov=Matrix_F'*TotalCovariance*Matrix_F;
[V,D]=eig(adjusted_Cov);
eig_value=zeros(NumberOfSpectrum,2);
for II=1:NumberOfSpectrum
    eig_value(II,1)=D(II,II);
    eig_value(II,2)=II;
end
%disp(eig_value);
V_sort_min_to_max=sortrows(eig_value,1);
Matrix_of_Vector_before=zeros(NumberOfSpectrum,NumberOfSpectrum);
for II=1:NumberOfSpectrum
    Matrix_of_Vector_before(:,II)=squeeze(V(:,V_sort_min_to_max(NumberOfSpectrum-II+1,2)));
end
Matrix_of_Vector=Matrix_F*Matrix_of_Vector_before;
% ——————————————— end of computing Matrix_of_Vector ———————————————— %
for II=1:height
    for JJ=1:width
        r=squeeze(HIM(II,JJ,:))-meanSpect;
        HIM_MNF(II,JJ,:)=Matrix_of_Vector'*r;
    end
end