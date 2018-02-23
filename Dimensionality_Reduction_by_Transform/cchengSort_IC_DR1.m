function [IC_sorted]=cchengSort_IC_DR1(HIM,M)
%% Code Comments:
%There were none, but here is my best interpretation
%
% There are two parameters,sort_IC_DR1(HIM,M) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% M is the number of independent components
%

%% Authors Comments
%	Algorithm name: ICA
%	Category: component analysis-based transform
%	Designed criteria: infinite-order statistics
%	Designed method: blind source separation via a linear mixture model
%	Typical use (LOI's addressed): no prior knowledge is required
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale images
%	Assumptions: no prior target knowledge required
%	Sensitivity to LOI (target knowledge): high
%	Sensitivity to noise: low
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		ICA (Hyvarinen et al., 2001) is widely used in blind source separation via a linear mixture
%		model, which assumes that there exists at most one Gaussian source. A fast algorithm developed
%		by Hyvarinen and Oja (1997), called FastICA, is the most used algorithm to implement ICA.
%		When ICA is implemented to perform DR, it suffers from a serious issue that the results are not
%		repeatable and inconsistent due to its use of random initial conditions. As a result, the ICA-generated
%		independent components (ICs) generally appear in a random order. In this case, the IC appearing
%		earlier does not necessarily imply that it is more important or significant than that generated
%		later. In order to resolve this issue, three versions of ICA, referred to as ICA-DR1,
%		ICA-DR2, and ICA-DR3, have been proposed by Wang and Chang (2006a, 2006b) and re-named
%		as statistics prioritized ICA-DR (SPICA-DR), random ICA (RICA-DR), and initialization driven
%		ICA-DR (IDICA-DR) in Sections 6.4.1–6.4.3, respectively. In all the three algorithms, an algorithm
%		developed by Hyvarinen and Oja (1997), called FastICA, is modified and implemented in
%		MATLAB codes to realize ICA. To distinguish it from the original FastICA, it is named “My
%		FastICA” and is provided as follows.

%% Code
clear J;
%   bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);
[ICs]=cchengFastICA_v5(HIM,M);
ICs=abs(ICs);
% %====Show the ICs in the original order===
% figure;
% for m=1:size(ICs,1)
%     s=reshape(abs(ICs(m,:)),xx,yy);
%     s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
%     -min(min(s)));
%     %   temp=mean(reshape(s,xx*yy,1));
%     subplot(6,8,m); imshow(uint8(s));
%     % title(m);
%     %
% end
%======Calculate the contrast function
for m=1:size(ICs,1)
    s=ICs(m,:);
    var1=var(s);
    mean1=mean(s);
    %   sita=sqrt(var1);
    skew_temp=sum((s-mean1).^3)/(xx*yy-1);
    kurt_temp=sum((s-mean1).^4)/(xx*yy-1);
    J(1,m)=(skew_temp.^2)/12+((kurt_temp-3).^2)/48;
end
%======IC sorting ========
[~,b]=sort(J);
b=flipud(b');
IC_sorted=ICs(b',:);

%=========Show the ICS after sorting...
figure;
for m=1:size(ICs,1)
    s=reshape(abs(IC_sorted(m,:)),xx,yy);
    s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
        -min(min(s)));
    %   temp=mean(reshape(s,xx*yy,1));
    subplot(6,8,m); imshow(uint8(s));
    % title(m);
    %
end

