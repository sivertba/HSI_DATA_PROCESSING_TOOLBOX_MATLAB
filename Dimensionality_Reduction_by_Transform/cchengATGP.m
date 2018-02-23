function [Loc,Sig]=cchengATGP(HIM,M)
%% Code Comments:
%There were none, but here is my best interpretation
%
% There are two parameters,cchengATGP(HIM,M) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% M is the number of independent components

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

bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);
r=reshape(HIM,xx*yy,bnd);
r=r';
%=====Find the first point
temp=sum(r.*r);
[a,b]=max(temp);
if (rem(b,xx)==0)
    Loc(1,1)=b/xx;
    Loc(1,2)=xx;
elseif (floor(b/xx)==0)
    Loc(1,1)=1;
    Loc(1,2)=b;
else
    Loc(1,1)=floor(b/xx)+1; % y
    Loc(1,2)=b-xx*floor(b/xx); % x
end
Sig(:,1)=r(:,b);
fprintf('1\n');
%==========
for m=2:M
    U=Sig;
    P_U_perl=eye(bnd)-U*inv(U'*U)*U';
    y=P_U_perl*r;
    temp=sum(y.*y);
    [a,b]=max(temp);
    if (rem(b,xx)==0)
        Loc(m,1)=b/xx;
        Loc(m,2)=xx;
    elseif (floor(b/xx)==0)
        Loc(m,1)=1;
        Loc(m,2)=b;
    else
        Loc(m,1)=floor(b/xx)+1; % y
        Loc(m,2)=b-xx*floor(b/xx); % x
    end
    Sig(:,m)=r(:,b);
    disp(m)
end
%
% figure; imagesc(HIM(:,:,30)); colormap(gray); hold on
% axis off
% axis equal
% for m=1:size(Loc,1)
% plot(Loc(m,1),Loc(m,2),'o','color','g');
% text(Loc(m,1)+2,Loc(m,2),num2str(m),'color','y','FontSize',12);
% end

