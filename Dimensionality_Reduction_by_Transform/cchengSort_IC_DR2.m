function [IC_selected]=cchengSort_IC_DR2(HIM,M,run_times)
%% Code Comments:
%There were none, but here is my best interpretation
%
% There are three parameters,cchengSort_IC_DR2(HIM,M,run_times) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% M is the number of independent components
%run_times is the number of iterations or something ...

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
% bnd=size(HIM,3);
% xx=size(HIM,1);
% yy=size(HIM,2);
fprintf('first ICA run \n');
[ICs]=cchengFastICA_v5(HIM,M*2);
set1=abs(ICs);
set_com=set1;
size1=ones(1);
for round=1:run_times
    fprintf('ICA run, order is %d \n',round+1);
    [ICs]=cchengFastICA_v5(HIM,2*M);
    set2=abs(ICs);
    set_com_new=[];
    distance=0;
    for m=1:size(set_com,1)
        for n=1:size(set2,1)
            temp1=sqrt(sum(set_com(m,:).^2)); % SAM
            temp2=sqrt(sum(set2(n,:).^2));
            distance(m,n)=acos(sum(set_com(m,:).*set2(n,:))/(temp1*temp2));
        end
        t=distance(m,:)<=0.5;
        if (sum(t)>=1)
            set_com_new=cat(1,set_com_new,set_com(m,:));
        end
    end
    set_com=set_com_new;
    fprintf('the size of set_Com is');
    size(set_com_new)
    size1(round,1)=size(set_com_new,1);
    if (size(set_com_new,1)<=M)
        break;
    end
end
IC_selected=set_com;
