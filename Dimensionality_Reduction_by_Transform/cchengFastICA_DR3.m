function [ICs]=cchengFastICA_DR3(HIM,M)
%% Code Comments:
%There were none, but here is my best interpretation
%
% There are two parameters,cchengFastICA_v5(HIM,M) where the HIM is the
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
bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);
x=reshape(HIM,xx*yy,bnd);
x=x';
L=size(x,1);
K=size(x,2);
%====Sphering =====
u=mean(x,2); %dimension of u is 1*L
x_hat=x-u*ones(1,K);
m=mean(x,2);
C=(x*x')/size(x,2)-m*m';
%===========
[V,D]=eig(C);
A=sqrtm(D)\V'; % A is the whitening matrix....
x_whitened=A*(x_hat);
%====for cuprite data===
clear x;
clear x_hat;
%=========for initialization
[~,Sig]=cchengATGP(reshape(x_whitened',xx,yy,L),M);
W_initial=Sig;
threshold = 0.0001;
B=zeros(L);
%===============find the first point=========
for round=1:M
    fprintf('IC %d', round);
    %===Initial condition
    w=W_initial(:,round);
    %===
    w=w-B*B'*w;
    w=w/norm(w);
    wOld=zeros(size(w));
    wOld2=zeros(size(w));
    i=1;
    while i<=1000
        w=w-B*B'*w;
        w=w/norm(w);
        fprintf('.');
        if norm(w-wOld)<threshold || norm(w+wOld)<threshold
            fprintf('Convergence after %d steps\n', i);
            B(:,round)=w;
            W(round,:)=w';
            break;
        end
        wOld2=wOld;
        wOld=w;
        w=(x_whitened*((x_whitened'*w).^3))/K-3*w;
        w=w/norm(w);
        i=i+1;
    end
    if (i>1000)
        fprintf('Warning! cannot converge after 1000 steps \n, no more components');
        break;
    end
end
ICs=W*x_whitened;

figure
for m=1:M
    s=reshape(abs(ICs(m,:)),xx,yy);
    s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
    -min(min(s)));
    %   temp=mean(reshape(s,xx*yy,1));
    subplot(5,6,m); imshow(uint8(s));
    %
end
% A=W;
% x_re_ICA=V*sqrtm(D)*(A*ICs)+u*ones(1,xx*yy);