function temp=KOSP(image,d,U,sig)
%% Code Comments:
% Kernel based LSOSP function
% Input:
% image = image cube input
% d = desired signature, example: [2;3;4]
% U = undesired signature matrix
% sig = parameter that control RBF kernel function
% output:
% temp = resulting map
%

%% Authors Comments
%	Algorithm name: orthogonal subspace projection (OSP)
%	Authors: J.C. Harsanyi and Chein-I Chang
%	Category: spectrally matched filter
%	Designed Criteria: SNR ratio
%	Designed Method: A priori, supervised and unconstrained least squares-based linear spectral mixture analysis
%	Typical use (LOI's addressed): detection, classification, discrimination, identification
%	Inputs: reflectance or radiance cube, complete target knowledge
%	Outputs: gray-scale abundance fractional images
%	Assumptions: complete prior target knowledge required
%	Sensitivity to LOI (target knowledge): high
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, real Time
%	Brief description:
%		The OSP approach was first developed in Harsanyi's Ph.D. dissertation in 1993 (Harsanyi, 1993)
%		and was later published in IEEE Transaction on Geoscinece and Remote Sensing, July 1994
%		(Harsanyi and Chang, 1994). It is a linear unmixng method that takes advantage of a linear mixture
%		model to detect, classify, and identify targets of interest. The idea is to separate target sources
%		into desired and undesired targets and then use an orthogonal project to reject the desired
%		targets before a matched filtration takes place. So, it operates two functions in sequence: an
%		undesired target rejecter followed by a spectral matched filtration. Since OSP was originally
%		designed as a detector and cannot accurately estimate signature abundance fractions, a least
%		squares version of OSP, referred to as least squares OSP (LSOSP) was further developed in
%		Chang et al. (1998b) for abundance fraction estimation. The only difference between OSP and
%		LSOSP is that LSOSP includes a normalization constant to account for estimation error incurred
%		in the OSP-derived detector (Chang, 2009). So, the MATLAB codes provided in the following is
%		a more general version of OSP, LSOSP

%% Code

[x, y, z]=size(image);
temp = zeros(x,y);
% perform least squares-based estimator on all image vectors
KdU = kernelized(d,U,sig,0);%disp(KdU),
KUU = kernelized(U,U,sig,0);%disp(KUU),
Kdd = kernelized(d,d,sig,0);
KUd = kernelized(U,d,sig,0);%disp(KUd),
for i = 1:x
    for j = 1:y
        for k = 1:z
            r(k,1) = image(i,j,k);
        end
        Kdr = kernelized(d,r,sig,0);%disp(Kdr),
        KUr = kernelized(U,r,sig,0);%disp(KUr),
        temp(i,j)=(Kdr-KdU*inv(KUU)*KUr);%/(Kdd-KdU*inv(KUU)*KUd);
    end
end
%% kernelization function
function results = kernelized(x,y,d,chk)
x_l = size(x,2);
y_l = size(y,2);
results = zeros(x_l,y_l);
for i = 1:x_l
    for j = 1:y_l
        results(i,j)= exp((-1/2)*(norm(x(:,i)-y(:,j))^2)/(d^2));
        %RBF kernel (can be changed)
    end
end
if chk == 1
    results = results-(sum(sum(results))/(x_l*y_l))*ones(x_l,y_l);
elseif chk == 2
    N = (1/(x_l*y_l))*ones(x_l,y_l);
    results = results-N*results-results*N+N*results*N;
end
