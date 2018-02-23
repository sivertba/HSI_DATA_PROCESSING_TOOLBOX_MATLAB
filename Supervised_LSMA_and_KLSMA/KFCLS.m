function [abundance,error_vector]= KFCLS(M,r1,delta,d)
%% Code Comments:
% input M is the signatures of endmembers. It is of size [ bands p].
% input r1 is the signature whose abundance is to be estimated.
% input delta is the parameter to control ASC: see explanation on Dr. Chang's first book , P 183
% usually, delta is 1/(10*max(max(A));
% output abundance is the abundance of each material in r1. It is of size [p 1].
% output error_vector is the error vector of size [bands 1].
%

%% Authors Comments
%	Algorithm name: fully constrained least squares (FCLS)
%	Authors: Daniel Heinz and Chein-I Chang
%	Category: least squares error-based spectral filter
%	Designed criteria: least squares error (LSE)
%	Designed Method: a priori, supervised and fully constrained LSMA
%	Typical use (LOI’s addressed): detection, classification, discrimination, identification
%	Inputs: reflectance or radiance cube, complete target knowledge
%	Outputs: gray-scale abundance fractional images
%	Assumptions: complete prior target knowledge required
%	Sensitivity to LOI (target knowledge): high
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, real time
%	Brief description:
%		Fully constrained least squares (FCLS) was developed in Heinz and Chang (2001) as an abundance
%		estimator to accurately estimate abundance fractions. Therefore, as far as unmixing is
%		concerned, FCLS is one of best designed linear estimator. But it does not imply that FCLS is
%		also best for other applications such as detection, discrimination, classification, etc. As a matter
%		of fact, NCLS generally outperforms FCLS in these applications, where NCLS does not comply
%		with the constraint of ASC, specifically, in signal detection, where signals to be detected are
%		corrupted by noise and ASC is certainly violated.


%% Code

tic
%Dan's: optimal
A=M;
numloop=size(A,2);
e=delta;
eA=e*A;
E=[ones(1,numloop);eA];
EtE=kernelized(E,E,d,0);
[m,~] = size(EtE);
One=ones(m,1);
iEtE=inv(EtE);
iEtEOne=iEtE*One;
sumiEtEOne=sum(iEtEOne);
weights=diag(iEtE);
c=0;
sample=r1;
er=e*sample;
f=[1;er];
Etf=kernelized(E,f,d,0);
tol=1e-7;
%fcls1a
%%%% THIS IS lamdiv2
ls=iEtE*Etf;
lamdiv2=-(1-(ls'*One))/sumiEtEOne;
x2=ls-lamdiv2*iEtEOne;
x2old=x2;
if (any(x2<-tol))
    Z=zeros(m,1);
    iter=0;
    while(any(x2<-tol) && iter >(m))
        Z(x2<-tol)=1;
        zz=find(Z);
        x2=x2old; % Reset x2
        L=iEtE(zz,zz);
        ab=size(zz);
        lastrow=ab(1)+1;
        lastcol=lastrow;
        L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';
        L(1:ab(1),lastcol)=iEtEOne(zz);
        L(lastrow,lastcol)=sumiEtEOne;
        xerow=x2(zz);
        xerow(lastrow,1)=0;
        lagra=L\xerow;
        while (any(lagra(1:ab(1))>0)) % Reset Lagrange multipliers
            maxneg=weights(zz).*lagra(1:ab(1));
            [~,iz]=max(maxneg); % Remove the most positive
            Z(zz(iz))=0;
            zz=find(Z); % Will always be at least one (prove)
            L=iEtE(zz,zz);
            ab=size(zz);
            lastrow=ab(1)+1;
            lastcol=lastrow;
            L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';
            L(1:ab(1),lastcol)=iEtEOne(zz);
            L(lastrow,lastcol)=sumiEtEOne;
            xerow=x2(zz);
            xerow(lastrow,1)=0;
            lagra=L\xerow;
        end
        %problem with lamscls zz may be null
        if ~isempty(zz)
            x2=x2-iEtE(:,zz)*lagra(1:ab(1))-lagra(lastrow)*iEtEOne;
        end
        iter=iter+1;
    end
end
abundance=x2;
error_vector=A*abundance-r1;
function results = kernelized(x,y,d,chk)
x_l = size(x,2);
y_l = size(y,2);
results = zeros(x_l,y_l);
for i = 1:x_l
    for j = 1:y_l
        results(i,j)= exp((-1/2)*(norm(x(:,i)-y(:,j))^2)/(d^2));
    end
end
if chk == 1
    results = results-(sum(sum(results))/(x_l*y_l))*ones(x_l,y_l);
elseif chk == 2
    N = (1/(x_l*y_l))*ones(x_l,y_l);
    results = results-N*results-results*N+N*results*N;
end
