function results = FCLS_v2(image,M,tol)

%% Code Comments
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% – Fully Constrain Least-Squares (FCLS) abundances estimate –
%
% function: results = FCLS(image,M,tol)
%
% input: image = image of size [X,Y,Z]
% M = endmember matrix of size [Z,p]
% tol = NCLS tolerance, e.g. -1e-6
%
% output: results = resulting abundance map of size [X,Y,p]
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

%% Authors Comments
%	Algorithm name: fully constrained least squares (FCLS)
%	Authors: Daniel Heinz and Chein-I Chang
%	Category: least squares error-based spectral filter
%	Designed criteria: least squares error (LSE)
%	Designed Method: a priori, supervised and fully constrained LSMA
%	Typical use (LOI's addressed): detection, classification, discrimination, identification
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

[x,y,z] = size(image);
num_p = size(M,2);
results = zeros(x,y,num_p);
for i = 1:x
    for j = 1:y
        r = reshape(image(i,j,:),z,1);
        delta = 1/(10*max(max(M)));
        s = [delta.*r;1];
        N = [delta.*M;ones(1,num_p)];
        [ab] = NCLS(s, N, tol);
        ab = reshape(ab,[1,1,num_p]);
        results(i,j,:) = ab;
    end
end
return;

function [abundance]=NCLS(x, MatrixZ, tol)
% input MatrixZ is the signatures of endmembers. It is of size [bands p].
% input x is the signature whose abundance is to be estimated.
% output abundance is the abundance of each material in r1. It is of size [p 1].
% This function is written according to Dr. Chang's first book , P 47
M=size(MatrixZ,2);
R=zeros(M,1);
P=ones(M,1);
invMtM=(MatrixZ'*MatrixZ)^(-1);
Alpha_ls=invMtM*MatrixZ'*x;
Alpha_ncls=Alpha_ls;
min_Alpha_ncls=min(Alpha_ncls);
j=0;
while(min_Alpha_ncls<-tol && j<500)
    j = j+1;
    for II=1:M
        if((Alpha_ncls(II)<0)&&(P(II)==1))
            R(II)=1;
            P(II)=0;
        end %%% end of if (Alpha_ncls(II)<0)
    end % end of for II=1:M
    S = R;
    goto_step6=1;
    counter = 0;
    while(goto_step6==1)
        index_for_Lamda = find(R==1);
        Alpha_R = Alpha_ls(index_for_Lamda);
        Sai = invMtM(index_for_Lamda,index_for_Lamda);
        inv_Sai = (Sai)^(-1); % remember inversion of Sai
        Lamda=inv_Sai*Alpha_R;
        [max_Lamda,index_Max_Lamda]=max(Lamda);
        counter = counter+1;
        if ( max_Lamda<=0 || counter == 200 )
            break;
        end
        temp_i = inv_Sai; % simplify the inversion of matrix
        temp_i(1,:) = inv_Sai(index_Max_Lamda,:);
        if index_Max_Lamda>1
            temp_i(2:index_Max_Lamda,:) = inv_Sai(1:index_Max_Lamda-1,:);
        end
        inv_Sai_ex = temp_i;
        inv_Sai_ex(:,1) = temp_i(:,index_Max_Lamda);
        if index_Max_Lamda>1
            inv_Sai_ex(:,2:index_Max_Lamda) = temp_i(:,1:index_Max_Lamda-1);
        end
        inv_Sai_next = inv_Sai_ex(2:end,2:end)...
            - inv_Sai_ex(2:end,1)*inv_Sai_ex(1,2:end)/inv_Sai_ex(1,1);
        P(index_for_Lamda(index_Max_Lamda))=1;
        R(index_for_Lamda(index_Max_Lamda))=0;
        index_for_Lamda(index_Max_Lamda) = [];
        Alpha_R = Alpha_ls(index_for_Lamda);
        Lamda=inv_Sai_next*Alpha_R;
        Phai_column = invMtM(:,index_for_Lamda);
        if (size(Phai_column,2)~=0)
            Alpha_s=Alpha_ls-Phai_column*Lamda;
        else
            Alpha_s=Alpha_ls;
        end
        goto_step6=0;
        for II=1:M
            if ((S(II)==1)&&(Alpha_s(II)<0))
                P(II)=0;
                R(II)=1;
                goto_step6=1;
            end
        end
    end % end of while (gotostep6==1)
    index_for_Phai = R==1;
    Phai_column = invMtM(:,index_for_Phai);
    if (size(Phai_column,2)~=0)
        Alpha_ncls=Alpha_ls-Phai_column*Lamda;
    else
        Alpha_ncls=Alpha_ls;
    end
    min_Alpha_ncls=min(Alpha_ncls);
end % end of while
abundance=zeros(M,1);
for II=1:M
    if (Alpha_ncls(II)>0)
        abundance(II)=Alpha_ncls(II);
    end
end
return;

