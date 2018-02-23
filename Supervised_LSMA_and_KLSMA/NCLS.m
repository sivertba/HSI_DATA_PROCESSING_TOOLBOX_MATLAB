function [abundance,error_vector]=NCLS(MatrixZ,r1)
%% Code Comments:
% input MatrixZ is the signatures of endmembers. It is of size [ bands p].
% input x is the signature whose abundance is to be estimated.
% output abundance is the abundance of each material in r1. It is of size [p 1].
% output error_vector is the error vector of size [bands 1].
% This function is written according to Dr. Chang's first book , P 47
%

%% Authors Comments
%	Algorithm name: nonnegativity least squares (NCLS)
%	Authors: Chein-I Chang and Daniel Heinz
%	Category: least squares error-based spectral filter
%	Designed criteria: least squares error
%	Designed method: A priori, supervised and ANC-constrained LSMA
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
%		Nonnegativity constrained least squares (NCLS) was developed to improve OSP in signal detection
%		(Chang and Heinz, 2000b). It imposes ANC on the linear mixing model to make sure that
%		the unmixed abundance fractions are nonnegative. Despite not being fully constrained, on many
%		occasions NCLS can perform abundance estimation as well as a fully abundance-constrained
%		method in. However, in signal detection, NCLS generally performs better than unconstrained
%		and fully constrained methods.

%% Code

x=r1; %rename r1 as x;
M=size(MatrixZ,2);
count_R=0;
count_P=M;
R=zeros(M,1);
P=ones(M,1);
%tolerance=0.000001;
d=zeros(M,1);
Alpha_ls=inv(MatrixZ'*MatrixZ)*MatrixZ'*x;
Alpha_ncls=Alpha_ls;
min_Alpha_ncls=min(Alpha_ncls);
M_t_r=MatrixZ'*x;
invMtM=inv(MatrixZ'*MatrixZ);
while(min_Alpha_ncls<-0.000000001)
    for II=1:M
        if((Alpha_ncls(II)<0)&&(P(II)==1))
            R(II)=1;
            P(II)=0;
        end %%% end of if (Alpha_ncls(II)<0)
    end % end of for II=1:M
    S=R;
    goto_step6=1;
    while(1)
        sum_R=sum(R);
        Alpha_R=zeros(sum_R,1);
        count_for_Alpha_R=0;
        for II=1:M
            if (R(II)==1)
                count_for_Alpha_R=count_for_Alpha_R+1;
                Alpha_R(count_for_Alpha_R)=Alpha_ls(II);
                index_for_Lamda(count_for_Alpha_R)=II;
            end
        end
        count_1_for_P=0;
        Sai_column=[];
        for II=1:M
            if (P(II)~=1)
                Sai_column=[Sai_column squeeze(invMtM(:,II)) ];
            end
        end
        Sai=[];
        for II=1:M
            if (P(II)~=1)
                Sai=[Sai
                    squeeze(Sai_column(II,:)) ];
            end
        end
        Lamda=inv(Sai)*Alpha_R;
        if(max(Lamda)<0)
            break;
        end
        [max_Lamda,index_Max_Lamda]=max(Lamda);
        P(index_for_Lamda(index_Max_Lamda))=1;
        R(index_for_Lamda(index_Max_Lamda))=0;
        sum_R=sum(R);
        Alpha_R=zeros(sum_R,1);
        count_for_Alpha_R=0;
        for II=1:M
            if (R(II)==1)
                count_for_Alpha_R=count_for_Alpha_R+1;
                Alpha_R(count_for_Alpha_R)=Alpha_ls(II);
                index_for_Lamda(count_for_Alpha_R)=II;
            end
        end
        Sai_column=[];
        for II=1:M
            if (P(II)~=1)
                Sai_column=[Sai_column squeeze(invMtM(:,II)) ];
            end
        end
        Sai=[];
        for II=1:M
            if (P(II)~=1)
                Sai=[Sai
                    squeeze(Sai_column(II,:)) ];
            end
        end
        Lamda=inv(Sai)*Alpha_R;
        Phai_column=[];
        for II=1:M
            if (P(II)~=1)
                Phai_column=[Phai_column squeeze(invMtM(:,II)) ];
            end
        end
        if (size(Phai_column,2)~=0)
            Alpha_s=Alpha_ls-Phai_column*Lamda;
        else
            Alpha_s=Alpha_ls;
        end
        goto_step6=0;
        find_smallest_in_S=zeros(M,2);
        find_smallest_in_S(:,1)=Alpha_s;
        find_smallest_in_S(:,2)=[1:M]';
        sort_find=sortrows(find_smallest_in_S,1);
        for II=1:M
            if ((S(II)==1)&(Alpha_s(II)<0))
                P(II)=0;
                R(II)=1;
                goto_step6=1;
            end
        end
    end % end of while (gotostep6==1)
    Phai_column=[];
    for II=1:M
        if (P(II)~=1)
            Phai_column=[Phai_column squeeze(invMtM(:,II)) ];
        end
    end
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
error_vector=MatrixZ*abundance-x;