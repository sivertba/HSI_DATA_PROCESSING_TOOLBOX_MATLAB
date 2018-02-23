function [ab] = KNCLS(r,M,d)
% Non-negative constrain Least-square abundance estimator
% input M = signature matrix to be estimated.
% input r = image pixel vector.
% output abundance = abundance vector correspondence to the signature
% matrix.

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

% initial stage k = 0
k = 0;
[~,y] = size(M);
P = 1:y;
R = [];
%least square estimate of abundance
a_ls = kernelized(M,M,d,0)\kernelized(M,r,d,0);
%initially set abundance to least square abundance
ab = a_ls;
%iterate until all ab is positive
while any(ab<-0.5) && k<50
    % disp(k)
    % clc,
    k = k+1;
    % find negative index and move them from P to R
    neg_ind = find(ab<0);
    [P,R] = move_index(P,R,neg_ind);
    S = R;
    %initialize lum so the for loop can run at least once
    lum = 1;
    % iterate until all lum is smaller or equal to zero
    j=0;
    while any(lum>=0) && j<50
        % disp(j)
        j = j+1;
        a_R = a_ls(R);
        % define the steer matrix by removing row and columns defined by P
        a_steer_mat = inv_MTM;
        a_steer_mat(:,P) = [];
        a_steer_mat(P,:) = [];
        % calculate Lagrange multiplier lum
        lum = ((a_steer_mat)^(-1))*a_R;
        % if lum contains positive value
        if any(lum>0) && j<1000
            % find the maximum lum move its index from R to P create a new
            % lum by removing the max lum
            j = j+1;
            lum_max_ind = R(lum == max(lum));
            [R,P]=move_index(R,P,lum_max_ind);
            a_steer_mat = inv_MTM;
            a_steer_mat(:,P) = [];
            a_steer_mat(P,:) = [];
            a_R = a_ls(R);
            lum_new = ((a_steer_mat)^(-1))*a_R;
            % form another lum steering matrix
            lum_steer_mat = inv_MTM;
            lum_steer_mat(:,P) = [];
            if isempty(lum_steer_mat)
                lum_steer_mat = 0;
            end
            % calcuate the abundance base on R and P and lum_new. If any
            % value specify by index S is negative move it from P to R.
            a_S = a_ls - lum_steer_mat*lum_new;
            s_neg_ind = S(a_S(S)<0);
            [P,R]= move_index(P,R,s_neg_ind);
        end
    end
    % calculate the new abundance base on the new lum found
    lum_steer_mat = inv_MTM;
    lum_steer_mat(:,P) = [];
    ab = a_ls - lum_steer_mat*lum;
end
return;
% Move index function
function [P_new,R_new] = move_index(del_ind, get_ind, move_ind)
n = length(move_ind);
if n>0
    for i = 1:n
        a = del_ind == move_ind(i);
        del_ind(a)=[];
        b = get_ind - move_ind(i);
        if any(b == 0)
            get_ind = get_ind;
        else
            get_ind = [get_ind,move_ind(i)];
        end
    end
end
P_new = sort(del_ind);
R_new = sort(get_ind);
return
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
