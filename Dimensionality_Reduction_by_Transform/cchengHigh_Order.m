function [ICs]=cchengHigh_Order(HIM,M,k,Initial)
%% Code Comments:
% HIM is the image cube.
% M is the number of components to generated using high order
% k is the order of statistics. For example, k=3, is the skewness, k=4, is
% the kurtosis, k=5 is the 5th moment, and so on...
% Initial condition preference, 0 is the random initial, 1 is the eigen
% initial, 2 is the unity initial. default is 0

%% Authors Comments
%	Algorithm name: high-order statistics DR (HOS-DR)
%	Authors: H. Ren and Chein-I Chang
%	Category: component analysis-based transform
%	Designed criteria: statistical moments higher than 2
%	Designed method: moment projection
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
%		ICA is widely used in blind source separation via a linear mixture model, which assumes that
%		there exists at most one Gaussian source. A fast algorithm developed by Hyvarinen and Oja
%		(1997), called FastICA, is the most used algorithm to implement ICA.
%
%	Furthere informationthat as given in a weird way:
%		Variance-based PCA and SNR-based MNF represent second-order statistics-based component
%		analysis transformations to perform DR. At the other extreme, ICA is an infinite-order statisticsbased
%		component analysis that makes use of mutual information to measure statistical independency
%		among all ICs. However, due to complexity of implementing mutual information, the commonly
%		used FastICA is actually developed by combining the third and fourth moments as a
%		criterion to measure the dependency among its generated components. So, technically speaking,
%		the FastICA-generated components are not really statistically independent. Instead, they can be
%		only considered as high-order statistically dependent. HOS-DR is developed to extend DR transformation
%		with order of statistics higher than 2 (Ren et al., 2006). In this context, the FastICA
%		(Hyvarinen and Ojha, 1997) can be considered as a special case of HOS-DR transformation.

%% Code

if nargin < 3
    fprintf('Please identify the order of statistics!');
    ICs=[];
else
    if nargin < 4
        Initial=0;
    end
    bnd=size(HIM,3);
    xx=size(HIM,1);
    yy=size(HIM,2);
    x=reshape(HIM,xx*yy,bnd);
    x=x';
    L=size(x,1);
    K=size(x,2);
    %===input x is a matrix with size=L*K;
    %====Sphering =====
    u=mean(x,2); %dimension of u is 1*L
    x_hat=x-u*ones(1,K);
    %===jing's code of cov
    m=mean(x,2);
    C=(x*x')/size(x,2)-m*m';
    %===========
    [V,D]=eig(C);
    A=inv(sqrtm(D))*V'; % A is the whitening matrix....
    x_whitened=A*(x_hat);
    clear x;
    clear x_hat;
    % Seperating , using high-order
    threshold = 0.01;
    B=zeros(L);
    y=x_whitened;
    W=ones(1,L);
    P_U_perl=eye(bnd);
    for round=1:M
        fprintf('IC %d', round);
        %===initial condition ===
        switch (Initial)
            case 0
                w=rand(L,1);
            case 1
                w=V(:,round); % Final version
            case 2
                w=ones(L,1);
            otherwise
                w=rand(L,1);
        end
        i=1;
        while i<=100 % maximum times of trying..
            a=(y.*repmat((w'*y).^(k-2),bnd,1))*y'; %skewness
            a=a/K; % get the sample mean as expectation
            [V,D]=eig(a);
            D=abs(D);
            [~,I]=max(diag(D));
            V1 = V(:,I);
            fprintf('.');
            distance(round,1,i)=norm(w-V1);
            distance(round,2,i)=norm(w+V1);
            if norm(w-V1)<threshold || norm(w+V1)<threshold
                fprintf('Convergence after %d steps\n', i);
                B(:,round)=w;
                W(round,:)=w';
                break;
            end
            w=V1;
            i=i+1;
        end
        %===if not converge. then use the results after 10 iterations
        B(:,round)=w;
        W(round,:)=w';
        %=======
        P_U_perl=eye(bnd)-W'*inv(W*W')*W;
        y=P_U_perl*y;
    end
    ICs=W*x_whitened;
    figure;
    for m=1:M
        s=reshape(abs(ICs(m,:)),xx,yy);
        s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
        -min(min(s)));
        temp=mean(reshape(s,xx*yy,1));
        subplot(5,6,m); imshow(uint8(s));
        %
    end
end
