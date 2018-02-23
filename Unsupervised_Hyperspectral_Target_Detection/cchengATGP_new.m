function [loc,Sig]=cchengATGP_new(HIM,num_targets)
%% Code Comments:
%
%
% input HIM is the Hyperspectral data cube of size [height width bands].
% input M is the number of target p that you need to extract.
% output Sig is the signatures corresponding to the extracted targets. It is of size [bands p].
% output loc is the positions of the targets. It is of size [p 2]
% with the first column being the vertical position, the second column being the horizontal position.


%% Authors Comments
%	Algorithm name: automatic target generation process (ATGP)
%	Authors: H. Ren and Chein-I Chang
%	Category: unsupervised least squares-based filter
%	Designed criteria: orthogonal projection (OP)
%	Designed method: maximum OP
%	Typical use (LOI's addressed): detection, classification, discrimination, identification
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale abundance fractional images
%	Assumptions: no target knowledge required
%	Sensitivity to LOI (target knowledge): high
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, real time
%	Brief description:
%		The ATGP is derived from an algorithm, called automatic target detection, and classification
%		(ATDCA) developed by Ren and Chang to find target pixels of interest for recognition without
%		having prior target knowledge (Ren and Chang, 2003). It performs successive orthogonal projections
%		to find a data sample vector, considered as a target of interest, which has the maximal projection
%		after each orthogonal projection. The potential of ATGP has been shown to have a wide
%		range of applications, such as VD determination, anomaly detection, endmember extraction,
%		unsupervised LSMA, etc.


%% Code

bnd=size(HIM,3);
xx=size(HIM,1);
yy=size(HIM,2);
r=reshape(HIM,xx*yy,bnd);
r=r';
%=====Find the first point
temp=sum(r.*r);
[~,b]=max(temp);
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
%==========
for m=2:num_targets
    U=Sig;
    P_U_perl=eye(bnd)-U*inv(U'*U)*U';
    y=P_U_perl*r;
    temp=sum(y.*y);
    [~,b]=max(temp);
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
    %disp(m)
end

% figure; imagesc(HIM(:,:,30)); colormap(gray); hold on
% axis off
% axis equal
% for m=1:size(Loc,1)
%     plot(Loc(m,1),Loc(m,2),'o','color','g');
%     text(Loc(m,1)+2,Loc(m,2),num2str(m),'color','y','FontSize',12);
% end

loc(:,1)=Loc(:,2);
loc(:,2)=Loc(:,1);