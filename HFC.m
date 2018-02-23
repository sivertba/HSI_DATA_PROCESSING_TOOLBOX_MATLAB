function number=HFC(HIM,t);
%
% HFC gives the VD number estimated by given false alarm property using HFC
% method.
%
% There are two parameters,HFC(HIM,t) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube.
% t is the false alarm probability.
%
% HFC uses the HFC algorithm developed by Dr. Chein-I Chang,
% see http://www.umbc.edu/rssipl/. The Matlab code was
% programmed by Jing Wang in Remote Sensing Signal and
% Image Processing Lab.
%
[XX,YY,bnd] = size(HIM);
pxl_no = XX*YY;
r = (reshape(HIM,pxl_no,bnd))’;
R = (r*r’)/pxl_no;
u = mean(r,2);
K = R-u*u’;
%======HFC=====
D1=sort(eig(R));
D2=sort(eig(K));
sita=sqrt((D1.^2+D2.^2)*2/pxl_no);
P_fa=t;
Threshold=(sqrt(2))*sita*erfinv(1-2*P_fa);
Result=0;
for m=1:bnd
if ((D1(m,1)-D2(m,1))>Threshold(m,1))
Result(m,1)=1;
else
Result(m,1)=0;
end
end
fprintf(’The VD number estimated is’);
disp(sum(Result));
number=sum(Result);
