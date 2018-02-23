function [location_pair,spectrum_of_targets]=cchengUNCLS(ImageCub, NumberOfClass)
%% Code Comments:
%
% input ImageCub is of size [height width bands].
% input NumberOfClass is the number of target p that you need to extract.
% output spectrum_of_targets is the signatures corresponding to the extracted 
% targets. It is of size [bands p].
% output location_pair is the positions of the targets. It is of size [p 2] 
% with the first column being the vertical position, the second column 
% being the horizontal position.

%% Authors Comments
%	Algorithm name: unsupervised non-negativity constrained least squares (UNCLS)
%	Authors: Chein-I Chang and Daniel Heinz
%	Category: unsupervised least-squares-based filter
%	Designed criteria: least squares error
%	Designed method: least squares constrained method
%	Typical use (LOI's addressed): detection, spectral unmixing, classification, quantification
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale abundance fractional images
%	Assumptions: no target knowledge required
%	Sensitivity to LOI (target knowledge): moderate
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		The UNCLS is an unsupervised version of the NCLS developed by Chang and Heinz (2000b).
%		Unlike NCLS, which requires complete a priori target knowledge to perform spectral unmixing,
%		UNCLS performs unsupervised target detection by using the NCLS to generate a set of potential
%		targets directly from the data to be processed. So, its main goal is to find targets of interest
%		without any prior knowledge. Because of that, UNCLS is primarily used for target detection and
%		endmember extraction as opposed to NCLS, whose main functionality is spectral unmixing.
%		That is, the NCLS and the UNCLS have rather different applications.


%% Code

NumberOfClass=NumberOfClass-1;
%[height, width, NumberOfSpectrum]=size(ImageCub);
[xx, yy, bnd]=size(ImageCub);
r=reshape(ImageCub,xx*yy,bnd);
r=r';
%height=size(ImageCub,1);
%width=size(ImageCub,2);
%MatrixH=zeros( NumberOfSpectrum, NumberOfClass+1);
location_pair=zeros(NumberOfClass,2);
Brightest=10^(-100);
count=0;
NumberOfExisted=1;
%=====Find the first point=====
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
%location_pair(1:NumberOfExisted,:)=squeeze(BrightForEachClass(1,:));
location_pair(:,1)=Loc(:,2);
location_pair(:,2)=Loc(:,1);
MatrixH=r(:,b);
for II=1:NumberOfClass
    % disp([ ' Class : ' int2str(II) ]);
    % disp(II+NumberOfExisted);
    minError=0.0000000000000000000000001;
    for b=1:xx*yy
        trypixel=r(:,b);
        max_value=max(trypixel);
        if (max_value>0)
            [abundance]= NCLS(MatrixH,trypixel);
            error_vector=MatrixH*abundance-trypixel;
            error=error_vector'*error_vector;
            if(error>minError)
                minError=error;
                if (rem(b,xx)==0)
                    location_pair(II+NumberOfExisted,1)=xx;
                    location_pair(II+NumberOfExisted,2)=b/xx;
                elseif (floor(b/xx)==0)
                    location_pair(II+NumberOfExisted,1)=b;
                    location_pair(II+NumberOfExisted,2)=1;
                else
                    location_pair(II+NumberOfExisted,1)=b-xx*floor(b/xx);
                    location_pair(II+NumberOfExisted,2)=floor(b/xx)+1;
                end
            end
        end
    end
    MatrixH=[MatrixH
        squeeze(ImageCub(location_pair(II+NumberOfExisted,1),...
            location_pair(II+NumberOfExisted,2),:))];
end
spectrum_of_targets=MatrixH;
