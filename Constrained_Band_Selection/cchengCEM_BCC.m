function [ band_select, newcube ] = cchengCEM_BCC(imagecube, num)

%% Code Comments:
% There were none
%

%% Authors Comments
%	Algorithm name: constrained band selection (CBS)
%	Authors: Chein-I Chang and Su Wang
%	Category: spectral filter
%	Designed criteria: least squares error
%	Designed method: linearly least squares constrained method
%	Typical use (LOI's addressed): BS
%	Inputs: reflectance or radiance cube
%	Outputs: gray-scale abundance fractional images
%	Assumptions: prior knowledge of the number of bands to be selected
%	Sensitivity to LOI (target knowledge): moderate
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation:
%	Brief description:
%		CBS was introduced in Chang and Wang (2006) to explore the idea of CEM that can be used for
%		BS. It interprets a band image as a desired target signature vector while considering other band
%		images as unknown signature vectors. As a result, the proposed CBS linearly constrains a band
%		image while also minimizing band correlation or dependence provided by other band images is
%		referred to as CEM-CBS. Four different criteria, referred to as band correlation constraint
%		(BCC), band correlation minimization (BCM), band dependence constraint (BDC), and band
%		dependence minimization (BDM), are derived for CEM-CBS. Since dimensionality resulting
%		from conversion of a band image to a vector may be huge, the CEM-CBS is further reinterpreted
%		as linearly constrained minimum variance (LCMV)-based CBS by constraining a band image as
%		a matrix, where the same four criteria BCM, BCC, BDC, and BDM can also be used for LCMVCBS.
%		In order to determine the number of bands required to select, p, a recently developed concept,
%		called virtual dimensionality (VD) is used to estimate the p. Once the p is determined, a set
%		of p desired bands can be selected by the CEM/LCMV-CBS. In what follows, only MATLAB
%		codes of four versions of CEM-CBS are provided, but these codes can be easily modified for
%		their counterparts of LCMV-CBS.


%% Code

close all;
[ xx, yy, band_num ] = size(imagecube);
%%%% get band image correlation %%%%
test_image = reshape(imagecube, xx*yy, band_num);
R = test_image * test_image'/band_num;
tt = inv(R);
%%%% get prioritization score for each band %%%%
for i = 1: band_num
    endmember_matrix = reshape(squeeze(imagecube(:, :, i)), xx*yy, 1);
    W = tt * endmember_matrix * inv(endmember_matrix' * tt * endmember_matrix);
    for j = 1: band_num
        if (i ~= j)
            test = reshape(squeeze(imagecube(:, :, j)), xx*yy, 1);
            score(i, j) = test' * W;
        else
            score(i,j) = 1;
        end
    end
end
clear endmember_matrix;
clear W;
clear i j;
%%%% select small subset of original bands %%%%
weight = zeros(1, band_num);
for i = 1: band_num
    test = score(i,:);
    scalar = sum(test) - score(i,i);
    weight(i) = scalar;
end
%weight = abs(weight);
original = 1:band_num;
coefficient_integer = weight * 100000;
band_select = zeros(1, num);
i = 1;
while (i <= num)
    max_coe = max(coefficient_integer(:));
    index = find(coefficient_integer == max_coe);
    if (length(index) == 1)
        band_select(i) = original(index);
        i = i + 1;
    else
        j = 1;
        while (j <= length(index))&&(i <= num )
            band_select(i) = original(index(j));
            i = i + 1;
            j = j + 1;
        end
    end
    coefficient_integer(index) = -10^10;
end
band_sort = sort(band_select);
clear coefficient_integer;
clear max_coe
clear index;
clear i;
clear j;
%get new imagecube
newcube = zeros(xx,yy, num);
for i = 1: xx
    for j = 1: yy
        test = squeeze(imagecube(i,j,:));
        newcube(i,j,:) = test(band_sort);
    end
end