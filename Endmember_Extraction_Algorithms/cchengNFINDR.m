function [endmemberindex, duration]=cchengNFINDR(imagecube,p)
%% Code Comments:
% The N-FINDR algorithm
%
% ————— Input variables ———————————
% 'imagecube' - The data transformed components [row column band]
% 'p' - The number of endmembers to be generated
%
% if band > p, then the program will
% automatically use Singular Value Decomposition to calculate the volume
%
% ————— Output variables —————————
% 'endmemberindex - The locations of the final endmembers (x,y)
% 'duration - The number of seconds used to run this program
%% Authors Comments
%	Algorithm name: N-finder algorithm (N-FINDR)
%	Authors: M.E. Winter
%	Category: convex geometry
%	Designed criteria: maximum simplex volume
%	Designed method: finding a simplex with maximum volume
%	Typical use (LOI's addressed): endmember extraction with number of endmembers to be known
%	Inputs: reflectance or radiance cube
%	Outputs: endmembers
%	Assumptions: All the vertices of a simplex with maximal volume should be specified by endmembers
%	Sensitivity to LOI (target knowledge): high to value of p
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		The N-FINDR (Winter, 1999a, 1999b, 2004) is a popular endmember extraction algorithm
%		(EEA) other than PPI. It is quite different from PPI in terms of its design criterion. The NFINDR
%		makes use of maximum simplex volume as a criterion as opposed to orthogonal projection
%		used by PPI. So, its design criteria make PPI an unconstrained EEA and N-FINDR a fully
%		1020 Appendix: Algorithm Compendium
%		constrained EEA. Like PPI, which requires prior knowledge of the number of skewers, K, NFINDR
%		also needs to know the number of endmembers, p, a priori. Similar to PPI, N-FINDR
%		also suffers the same issue encountered in PPI, which is its use of random initial conditions that
%		result in unrepeatable and inconsistent endmember results. This issue is also addressed by
%		Chang et al. (Plaza and Chang, 2006; Chang et al., 2011b). Another serious issue arising in NFINDR
%		implementation that is not encountered in PPI is its very high computational complexity.
%		Many research efforts have been reported to address this issue. Details can be found in Xiong
%		et al. (2011) and Chang (2013). Also, real-time implementation of N-FINDR has also been
%		proposed in Wu et al. (2010) with details in Chang (2013).

%% Code
% Set initial condition
endmemberindex=[];
newvolume = 0;
prevolume = -1;
[row, column, band]=size(imagecube);
switch_results=1;
% Determine to use SVD to calculate the volume or not
if(band > p)
    use_svd=1;
else
    use_svd=0;
end
% Start to count the CPU computing time
start=cputime();
% Randomly select p initial endmembers
rand('state',sum(100*clock));
for i=1:p
    while(1)
        temp1=round(row*rand);
        temp2=round(column*rand);
        if(temp1>0 && temp2>0)
            break;
        end
    end
    endmemberindex=[endmemberindex;[temp1 temp2]];
end
endmemberindex=endmemberindex';
% Generate endmember vector from reduced cub
display(endmemberindex);
endmember=[];
for i=1:p
    if(use_svd)
        endmember=[endmember squeeze(...
            imagecube(endmemberindex(1,i),endmemberindex(2,i),:))];
    else
        endmember=[endmember squeeze(...
            imagecube(endmemberindex(1,i), endmemberindex(2,i),1:p-1))];
    end
end
% calculate the endmember's volume
if(use_svd)
    s=svd(endmember);
    endmembervolume=1;
    for i=1:p
        endmembervolume=endmembervolume*s(i);
    end
else
    jointmatrix=[ones(1,p) ; endmember];
    endmembervolume=abs(det(jointmatrix))/factorial(p-1);
end
% The main algorithm
while newvolume > prevolume % if the new generated endmember volume is larger than the old one, continue the algorithm
    % Use each sample vector to replace the original one, and calculate new volume
    for i=1:row
        for j=1:column
            for k=1:p
                caculate=endmember;
                if(use_svd)
                    caculate(:,k)=squeeze(imagecube(i, j, :));
                    s=svd(caculate);
                    volume=1;
                    for z=1:p
                        volume=volume*s(z);
                    end
                else
                    caculate(:,k)=squeeze(imagecube(i, j, 1:p-1));
                    jointmatrix=[ones(1,p);calculate];
                    volume=abs(det(jointmatrix))/factorial(p-1); % The formula of
                    Simplex volume
                end
                if volume > endmembervolume
                    endmemberindex(:,k)=[i;j];
                    endmember=calculate;
                    endmembervolume=volume;
                end
            end
        end
    end
    prevolume=newvolume;
    newvolume=endmembervolume;
end
stop=cputime();
duration=stop-start;
% Switch results for the standard
if(switch_results)
    endmemberindex(3,:)=endmemberindex(1,:);
    endmemberindex(1,:)=[];
    endmemberindex=endmemberindex';
end

