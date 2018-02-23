function [eeindex, score, duration]=cchengPPI(imagecub,skewer_no)
%% Code Comments:
% The Matlab PPI algorithm
% ————— Input variables ———————————
% 'imagecub' - The hyperspectral image cube
% 'skewer_no' - The number of skewers
%
% ————— Output variables —————————
% 'eeindex' - The locations of the final endmembers (x,y)
% 'score' - The PPI score of each pixel
% 'duration - The number of seconds used to run this program

%% Authors Comments
%	Algorithm name: PPI
%	Authors: J. Boardman
%	Category: convex geometry
%	Designed criteria: orthogonal projection
%	Designed method: randomly generated vectors, called skewers
%	Typical use (LOI’s addressed): endmember extraction with number of endmembers to be known
%	Inputs: reflectance or radiance cube
%	Outputs: endmembers
%	Assumptions: the number of skewers to be generated, K must be sufficiently large
%	Sensitivity to LOI (target knowledge): high to K as well as skewers
%	Sensitivity to noise: moderate
%	Operating bands: VNIR through LWIR
%	Maturity: mature/operational
%	Effectiveness: high
%	Implementation: simple, easy to use, but not real time
%	Brief description:
%		PPI was first developed by Boardman to extract endmembers (Boardman, 1994) and is available
%		in ENVI software. Unfortunately, its detailed steps in implementation are not available for users
%		who would like to make modifications or changes at their discretion. Its MATLAB codes provided
%		in the following serve this purpose. It is developed based on the concept of the convex
%		geometry and the criterion of orthogonal projection. It first generates a set of K random unit
%		vectors, called skewers, to cover all possible projection directors and then orthogonally projects
%		all data sample vectors on these skewers to find the maximal and minimal orthogonal projections
%		of each skewer. For each data sample vector, it counts the number of skewers on which its
%		orthogonal projections yield either maximal or minimal projections. This count is referred to as
%		the PPI count, which will be used to determine whether or not a particular data sample vector is
%		an endmember. In doing so, a threshold needs to be specified in advance. In ENVI, this is done
%		manually. In addition, since its skewers are generated randomly, the results are not repeatable. In
%		other words, the same user running PPI in different times or different users running PPI at
%		the same time will all have different results. This serious drawback can further be fixed by
%		initialization-driven PPI (ID-PPI), such as a fast iterative PPI (FIPPI) algorithm (Chang and
%		Plaza, 2006), and random PPI (RPPI) (Chang et al., 2010).


%% Code
% Initial Variables
[rows, columns, bands]=size(imagecub);
score=zeros(rows*columns,1);
switch_results=1;
% Record the start CPU time
start=cputime();
% Separate the total number of skewers into several sets a
%nd each set uses 500 skewers

skewer_sets=floor(skewer_no/500)+1;
last_skewer_no=mod(skewer_no,500);
for i=1:skewer_sets
    if (skewer_sets-i) == 0
        skewer_no=last_skewer_no;
    else
        skewer_no=500;
    end
    % Generate skewers
    rand('state',sum(100*clock));
    skewers=rand(bands,skewer_no)-0.5;
    % Normalize skewers
    for j=1:skewer_no
        skewers(:,j)=skewers(:,j)/norm(skewers(:,j));
    end
    % project every sample vector to the skewers
    projcub=reshape(imagecub, rows*columns, bands);
    proj_result=projcub*skewers;
    % Find the extrema set for each skewer and add 1 to their score
    for k=1:skewer_no
        max_pos=find(proj_result(:,k)==max(proj_result(:,k)));
        min_pos=find(proj_result(:,k)==min(proj_result(:,k)));
        score(max_pos)=score(max_pos)+1;
        score(min_pos)=score(min_pos)+1;
    end
end
% Find the pixel which has score larger than 0
result=find(score>0);
% Find the position of the p highest scores
%result=[];
%for i=1:skewer_no,
% result=[result find(max(score)==score,1)];
% score(find(max(score)==score,1))=0;
%end
% Convert one dimension to two dimension index
if(switch_results)
    eeindex=translate_index(result,rows,columns,1);
else
    if(mod(result,rows)==0)
        eeindex(2,:)=floor(result./rows);
    else
        eeindex(2,:)=floor(result./rows)+1;
    end
    eeindex(1,:)=mod(result-1,rows)+1;
end
duration=cputime()-start;