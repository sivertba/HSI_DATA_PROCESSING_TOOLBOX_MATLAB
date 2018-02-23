function [FinalPositions, running_time]=cchengFIPPIoptimized(Image, InitialSkewers)%% Code Comments:
%% Code Comments:
% The Matlab PPI algorithm
% Fast Iterative Pixel Purity Index Algorithm
%
% Input parameters:
% ———————————————
% Image: Hyperspectral image data after MNF dimensionality reduction
% InitialSkewers: Positions of ATGP-generated pixels
%
% Output parameter:
% ———————————————
% FinalPositions: Positions of FIPPI-generated endmember pixels
% running_time - The total running time used by this run
%
% Authors: Chein-I Chang and Antonio Plaza
% Minor Modified by Chao-Cheng Wu

%% Authors Comments
%	Algorithm name: PPI
%	Authors: J. Boardman
%	Category: convex geometry
%	Designed criteria: orthogonal projection
%	Designed method: randomly generated vectors, called skewers
%	Typical use (LOI's addressed): endmember extraction with number of endmembers to be known
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

% Check CPU time at the beginning
start=cputime;
% Code initialization for data and visualization
[ns,nl,~]=size(Image);
[~,~]=size(InitialSkewers);
Extrema=zeros(ns,nl);
ProjectionScores=zeros(ns,nl);
subplot(2,1,1);
imagesc(Image(:,:,1)); colormap(gray);
title('Pixels extracted by FPPI:');
set(gca,'DefaultTextColor','black','xtick',[],...
    'ytick',[],'dataaspectratio',[1 1 1]);
po1 = get(gca,'position');
% Use ATGP-generated pixels as the initial skewers
NewSkewers=InitialSkewers;
% Begin iterative process
Other = 1;
SkewersUsed = [];
while (Other >= 1)
    [ne,~]=size(NewSkewers);
    disp(['Iteration: ' int2str(Other)]);
    disp(['Skewers: ' int2str(ne)]);
    for k = 1:ne
        [ne_old,~]=size(SkewersUsed);
        skewer=squeeze(Image(NewSkewers(k,1),NewSkewers(k,2),:));
        skewer=skewer/norm(skewer);
        SkewersUsed = union(SkewersUsed,skewer);
        [ne_new,~]=size(SkewersUsed);
        subplot(2,1,2);
        drawnow;
        plot(skewer);
        title(['Current skewer: ' int2str(k)]);
        if (ne_new~=ne_old)
            % Project all the sample data vectors onto this particular skewer
            for i=1:ns
                for j=1:nl
                    pixel = squeeze(Image(i,j,:));
                    ProjectionScores(i,j) = dot(skewer,pixel);
                end
            end
            % Obtain the extrema set for each skewer (maximum and minimum
            % projection)
            [~,mpos] = max(ProjectionScores(:));
            [~,pos] = min(ProjectionScores(:));
            mposx = floor((mpos-1)/ns)+1; mposy = mod(mpos-1,ns)+1;
            posx = floor((pos-1)/ns)+1; posy = mod(pos-1,ns)+1;
            % Display the pixel positions of the pixels in the extrema set
            drawnow;
            subplot(2,1,1);
            text(mposx,mposy,'o','Margin',1,...
                'HorizontalAlignment','center',...
                'FontSize',22,'FontWeight','light',...
                'FontName', 'Garamond','Color','yellow');
            drawnow;
            subplot(2,1,1);
            text(posx,posy,'o','Margin',1,...
                'HorizontalAlignment','center',...
                'FontSize',22,'FontWeight','light',...
                'FontName','Garamond','Color','yellow');
            % Increase PPI count of extrema pixels
            Extrema(posy,posx)=Extrema(posy,posx)+1;
            Extrema(mposy,mposx)=Extrema(mposy,mposx)+1;
            % Incorporate sample vectors with PPI count greater than zero to
            % the skewer set
            vnew = [ mposx mposy ; posx posy ];
            NewSkewers = union(NewSkewers,vnew,'rows');
        end
    end
    % Check stopping rule
    [ne2,~]=size(NewSkewers);
    if (ne2==ne)
        Other = 0;
    else
        Other = Other+1;
        % Extrema=zeros(ns,nl);
    end
end
% Produce the positions of the final endmember set
Binary = Extrema>0;
ne = sum(Binary(:));
disp(['Extracted endmembers: ' int2str(ne)]);
FinalPositions = zeros(ne,2);
Current = 1;
for i=1:ns
    for j=1:nl
        if Binary(i,j)>0
            FinalPositions(Current,1)=i;
            FinalPositions(Current,2)=j;
            Current = Current+1;
        end
    end
end
% Check CPU time at the end
stop=cputime;
running_time=stop-start;
