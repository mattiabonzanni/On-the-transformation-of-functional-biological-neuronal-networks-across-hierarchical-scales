function [moduleTop,numberModulesTop,xcentroidTop,ycentroidTop,zcentroidTop]=clusterNodes(numberROIs,xcoord,ycoord)
%% To compute the module compositions and characteristics 
% Input:
% 1. numberROIs= number of even ROIs in which you want to divide the
%                space.This value will be rounded.
% 2. xcoord= x-axis node coordinates
% 3. ycoord= y-axis node coordinates

% Output:
% 1. module= structure containing the composition of each module in the
%            first column and the number of nodes in each module in the 
%            second column
% 2. numberModules= number of ROI in which more than one node is identified
% 3. xcentroid= x-axis module coordinates
% 4. ycentroid= y-axis module coordinates

%% To create the scatter plot of the nodes
% scatter(xcoord,ycoord,'filled');
% hold on
%% Divide the space in even surfaces and plot them
side=round(sqrt(numberROIs));       % to calculate the side of each square based on the number of desired ROI; if the sqrt is not an integer, the number is rounded
spaceDivision=linspace(0,1,side+1); % to divide each axis in equal spaces
fprintf('%d modules will be computed\n',side^2)
%% Identify the modules and plot the the polygons of the ROI with at least 1 node
numberModulesTop=0;
nodey=find(ycoord>0.8);
for i=1:side
    xvalues=xcoord;                                     % to load at each iteration the original coordinates and not the logic vectors
    spacecoord1=spaceDivision(i);       
    spacecoord2=spaceDivision(i+1);
    xvalues(xvalues<=spacecoord1)=0;
    xvalues(xvalues>spacecoord2)=0;
    nodex=find(xvalues);
    nodexy=intersect(nodex, nodey);
    for j=1:side
        zvalues=zcoord;
        spacecoord3=spaceDivision(j);
        spacecoord4=spaceDivision(j+1);
        zvalues(zvalues<=spacecoord3)=0;
        zvalues(zvalues>spacecoord4)=0;
        nodez=find(zvalues);                
        findnode=intersect(nodexy,nodez);                % to identify the nodes that meet the geometrical criteria
        if length(findnode)>0                           % to apply the following lines if and only if the intersection led to a non-empty result
           numberModulesTop=numberModulesTop+1;               % to update the number of modules with more than one node
           moduleTop{numberModulesTop,1}=findnode;            % to store the module composition, namely which nodes belong to which module
           moduleTop{numberModulesTop,2}=length(findnode);    % to store how many nodes in each module
        end     
    end
end
for i=1:numberModulesTop
    xcentroidTop(i)=median(xcoord(moduleTop{i,1}));           % median of the x-axis coordinates of the ith module
    ycentroidTop(i)=median(ycoord(moduleTop{i,1}));           % median of the y-axis coordinates of the ith module
    zcentroidTop(i)=median(zcoord(moduleTop{i,1}));
end
fprintf('Modules computation completed\n')
end