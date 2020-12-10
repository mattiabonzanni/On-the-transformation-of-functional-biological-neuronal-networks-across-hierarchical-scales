function [module,numberModules,xcentroid,ycentroid]=clusterNodes(numberROIs,xcoord,ycoord)
%% To compute the module compositions 
% Input:
% 1. numberROIs= number of even ROIs in which you want to divide the surface.This value will be rounded.
% 2. xcoord= x-axis node coordinates
% 3. ycoord= y-axis node coordinates

% Output:
% 1. module= structure containing the IDs of each node composing the module in the
%            first column and the number of nodes in each module in the 
%            second column
% 2. numberModules= number of ROI in which more than one node is identified
% 3. xcentroid= x-axis module coordinates
% 4. ycentroid= y-axis module coordinates

% written by Mattia Bonzanni 

%% Divide the space in even surfaces 
side=round(sqrt(numberROIs));                   % to calculate the side of each square based on the number of desired ROI; if the sqrt is not an integer, the number is rounded
spaceDivision=linspace(0,1,side+1);             % to equally divide each axis 
fprintf('%d modules will be computed\n',side^2)
%% Identify the modules and plot the the polygons of the ROI with at least 1 node
numberModules=0;
for i=1:side
    xvalues=xcoord;                                                                     % to load at each iteration the original coordinates and not the logic vectors derived from following steps
    spacecoord1=spaceDivision(i);       
    spacecoord2=spaceDivision(i+1);
    xvalues(xvalues<=spacecoord1)=0;
    xvalues(xvalues>spacecoord2)=0;
    nodex=find(xvalues);                                                                % to find the nodes in the given x-range             
    for j=1:side
        yvalues=ycoord;
        spacecoord3=spaceDivision(j);
        spacecoord4=spaceDivision(j+1);
        yvalues(yvalues<=spacecoord3)=0;
        yvalues(yvalues>spacecoord4)=0;
        nodey=find(yvalues);                
        findnode=intersect(nodex,nodey);                                                    % to identify the nodes that meet the geometrical criteria
        if length(findnode)>0                                                               % to apply the following lines if and only if the intersection is non-empty
           numberModules=numberModules+1;                                                   % to update the number of modules with more than one node                                       
           module{numberModules,1}=findnode;                                                % to store the module composition, namely which nodes belong to which module
           module{numberModules,2}=length(findnode);                                        % to store how many nodes in each module
        end     
    end
end
%% To calculate the coordinates of each module
for i=1:numberModules
    xcentroid(i)=median(xcoord(module{i,1}));                                               % median of the x-axis coordinates of the ith module
    ycentroid(i)=median(ycoord(module{i,1}));                                               % median of the y-axis coordinates of the ith module
end
fprintf('Modules computation completed\n')
end