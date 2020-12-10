function [moduleFront,numberModulesFront,xcentroidFront,ycentroidFront,zcentroidFront]=clusterNodesFRONT(numberROIs,xcoord,ycoord,zcoord,NG,p)
%% To compute the module compositions and characteristics 
% Front Graph --> z=[0;0.2]
% Input:
% 1. numberROIs= number of even ROIs in which you want to divide the
%                space.This value will be rounded.
% 2. xcoord= x-axis node coordinates
% 3. ycoord= y-axis node coordinates
% 4. zcoord= z-axis node coordinates
% 5. NG= number of nodes of the global graph
% 6. p= penetration fraction of the recording

% Output:
% 1. module= structure containing the composition of each module in the
%            first column and the number of nodes in each module in the 
%            second column
% 2. numberModules= number of ROI in which more than one node is identified
% 3. xcentroid= x-axis module coordinates
% 4. ycentroid= y-axis module coordinates
% 5. zcentroid= z-axis module coordinates
% 6. HiddenNodes= list of nodes outside the modular region

%% To create the scatter plot of the nodes
% scatter(xcoord,ycoord,'filled');
% hold on
%% Divide the space in even surfaces and plot them
side=round(sqrt(numberROIs));       % to calculate the side of each square based on the number of desired ROI; if the sqrt is not an integer, the number is rounded
spaceDivision=linspace(0,1,side+1); % to divide each axis in equal spaces
fprintf('%d modules will be computed\n',side^2)
%% Identify the modules and plot the the polygons of the ROI with at least 1 node
HiddenNodes=[1:NG];
numberModulesFront=0;
nodez=find(zcoord<p);
for i=1:side
    yvalues=ycoord;                                     % to load at each iteration the original coordinates and not the logic vectors
    spacecoord1=spaceDivision(i);       
    spacecoord2=spaceDivision(i+1);
    yvalues(yvalues<=spacecoord1)=0;
    yvalues(yvalues>spacecoord2)=0;
    nodey=find(yvalues);
    nodeyz=intersect(nodey, nodez);
    for j=1:side
        xvalues=xcoord;
        spacecoord3=spaceDivision(j);
        spacecoord4=spaceDivision(j+1);
        xvalues(xvalues<=spacecoord3)=0;
        xvalues(xvalues>spacecoord4)=0;
        nodex=find(xvalues);                
        findnode=intersect(nodeyz,nodex);                % to identify the nodes that meet the geometrical criteria
        if length(findnode)>0                           % to apply the following lines if and only if the intersection led to a non-empty result
           HiddenNodes(findnode)=0;                                 % to eliminate the nodes identified in the findnode line to reveal the hidden nodes
           numberModulesFront=numberModulesFront+1;               % to update the number of modules with more than one node
           moduleFront{numberModulesFront,1}=findnode;            % to store the module composition, namely which nodes belong to which module
           moduleFront{numberModulesFront,2}=length(findnode);    % to store how many nodes in each module
        end     
    end
end
for i=1:numberModulesFront
    xcentroidFront(i)=median(xcoord(moduleFront{i,1}));           % median of the x-axis coordinates of the ith module
    ycentroidFront(i)=median(ycoord(moduleFront{i,1}));           % median of the y-axis coordinates of the ith module
    zcentroidFront(i)=median(zcoord(moduleFront{i,1}));
end
HiddenNodes=find(HiddenNodes);
fprintf('Modules FRONT computation completed\n')
end