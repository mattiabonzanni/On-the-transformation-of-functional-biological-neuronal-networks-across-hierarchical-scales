function [ModularGraphRight]=RightAnalysis(RangeROIs,xcoord,ycoord,zcoord,NG,distMatrix,p)
% Right Graph --> x=[0.8;0.1], y=[0;1], z=[0;1]
% To isolate the right section, compute the modular graph (Visible) and isolate the remaining graph (Hidden)
% % Input:
% 1. RangeROIs= number of ROIs to test at each iteration
% 2. xcoord= x-axis node coordinates
% 3. ycoord= y-axis node coordinates
% 4. zcoord= z-axis coordinates
% 5. NG= number of nodes of the global graph
% 6. distMatrix= distance-dependent matrix (weighted,symmetric)
% 7. p= penetration fraction of the recording

% Output:
% 1. ModularGraphBottom= a structure variable containing:
% a. numberROIs= number of even ROIs in which you want to divide the space.This value will be rounded.
% b. numberModules= number of ROI in which more than one node is identified
% c. ModulesComposition= structure containing the composition of each module in the first column and the number of nodes in each module in the second column
% d. GraphVisible= modular graph
% e. MatrixVisible= modular matrix (weighted)
% f. xcentroid= x-axis module coordinates
% g. ycentroid= y-axis module coordinates
% h. zcentroid= z-axis module coordinates
% i. AverageNodePerModule= average number of nodes in each module
% j. NodeFractionPerModule= average number of nodes normalized by NG per module
% k. NumberModulePerNG= number of modules normalized by NG
% l. Surface Percentage= surface percentage per module
% m. NetPLmodular= average path lenght of the modular graph
% n. NetCCmodular= average clustering coef of the modular graph
% o. HiddenNodes= list of nodes outside the modular region
% p. HiddenMatrix= submatrix of the hidden nodes
% q. NetCChidden= average clustering coef of the hidden graph
% r. NetPLhidden= avergae path length of the hidden graph

for i=1:length(RangeROIs)
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid,HiddenNodes]=clusterNodesRIGHT(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    net_path_Modular = avg_path_matrix(1./ModularMatrix);            % average path of the network
    net_clus_Modular = full(avg_clus_matrix(ModularMatrix,'O'));     % average clustering coef. with the Onella method ('O')
    ModularGraphRight(i).numberROIs=numberROIs;
    ModularGraphRight(i).numberModules=numberModules;
    ModularGraphRight(i).ModulesComposition=module;
    ModularGraphRight(i).GraphVisible=ModularG;
    ModularGraphRight(i).MatrixVisible=ModularMatrix;
    ModularGraphRight(i).xcentroid=xcentroid;
    ModularGraphRight(i).ycentroid=ycentroid;
    ModularGraphRight(i).zcentroid=zcentroid;
    ModularGraphRight(i).AverageNodePerModule=mean(cell2mat(module(:,2)));
    ModularGraphRight(i).NodeFractionPerModule=(mean(cell2mat(module(:,2))))/NG;   % average number of nodes/modules divided by NG
    ModularGraphRight(i).NumberModulePerNG=ModularGraphRight(i).numberModules/NG; % number of modules divided by NG. 
    ModularGraphRight(i).SurfacePercentage=100/numberROIs;
    ModularGraphRight(i).NetCCmodular=net_clus_Modular;
    ModularGraphRight(i).NetPLmodular=net_path_Modular;
end
HiddenMatrix=distMatrix(HiddenNodes,HiddenNodes);
net_path_Hidden = avg_path_matrix(1./HiddenMatrix);             % average path of the hidden network
net_clus_Hidden = full(avg_clus_matrix(HiddenMatrix,'O'));      % average clustering coef. with the Onella method ('O')
for i=1:length(RangeROIs)                                       % since the hidden nodes are the same independently of the size of the modules, I compute them just once and then populate the structure variable 
    ModularGraphRight(i).HiddenNodes=HiddenNodes;
    ModularGraphRight(i).HiddenMatrix=HiddenMatrix;
    ModularGraphRight(i).NetCChidden=net_clus_Hidden;
    ModularGraphRight(i).NetPLhidden=net_path_Hidden;
end
end