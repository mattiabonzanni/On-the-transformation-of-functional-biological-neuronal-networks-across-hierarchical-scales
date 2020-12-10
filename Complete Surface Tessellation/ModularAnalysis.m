function Results=ModularAnalysis(RangeROIs,Results)
%% Required codes
% 1. ModularAnalysis.m
    % 1.1 clusterNodes.m;
    % 1.2 createModularGraph.m;
    % 1.3 avg_clus_matrix.m;
        % 1.3.1 clustering_coef_matrix.m;
    % 1.4 randomize_matrix.m;
    % 1.5 spatialLattice_matrix.m;
    % 1.5 avg_path_matrix.m;
    % 1.6 modularity_louvain_und.m;
%% Inputs-Outputs
        %Inputs:
% 1. RangeROIs= number of Region of Interests (ROIs) used to tessellate the whole graph;
% 2. xcoord: random x-axis coordinates for each node [0;1];
% 3. ycoord: random y-axis coordinates for each node [0;1];
% 4. Matrix: sparse Global matrix. 
        %Outputs:
% 1.  numberROIs: number of even ROIs in which you want to divide the surface.This value will be rounded.
% 2.  Graph: the Modular graph;
% 3.  Matrix: the Modular matrix;
% 4.  xcentroid: x-axis coordinates of each module;
% 5.  ycentroid: y-axis coordinates of each module;
% 6.  numberModules: number of modules in which more than one node is identified;
% 7.  ModulesComposition: structure containing the IDs of each node composing the module in the first column and the number of nodes in each module in the second column;
% 8.  AverageNodePerModule: average number of nodes per module;
% 9.  NetCC: average Clustering Coef. value of the Modular graph;
% 10. NetPL: average Path length value of the Modular graph;
% 11. RandomMatrix: the regular matrix of the Modular matrix;
% 12. RegulalrMatrix: the regular matrix of the Modular matrix;
% 13. RandCC: average Clustering Coef. value of the Random model of the Modular graph;
% 14. RegularCC: average Clustering Coef. value of the Lattice model of the Modular graph;
% 15. RandPL: average Path length value of the Random model of the Modular graph; 
% 16. RegularPL: average Path length value of the Lattice model of the Modular graph;
% 17. SWP: Small-world propensity of the Modular graph;
% 18. Q: modularity of the Modular graph;
% 19. SumDegree: the sum of the weighted edges of the Modular matrix;
% 20. SumMST: 
% 21. secondsmallest: second smallest eigenvalue of the Modular matrix;
% 22. largest: largest eigenvalue of the Modular Matrix;
% 23. syncMod: synchronizabilty of the Modular graph.

%% References codes

% 1) avg_clus_matrix (written by Eric Bridgeford);
% 2) avg_path_matrix (written by Eric Bridgeford);
% 3) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 4) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);
% 5) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);

% written by Mattia Bonzanni


for i=1:length(RangeROIs)
    numberROIs=RangeROIs(i);
%% Create Modules
    [module,numberModules,xcentroid,ycentroid]=clusterNodes(numberROIs,Results.Input.xcoord,Results.Input.ycoord);
%% Create the modular graph
    [ModularG,ModularMatrix]=createModularGraph(Results.Output.GlobalGraph.Matrix,numberModules,module);
    Results.Output.ModularGraph(i).numberROIs=numberROIs;
    Results.Output.ModularGraph(i).Graph=ModularG;
    Results.Output.ModularGraph(i).Matrix=ModularMatrix;
    Results.Output.ModularGraph(i).xcentroid=xcentroid;
    Results.Output.ModularGraph(i).ycentroid=ycentroid;
    Results.Output.ModularGraph(i).numberModules=numberModules;
    Results.Output.ModularGraph(i).ModulesComposition=module;
    Results.Output.ModularGraph(i).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC,PL and SWP
    Results.Output.ModularGraph(i).NetCC=avg_clus_matrix(ModularMatrix,'O');                                        % to compute NetCC
    Results.Output.ModularGraph(i).NetPL=avg_path_matrix(1./ModularMatrix);                                         % to compute NetPL
    Results.Output.ModularGraph(i).RandomMatrix=randomize_matrix(Results.Output.ModularGraph(i).Matrix);            % to compute Random matrix
    Results.Output.ModularGraph(i).RegularMatrix=spatialLattice_matrix(Results.Output.ModularGraph(i).numberModules,Results.Output.ModularGraph(i).Matrix,xcentroid,ycentroid); % to compute Lattice Matrix based on spatial arguments
    Results.Output.ModularGraph(i).RandCC= avg_clus_matrix(Results.Output.ModularGraph(i).RandomMatrix, 'O');       % to compute RandCC
    Results.Output.ModularGraph(i).RegularCC= avg_clus_matrix(Results.Output.ModularGraph(i).RegularMatrix, 'O');   % to compute RegularCC
    Results.Output.ModularGraph(i).RandPL=avg_path_matrix(1./Results.Output.ModularGraph(i).RandomMatrix);          % to compute RandPL
    Results.Output.ModularGraph(i).RegularPL=avg_path_matrix(1./Results.Output.ModularGraph(i).RegularMatrix);      % to compute RegularPL
    z = (Results.Output.ModularGraph(i).NetPL - Results.Output.ModularGraph(i).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (Results.Output.ModularGraph(i).RegularPL - Results.Output.ModularGraph(i).RandPL);
    if Results.Output.ModularGraph(i).NetPL == Inf || Results.Output.ModularGraph(i).RandPL == Inf || Results.Output.ModularGraph(i).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (Results.Output.ModularGraph(i).RegularCC - Results.Output.ModularGraph(i).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (Results.Output.ModularGraph(i).RegularCC - Results.Output.ModularGraph(i).RandCC);
    if isnan(Results.Output.ModularGraph(i).RegularCC) || isnan(Results.Output.ModularGraph(i).RandCC) || isnan(Results.Output.ModularGraph(i).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    Results.Output.ModularGraph(i).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
%% Modularity
    Results.Output.ModularGraph(i).Q=modularity_louvain_und(Results.Output.ModularGraph(i).Matrix);
%% Sum Degree    
    Results.Output.ModularGraph(i).SumDegree=sum(Results.Output.ModularGraph(i).Matrix,'all');
%% Sum of the MST    
    H=Results.Output.ModularGraph(i).Graph;                             % to assign the modular graph to variable H
    H.Edges.Weight=1./ H.Edges.Weight;                                  % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T1= minspantree(H);                                                 % to compute the MST
    edge=T1.Edges.EndNodes;                                             % to store the edges of the spanning tree
    for p=1:length(edge)
        s=edge(p,1);
        t=edge(p,2);
        position(p)=findedge(Results.Output.ModularGraph(i).Graph,s,t); % to find the edge values in the original Modular graph
        weight=Results.Output.ModularGraph(i).Graph.Edges(position,2);  % to store the edge values
    end
    Results.Output.ModularGraph(i).sumMST=sum(table2array(weight));     % to compute the sum of the edges composing the MST
%% Sync
    A=Results.Output.ModularGraph(i).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                                  % convert it into a diagonal matrix D
    L = D - A;                                                          % Calculate Laplacian L
    e = eig(L);                                                         %Compute the eigenspectrum of L(t)
    e = sort(e);                                                        % Sort smallest to largest eigenvalue
    Results.Output.ModularGraph(i).secondsmallest=e(2);                 % to store second smallest eigenvalue
    Results.Output.ModularGraph(i).largest=e(end);                      % to store largest eigenvalue
    Results.Output.ModularGraph(i).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty 
end
end