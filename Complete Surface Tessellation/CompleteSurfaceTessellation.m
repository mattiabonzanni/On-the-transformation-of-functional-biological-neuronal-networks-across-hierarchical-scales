     %% Complete Surface Tessellation
% Aim1: To generate a 2D network (Global Network) and, based on spatial
% arguments,completelly renormalize it (Modular Network). The complete
% renormalization implies a tessellation of the entire network. 
% Aim2: To study the network properties (average clustering coefficient, average path length,
% small-world propensity, average edge weight, modularity,
% synchronizability and fraction of long-term connections) between the
% Global vs Modular networks.
%% Nomenclature:
% - Global (or Whole) graph: to mimic microscale neuronal networks;
% - Modular Graph: to mimic mesoscale neuronal networks (computed from whole graph);
% - Neighbor  nodes: nodes which are connected after the pruning step in the DistanceDependentMatrix.m script; they are spatially closed to each other
% - Region of Interests (ROIs) and modules are synonymous.  
% - Regular (or Lattice): regular null model (each node is attached to its neighbours);
% - Random: random null model (random wiring);
% - Net: the studied network.
%% Required codes
% 1. DistanceDependentMatrix.m
    % 1.1 avg_clus_matrix.m;
        % 1.1.1 clustering_coef_matrix.m;
    % 1.2 randomize_matrix.m
    % 1.3 avg_path_matrix.m;
    % 1.4 modularity_louvain_und.m;
% 2. ModularAnalysis.m
    % 2.1 clusterNodes.m;
    % 2.2 createModularGraph.m;
    % 2.3 avg_clus_matrix.m;
        % 2.3.1 clustering_coef_matrix.m;
    % 2.4 randomize_matrix.m;
    % 2.5 spatialLattice_matrix.m
    % 2.5 avg_path_matrix.m;
    % 2.6 modularity_louvain_und.m;
%% Inputs-Outputs
        %Inputs:
% 1. NG= total number of seeded nodes in the whole graph;
% 2. beta= probability that an edge is randomly re-wired to a non-neighbor  node;
% 3. RangeROIs= number of Region of Interests (ROIs) used to tessellate the whole graph;
% 4. edgeDensity= the average percentage of neighbor nodes in the whole graph; 
        %Outputs:
% A Results structure containing;
% Field1 - Output
    % GlobalGraph structure
    % 1.  xcoord: random x-axis coordinates for each node [0;1];
    % 2.  ycoord: random y-axis coordinates for each node [0;1];
    % 3.  LatticeMatrix: the regular matrix of the Global matrix computed in the DistanceDependentMatrix.m before the pruning step;
    % 4.  wt: calculated weight threshold to achieve the desired edge density; 
    % 5.  Graph: the Global graph;
    % 6.  Matrix: sparse Global matrix;
    % 7.  NetCC: average Clustering Coef. value of the Global graph;
    % 8.  NetPL: average Path length value of the Global graph;
    % 9.  RegCC: average Clustering Coef. value of the Lattice model of the Global graph;
    % 10. RegPL: average Path length value of the Lattice model of the Global graph;
    % 11. RandomMatrix: the regular matrix of the Global matrix;
    % 12. RandCC: average Clustering Coef. value of the Random model of the Global graph;
    % 13. RandPL: average Path length value of the Random model of the Global graph;
    % 14. SWP: Small-world propensity of the Global graph;
    % 15. Q: modularity of the Global graph;
    % 16. SumDegree: the sum of the weighted edges of the Global matrix;
    % 17. secondsmallest: second smallest eigenvalue of the Global matrix;
    % 18. largest: largest eigenvalue of the Global Matrix;
    % 19. sync: synchronizabilty of the Global graph.  
    
    % ModularGraph structure
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
% Field2 - Input
    % A structure variable containing NG, beta, RangeROIs, edgeDensity(all input by the user) and the generated random coordinates xcoord and ycoord (from DistanceDependentMatrix.m)
% Field3 - FileName
    % The character variable is the generated name for the file containing information about the number of nodes, the beta and the time. 
%% References codes
% 1) DistanceDependentMatrix (written by Mattia Bonzanni);
% 2) ModularAnalysis (written by Mattia Bonzanni);
% 3) avg_clus_matrix (written by Eric Bridgeford);
% 4) avg_path_matrix (written by Eric Bridgeford);
% 5) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 6) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);
% 7) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);

% Written by Mattia Bonzanni 

%% Inputs
clear
clc
NG=200;
beta=0.05;
RangeROIs=[4:12].^2;        % to test a range of ROIs 
edgeDensity=8;              % the percentage of NG nodes attached to a given node in the Global graph. If not specified, use a value in the range [8:15] (She et al., 2016);

%% Save the inputs in the Results
Results.Input.NG=NG;
Results.Input.RangeROIs=RangeROIs;
Results.Input.edgeDensity=edgeDensity;
Results.Input.betarewiring=beta;
%% Compute global and modular graphs
Results=DistanceDependentMatrix(NG,edgeDensity,beta);   % compute also Regular and Random matrices
Results=ModularAnalysis(RangeROIs,Results);
%% Save the output
time=datestr(now,'mmddyy_HHMM');
s=sprintf('NG%db%s_',NG,num2str(beta*1000));
filename=[s,time];
Results.FileName=filename;
save(filename,'Results')    % save the file with current date 
clearvars -except Results
disp('End of computation')

