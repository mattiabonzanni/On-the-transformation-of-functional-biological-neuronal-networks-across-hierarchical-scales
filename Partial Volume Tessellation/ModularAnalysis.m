function [ModularGraph]=ModularAnalysis(RangeROIs,xcoord,ycoord,zcoord,NG,distMatrix,p)
%% Graphs Orientation
% Bottom Graph -->  y=[0;p], x=[0;1], z=[0;1].
% Top Graph -->     y=[1-p;1], x=[0;1], z=[0;1].
% Left Graph -->    y=[0:1], x=[0;p], z=[0;1].
% Right Graph -->   y=[0:1], x=[1-p;1], z=[0;1].
% Front Graph -->   y=[0:1], x=[0;1], z=[0;p].
% Back Graph -->    y=[0:1], x=[0;1], z=[1-p;1].

% To isolate the bottom section, compute the modular graph (Visible) and isolate the remaining graph (Hidden)
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



n=length(RangeROIs);
%% Bottom
for i=1:n
    numberROIs=RangeROIs(i);
%% Create Bottom Modules
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesBOTTOM(numberROIs,xcoord,ycoord,zcoord,NG,p);
%% Create the Bottom modular graph    
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    ModularGraph(i).numberROIs=numberROIs;
    ModularGraph(i).Graph=ModularG;
    ModularGraph(i).Matrix=ModularMatrix;
    ModularGraph(i).xcentroid=xcentroid;
    ModularGraph(i).ycentroid=ycentroid;
    ModularGraph(i).zcentroid=zcentroid;    
    ModularGraph(i).numberModules=numberModules;
    ModularGraph(i).ModulesComposition=module;
    ModularGraph(i).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(i).NetCC=avg_clus_matrix(ModularMatrix,'O');                         % to compute NetCC
    ModularGraph(i).NetPL=avg_path_matrix(1./ModularMatrix);                          % to compute NetPL
    ModularGraph(i).RandomMatrix=randomize_matrix(ModularGraph(i).Matrix);            % to compute Random matrix
    ModularGraph(i).RegularMatrix=spatialLattice_matrix3d(ModularGraph(i).numberModules,ModularGraph(i).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(i).RandCC= avg_clus_matrix(ModularGraph(i).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(i).RegularCC= avg_clus_matrix(ModularGraph(i).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(i).RandPL=avg_path_matrix(1./ModularGraph(i).RandomMatrix);          % to compute RandPL
    ModularGraph(i).RegularPL=avg_path_matrix(1./ModularGraph(i).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(i).NetPL - ModularGraph(i).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(i).RegularPL - ModularGraph(i).RandPL);
    if ModularGraph(i).NetPL == Inf || ModularGraph(i).RandPL == Inf || ModularGraph(i).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(i).RegularCC - ModularGraph(i).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(i).RegularCC - ModularGraph(i).RandCC);
    if isnan(ModularGraph(i).RegularCC) || isnan(ModularGraph(i).RandCC) || isnan(ModularGraph(i).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(i).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(i).Q=modularity_louvain_und(ModularGraph(i).Matrix);
%% Sum Degree    
    ModularGraph(i).SumDegree=sum(ModularGraph(i).Matrix,'all');
%% Sum of the MST    
    H=ModularGraph(i).Graph;                                % to assign the modular graph to variable H
    H.Edges.Weight=1./ H.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T1= minspantree(H);                                     % to compute the MST
    edge=T1.Edges.EndNodes;                                 % to store the edges of the spanning tree
    for p1=1:length(edge)
        s=edge(p1,1);
        t=edge(p1,2);
        position(p1)=findedge(ModularGraph(i).Graph,s,t);    % to find the edge values in the original Modular graph
        weight=ModularGraph(i).Graph.Edges(position,2);    % to store the edge values
    end
    ModularGraph(i).sumMST=sum(table2array(weight));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(i).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(i).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(i).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(i).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty 
end
%% Top
for i=1:n
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesTOP(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    j=n+i;
    ModularGraph(j).numberROIs=numberROIs;
    ModularGraph(j).Graph=ModularG;
    ModularGraph(j).Matrix=ModularMatrix;
    ModularGraph(j).xcentroid=xcentroid;
    ModularGraph(j).ycentroid=ycentroid;
    ModularGraph(j).zcentroid=zcentroid;    
    ModularGraph(j).numberModules=numberModules;
    ModularGraph(j).ModulesComposition=module;
    ModularGraph(j).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(j).NetCC=avg_clus_matrix(ModularGraph(j).Matrix,'O');                % to compute NetCC
    ModularGraph(j).NetPL=avg_path_matrix(1./ModularGraph(j).Matrix);                 % to compute NetPL
    ModularGraph(j).RandomMatrix=randomize_matrix(ModularGraph(j).Matrix);            % to compute Random matrix
    ModularGraph(j).RegularMatrix=spatialLattice_matrix3d(ModularGraph(j).numberModules,ModularGraph(j).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(j).RandCC= avg_clus_matrix(ModularGraph(j).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(j).RegularCC= avg_clus_matrix(ModularGraph(j).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(j).RandPL=avg_path_matrix(1./ModularGraph(j).RandomMatrix);          % to compute RandPL
    ModularGraph(j).RegularPL=avg_path_matrix(1./ModularGraph(j).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(j).NetPL - ModularGraph(j).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(j).RegularPL - ModularGraph(j).RandPL);
    if ModularGraph(j).NetPL == Inf || ModularGraph(j).RandPL == Inf || ModularGraph(j).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(j).RegularCC - ModularGraph(j).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(j).RegularCC - ModularGraph(j).RandCC);
    if isnan(ModularGraph(j).RegularCC) || isnan(ModularGraph(j).RandCC) || isnan(ModularGraph(j).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(j).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(j).Q=modularity_louvain_und(ModularGraph(j).Matrix);
%% Sum Degree    
    ModularGraph(j).SumDegree=sum(ModularGraph(j).Matrix,'all');
%% Sum of the MST    
    H2=ModularGraph(j).Graph;                                % to assign the modular graph to variable H
    H2.Edges.Weight=1./ H2.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T12= minspantree(H2);                                     % to compute the MST
    edge2=T12.Edges.EndNodes;                                 % to store the edges of the spanning tree
    for p2=1:length(edge2)
        s2=edge2(p2,1);
        t2=edge2(p2,2);
        position2(p2)=findedge(ModularGraph(j).Graph,s2,t2);    % to find the edge values in the original Modular graph
        weight2=ModularGraph(j).Graph.Edges(position2,2);    % to store the edge values
    end
    ModularGraph(j).sumMST=sum(table2array(weight2));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(j).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(j).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(j).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(j).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty 
end
%% Left
for i=1:n
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesLEFT(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    j=2*n+i;
    ModularGraph(j).numberROIs=numberROIs;
    ModularGraph(j).Graph=ModularG;
    ModularGraph(j).Matrix=ModularMatrix;
    ModularGraph(j).xcentroid=xcentroid;
    ModularGraph(j).ycentroid=ycentroid;
    ModularGraph(j).zcentroid=zcentroid;    
    ModularGraph(j).numberModules=numberModules;
    ModularGraph(j).ModulesComposition=module;
    ModularGraph(j).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(j).NetCC=avg_clus_matrix(ModularMatrix,'O');                         % to compute NetCC
    ModularGraph(j).NetPL=avg_path_matrix(1./ModularMatrix);                          % to compute NetPL
    ModularGraph(j).RandomMatrix=randomize_matrix(ModularGraph(j).Matrix);            % to compute Random matrix
    ModularGraph(j).RegularMatrix=spatialLattice_matrix3d(ModularGraph(j).numberModules,ModularGraph(j).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(j).RandCC= avg_clus_matrix(ModularGraph(j).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(j).RegularCC= avg_clus_matrix(ModularGraph(j).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(j).RandPL=avg_path_matrix(1./ModularGraph(j).RandomMatrix);          % to compute RandPL
    ModularGraph(j).RegularPL=avg_path_matrix(1./ModularGraph(j).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(j).NetPL - ModularGraph(j).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(j).RegularPL - ModularGraph(j).RandPL);
    if ModularGraph(j).NetPL == Inf || ModularGraph(j).RandPL == Inf || ModularGraph(j).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(j).RegularCC - ModularGraph(j).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(j).RegularCC - ModularGraph(j).RandCC);
    if isnan(ModularGraph(j).RegularCC) || isnan(ModularGraph(j).RandCC) || isnan(ModularGraph(j).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(j).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(j).Q=modularity_louvain_und(ModularGraph(j).Matrix);
%% Sum Degree    
    ModularGraph(j).SumDegree=sum(ModularGraph(j).Matrix,'all');
%% Sum of the MST    
    H3=ModularGraph(j).Graph;                                % to assign the modular graph to variable H
    H3.Edges.Weight=1./ H3.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T13= minspantree(H3);                                     % to compute the MST
    edge3=T13.Edges.EndNodes;                                 % to store the edges of the spanning tree
    for p3=1:length(edge3)
        s3=edge3(p3,1);
        t3=edge3(p3,2);
        position3(p3)=findedge(ModularGraph(j).Graph,s3,t3);    % to find the edge values in the original Modular graph
        weight3=ModularGraph(j).Graph.Edges(position3,2);    % to store the edge values
    end
    ModularGraph(j).sumMST=sum(table2array(weight3));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(j).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(j).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(j).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(j).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty end
end
%% Rigth
for i=1:n
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesRIGHT(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    j=3*n+i;
    ModularGraph(j).numberROIs=numberROIs;
    ModularGraph(j).Graph=ModularG;
    ModularGraph(j).Matrix=ModularMatrix;
    ModularGraph(j).xcentroid=xcentroid;
    ModularGraph(j).ycentroid=ycentroid;
    ModularGraph(j).zcentroid=zcentroid;    
    ModularGraph(j).numberModules=numberModules;
    ModularGraph(j).ModulesComposition=module;
    ModularGraph(j).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(j).NetCC=avg_clus_matrix(ModularMatrix,'O');                         % to compute NetCC
    ModularGraph(j).NetPL=avg_path_matrix(1./ModularMatrix);                          % to compute NetPL
    ModularGraph(j).RandomMatrix=randomize_matrix(ModularGraph(j).Matrix);            % to compute Random matrix
    ModularGraph(j).RegularMatrix=spatialLattice_matrix3d(ModularGraph(j).numberModules,ModularGraph(j).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(j).RandCC= avg_clus_matrix(ModularGraph(j).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(j).RegularCC= avg_clus_matrix(ModularGraph(j).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(j).RandPL=avg_path_matrix(1./ModularGraph(j).RandomMatrix);          % to compute RandPL
    ModularGraph(j).RegularPL=avg_path_matrix(1./ModularGraph(j).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(j).NetPL - ModularGraph(j).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(j).RegularPL - ModularGraph(j).RandPL);
    if ModularGraph(j).NetPL == Inf || ModularGraph(j).RandPL == Inf || ModularGraph(j).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(j).RegularCC - ModularGraph(j).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(j).RegularCC - ModularGraph(j).RandCC);
    if isnan(ModularGraph(j).RegularCC) || isnan(ModularGraph(j).RandCC) || isnan(ModularGraph(j).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(j).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(j).Q=modularity_louvain_und(ModularGraph(j).Matrix);
%% Sum Degree    
    ModularGraph(j).SumDegree=sum(ModularGraph(j).Matrix,'all');
%% Sum of the MST    
    H4=ModularGraph(j).Graph;                                % to assign the modular graph to variable H
    H4.Edges.Weight=1./ H4.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T14= minspantree(H4);                                     % to compute the MST
    edge4=T14.Edges.EndNodes                                 % to store the edges of the spanning tree
    for p4=1:length(edge4)
        s4=edge4(p4,1)
        t4=edge4(p4,2)
        position4(p4)=findedge(ModularGraph(j).Graph,s4,t4)    % to find the edge values in the original Modular graph
        weight4=ModularGraph(j).Graph.Edges(position4,2);    % to store the edge values
    end
    ModularGraph(j).sumMST=sum(table2array(weight4));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(j).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(j).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(j).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(j).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty 
end
%% Front
for i=1:n
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesFRONT(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    j=4*n+i;
    ModularGraph(j).numberROIs=numberROIs;
    ModularGraph(j).Graph=ModularG;
    ModularGraph(j).Matrix=ModularMatrix;
    ModularGraph(j).xcentroid=xcentroid;
    ModularGraph(j).ycentroid=ycentroid;
    ModularGraph(j).zcentroid=zcentroid;    
    ModularGraph(j).numberModules=numberModules;
    ModularGraph(j).ModulesComposition=module;
    ModularGraph(j).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(j).NetCC=avg_clus_matrix(ModularMatrix,'O');                         % to compute NetCC
    ModularGraph(j).NetPL=avg_path_matrix(1./ModularMatrix);                          % to compute NetPL
    ModularGraph(j).RandomMatrix=randomize_matrix(ModularGraph(j).Matrix);            % to compute Random matrix
    ModularGraph(j).RegularMatrix=spatialLattice_matrix3d(ModularGraph(j).numberModules,ModularGraph(j).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(j).RandCC= avg_clus_matrix(ModularGraph(j).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(j).RegularCC= avg_clus_matrix(ModularGraph(j).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(j).RandPL=avg_path_matrix(1./ModularGraph(j).RandomMatrix);          % to compute RandPL
    ModularGraph(j).RegularPL=avg_path_matrix(1./ModularGraph(j).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(j).NetPL - ModularGraph(j).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(j).RegularPL - ModularGraph(j).RandPL);
    if ModularGraph(j).NetPL == Inf || ModularGraph(j).RandPL == Inf || ModularGraph(j).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(j).RegularCC - ModularGraph(j).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(j).RegularCC - ModularGraph(j).RandCC);
    if isnan(ModularGraph(j).RegularCC) || isnan(ModularGraph(j).RandCC) || isnan(ModularGraph(j).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(j).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(j).Q=modularity_louvain_und(ModularGraph(j).Matrix);
%% Sum Degree    
    ModularGraph(j).SumDegree=sum(ModularGraph(j).Matrix,'all');
%% Sum of the MST    
    H5=ModularGraph(j).Graph;                                % to assign the modular graph to variable H
    H5.Edges.Weight=1./ H5.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T15= minspantree(H5);                                     % to compute the MST
    edge5=T15.Edges.EndNodes;                                 % to store the edges of the spanning tree
    for p5=1:length(edge5)
        s5=edge5(p5,1);
        t5=edge5(p5,2);
        position5(p5)=findedge(ModularGraph(j).Graph,s5,t5);    % to find the edge values in the original Modular graph
        weight5=ModularGraph(j).Graph.Edges(position5,2);    % to store the edge values
    end
    ModularGraph(j).sumMST=sum(table2array(weight5));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(j).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(j).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(j).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(j).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty end
end
%% Back
for i=1:n
    numberROIs=RangeROIs(i);
    [module,numberModules,xcentroid,ycentroid,zcentroid]=clusterNodesBACK(numberROIs,xcoord,ycoord,zcoord,NG,p);
    [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module);
    j=5*n+i;
    ModularGraph(j).numberROIs=numberROIs;
    ModularGraph(j).Graph=ModularG;
    ModularGraph(j).Matrix=ModularMatrix;
    ModularGraph(j).xcentroid=xcentroid;
    ModularGraph(j).ycentroid=ycentroid;
    ModularGraph(j).zcentroid=zcentroid;    
    ModularGraph(j).numberModules=numberModules;
    ModularGraph(j).ModulesComposition=module;
    ModularGraph(j).AverageNodePerModule=mean(cell2mat(module(:,2)));
%% CC, PL and SWP
    ModularGraph(j).NetCC=avg_clus_matrix(ModularMatrix,'O');                         % to compute NetCC
    ModularGraph(j).NetPL=avg_path_matrix(1./ModularMatrix);                          % to compute NetPL
    ModularGraph(j).RandomMatrix=randomize_matrix(ModularGraph(j).Matrix);            % to compute Random matrix
    ModularGraph(j).RegularMatrix=spatialLattice_matrix3d(ModularGraph(j).numberModules,ModularGraph(j).Matrix,xcentroid,ycentroid,zcentroid); % to compute Lattice Matrix based on spatial arguments
    ModularGraph(j).RandCC= avg_clus_matrix(ModularGraph(j).RandomMatrix, 'O');       % to compute RandCC
    ModularGraph(j).RegularCC= avg_clus_matrix(ModularGraph(j).RegularMatrix, 'O');   % to compute RegularCC
    ModularGraph(j).RandPL=avg_path_matrix(1./ModularGraph(j).RandomMatrix);          % to compute RandPL
    ModularGraph(j).RegularPL=avg_path_matrix(1./ModularGraph(j).RegularMatrix);      % to compute RegularPL
    z = (ModularGraph(j).NetPL - ModularGraph(j).RandPL);
    if z < 0
       z = 0;
    end
    diff_path =  z/ (ModularGraph(j).RegularPL - ModularGraph(j).RandPL);
    if ModularGraph(j).NetPL == Inf || ModularGraph(j).RandPL == Inf || ModularGraph(j).RegularPL == Inf
        diff_path = 1;
    end
    if diff_path > 1
        diff_path = 1;
    end

    B = (ModularGraph(j).RegularCC - ModularGraph(j).NetCC);
    if B < 0
        B = 0;
    end
    diff_clus = B / (ModularGraph(j).RegularCC - ModularGraph(j).RandCC);
    if isnan(ModularGraph(j).RegularCC) || isnan(ModularGraph(j).RandCC) || isnan(ModularGraph(j).NetCC)
        diff_clus = 1;
    end
    if diff_clus > 1
        diff_clus = 1;
    end
    ModularGraph(j).SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
    %% Modularity
    ModularGraph(j).Q=modularity_louvain_und(ModularGraph(j).Matrix);
%% Sum Degree    
    ModularGraph(j).SumDegree=sum(ModularGraph(j).Matrix,'all');
%% Sum of the MST    
    H6=ModularGraph(j).Graph;                                % to assign the modular graph to variable H
    H6.Edges.Weight=1./ H6.Edges.Weight;                      % to compute the inverse of each edge (in this way, we compute the maximum rather than the minimum spanning tree)
    T16= minspantree(H6);                                     % to compute the MST
    edge6=T16.Edges.EndNodes;                                 % to store the edges of the spanning tree
    for p6=1:length(edge6)
        s6=edge6(p6,1);
        t6=edge6(p6,2);
        position6(p6)=findedge(ModularGraph(j).Graph,s6,t6);    % to find the edge values in the original Modular graph
        weight6=ModularGraph(j).Graph.Edges(position6,2);     % to store the edge values
    end
    ModularGraph(j).sumMST=sum(table2array(weight6));     	% to compute the sum of the edges composing the MST
%% Sync
    A=ModularGraph(j).Matrix;
    dvector = sum(A,1);
    D = diag(dvector);                                   % convert it into a diagonal matrix D
    L = D - A;                                           % Calculate Laplacian L
    e = eig(L);                                          % Compute the eigenspectrum of L(t)
    e = sort(e);                                         % Sort smallest to largest eigenvalue
    ModularGraph(j).secondsmallest=e(2);                 % to store second smallest eigenvalue
    ModularGraph(j).largest=e(end);                      % to store largest eigenvalue
    ModularGraph(j).syncMod = abs(e(2)/e(end));          % Calculate synchronizabilty end
end
end