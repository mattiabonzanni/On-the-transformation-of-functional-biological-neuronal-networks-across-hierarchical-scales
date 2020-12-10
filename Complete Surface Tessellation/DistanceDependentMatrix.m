% Required codes
% 1.1 avg_clus_matrix.m;
    % 1.1.1 clustering_coef_matrix.m;
% 1.2 randomize_matrix.m
% 1.3 avg_path_matrix.m;
% 1.4 modularity_louvain_und.m;

%Inputs:
% 1. NG= total number of seeded nodes in the whole graph;
% 2. beta= probability that an edge is randomly re-wired to a non-neighbor  node;
% 3. RangeROIs= number of Region of Interests (ROIs) used to tessellate the whole graph;
% 4. edgeDensity= the average percentage of neighbor nodes in the whole graph; 

% A Results structure containing;
% Field1 - Output
    % GlobalGraph
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


function Results=DistanceDependentMatrix(NG,edgeDensity,beta)
xcoord=rand(1,NG);                                              % x-axis coordinates between [0;1]
ycoord=rand(1,NG);                                              % y-axis coordinates between [0;1]
distMatrix=zeros(NG);                                           % create an empty matrix
for i=1:1:NG             
    for j=2:1:NG         
        xDist= (xcoord(i)-xcoord(j))^2;                         % Calcuate squared distance between x values 
        yDist= (ycoord(i)-ycoord(j))^2;                         % Calcuate squared distance between y values
        distMatrix(i, j)= sqrt(xDist+yDist);                    % Calculate total distance using the pythagorean theorem
        distMatrix(j, i)= sqrt(xDist+yDist); 
    end
end
distMatrix=distMatrix-(max(distMatrix,[],'all'));               % to have closer nodes with higher weight
distMatrix=distMatrix./min(distMatrix,[],'all');                % to range every between [0;1]. It follows that we can have the threshold's range [0;1].
%% Pruning
wt=0.6;                                                         % initial value of wt; the value will be updated until reaching the desired edge density
Degree=NG*(NG-1);
wtincrement=0.001;                                              % arbitrarly defined increment of wt for the iteration
while Degree>(2*edgeDensity*NG)/100                             % to identify the value of wt necessary to achieve the desired edge density
    distMatrix(distMatrix<=wt)=0;                               % to threshold the matrix in order to eliminate edges from nodes ij given the distance(ij)<threshold. 
    distMatrix = distMatrix - diag(diag(distMatrix));           % to eliminate self-loops (0 values on the diagonal)
    A=graph(distMatrix);                                    
    [s,t]=findedge(A);                                          % to extract the s and t table of the nodes pair
    weights=A.Edges.Weight(findedge(A,s,t));                    % to extract the edgeWeights table of the nodes pair
    NodeDegree=degree(A);                                       
    Degree=mean(NodeDegree);
    wt=wt+wtincrement;
    fprintf('wt=%.3f\n',wt);
end
LatticeGraph = graph(s,t,weights);
Results.Output.GlobalGraph.LatticeMatrix=full(adjacency(LatticeGraph,'weighted'));   % Spatial Lattice Matrix
Results.Output.GlobalGraph.wt=wt-wtincrement;
fprintf('Pruning step complete\n');
%% Rewiring
source=rand(size(t,1),1);                                       % to generate random numbers between [0;1]; we need to generate as many as are the edges, thus I need to have size of table t
source(source<=beta)=1;                                         % to generate logic 1 values for numbers with probability equal or smaller than beta
source(source<1)=0;                                             % to generate logic 0 values for numbers with probability larger than beta
nodesToRewire=find(source);                                     % to find the index positions of the logic 1 values
howManyToRewire=size(nodesToRewire,1);                          % to find how many edges will be rewired
for q=1:howManyToRewire
    g=graph(s,t);                                               % creating the graph in the loop allows to consider at every iteration the new targets for s in order to avoid multiple edges between pairs of nodes
    AllTargets=(1:NG);                                          % to list all the nodes
    position=nodesToRewire(q,1);                                % to substitute an edge in t using the position based on nodesToRewire
    nodeS=s(position,1);
    oldTargets=[transpose(neighbors(g,nodeS))];                 % to find the node in table s
    oldTargets=horzcat(oldTargets,nodeS);                       % to exclude self loops
    AllTargets(oldTargets)=[];                                  % to exclude old targets and thus generate a new list of avaible targets for rewiring
    newTargetidx=randi(size(AllTargets,2),1);                   % to generate random numbers to choose the new targets for rewiring from the reduced AllTargets array
    t(position,1)=AllTargets(1,newTargetidx);                   % to rewire in table t without multiple edges using a random value from the reduced AllTargets array  
end
Results.Output.GlobalGraph.Graph = graph(s,t,weights);
Results.Output.GlobalGraph.Matrix=adjacency(Results.Output.GlobalGraph.Graph,'weighted');
fprintf('Re-wiring step complete\n');
%% CC and PL computation
Results.Output.GlobalGraph.NetCC=avg_clus_matrix(Results.Output.GlobalGraph.Matrix, 'O');           % to compute NetCC
Results.Output.GlobalGraph.NetPL=avg_path_matrix(1./Results.Output.GlobalGraph.Matrix);             % to compute NetPL
Results.Output.GlobalGraph.RegCC=avg_clus_matrix(Results.Output.GlobalGraph.LatticeMatrix, 'O');    % to compute RegCC
Results.Output.GlobalGraph.RegPL=avg_path_matrix(1./Results.Output.GlobalGraph.LatticeMatrix);      % to compute RegPL
Results.Output.GlobalGraph.RandomMatrix = randomize_matrix(Results.Output.GlobalGraph.Matrix);      % to generate the Random null model
Results.Output.GlobalGraph.RandCC=avg_clus_matrix(Results.Output.GlobalGraph.RandomMatrix, 'O');    % to compute RandCC
Results.Output.GlobalGraph.RandPL=avg_path_matrix(1./Results.Output.GlobalGraph.RandomMatrix);      % to compute RandPL
%% SWP computation
D = (Results.Output.GlobalGraph.NetPL - Results.Output.GlobalGraph.RandPL);
if D < 0
    D = 0;
end
diff_path =  D/ (Results.Output.GlobalGraph.RegPL - Results.Output.GlobalGraph.RandPL);
if Results.Output.GlobalGraph.NetPL == Inf || Results.Output.GlobalGraph.RandPL == Inf || Results.Output.GlobalGraph.RegPL == Inf
    diff_path = 1;
end
if diff_path > 1
    diff_path = 1;
end

B = (Results.Output.GlobalGraph.RegCC - Results.Output.GlobalGraph.NetCC);
if B < 0
    B = 0;
end
    
diff_clus = B / (Results.Output.GlobalGraph.RegCC - Results.Output.GlobalGraph.RandCC);
if isnan(Results.Output.GlobalGraph.RegCC) || isnan(Results.Output.GlobalGraph.RandCC) || isnan(Results.Output.GlobalGraph.NetCC)
    diff_clus = 1;
end
if diff_clus > 1
    diff_clus = 1;
end
Results.Output.GlobalGraph.SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));             % to compute SWP
%% Modularity computation
Results.Output.GlobalGraph.Q=modularity_louvain_und(Results.Output.GlobalGraph.Matrix);     % to compute modularity
%% Sum Degree computation 
Results.Output.GlobalGraph.SumDegree=sum(sum(Results.Output.GlobalGraph.Matrix));           % to compute the sum of the degree values
%% Sync computation
A=Results.Output.GlobalGraph.Matrix;
dvector = sum(A,1);
Dia = diag(dvector);                                    % convert it into a diagonal matrix D
L = Dia - A;                                            % Calculate Laplacian L
e = eig(L);                                             % Compute the eigenspectrum of L(t)
e = sort(e);                                            % Sort smallest to largest eigenvalue
Results.Output.GlobalGraph.secondsmallest=e(2);         % to store second smallest eigenvalue
Results.Output.GlobalGraph.largest=e(end);              % to store largest eigenvalue
Results.Output.GlobalGraph.sync = abs(e(2)/e(end));     % Calculate synchronizabilty 
% to store the node coordinates
Results.Input.xcoord=xcoord;
Results.Input.ycoord=ycoord;
end