function LattMatrix=spatialLattice_matrix(numberModules,RewiredMatrix,xcentroid,ycentroid)
for i=1:numberModules             
    for j=2:numberModules         
        xDist= (xcentroid(i)-xcentroid(j))^2;                         % Calcuate squared distance between x values 
        yDist= (ycentroid(i)-ycentroid(j))^2;                         % Calcuate squared distance between y values
        distMatrix(i, j)= sqrt(xDist+yDist);                    % Calculate total distance using the pythagorean theorem
        distMatrix(j, i)= sqrt(xDist+yDist); 
    end
end
for nodes=1:numberModules
distvector=distMatrix(nodes,[nodes+1:numberModules]);  % to extract all the distance values for all the nodes with a index higher than node
orderVector=sort(distvector);               % to order distances in descending order
for idx=1:length(orderVector)
    closer(idx)=find(distvector==orderVector(idx)); %to identify the index of the nodes from the closer to the further
end
closer=closer+nodes;                                 % since distvector starts at nodes+1, the index of the nodes is shifted by the value of nodes
matrixvector=RewiredMatrix(nodes,[nodes+1:numberModules]);      % to extract all the edges  for all the nodes with a index higher than node; in this way, I am not going to break symmetry and the previous assignments 
matrixvector(matrixvector==0)=[];                    
orderMatrix=sort(matrixvector,'descend');            % to sort the edges from the stronger to the weaker
lattice=zeros(1,numberModules);                                % preallocation vector
closerToPlace=closer(1:length(orderMatrix));        % extract the n closer nodes, where n is the number of edges for the row
lattice(closerToPlace)=orderMatrix;                 % assign the stronger edges to the closer nodes
if nodes==numberModules
latticetocopy=lattice(numberModules);
else
    latticetocopy=lattice(nodes+1:numberModules);
end
LattMatrix(nodes,[nodes+1:numberModules])=latticetocopy;    % create row
LattMatrix([nodes+1:numberModules],nodes)=latticetocopy'; % create symmetric column
end
end