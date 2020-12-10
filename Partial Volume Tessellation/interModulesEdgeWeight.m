function finaledge=interModulesEdgeWeight(G0,G,module1,module2);
%% To compute the edge weight between two modules based on the shorthest path
% The assumption is that the edge between two modules is equal to the sum
% of the product of the edges in the ij pair shorthest paths (where i and j belongs to module1 and module2,
% respectively) divided by the number of the total number of paths (which
% is the product between the number of nodes in module1 and in module2).

% Input:
% 1.
% 2. module1
% 3. module2

% Output:
% 1. finaledge= edge weight between module1 and module2

weightShortPath=ones(length(module1),length(module2));      % matrix preallocation
for source=1:length(module1)                                % to scan trought the nodes in module1
    for target=1:length(module2)                            % to scan through the nodes in module2
        [P,d,edgepath] = shortestpath(G0, module1(source), module2(target));          % to identify the shortest path between source-target
        IndirectEdgeWeightVector=G.Edges.Weight(edgepath);                            % to create a vector with all the edge weigths that connect source and target. 
        weightShortPath(source,target)=prod(IndirectEdgeWeightVector);                % to multiply the edge weight to calculate the weight of the short path between source and target. You substitute the product of the weights of the path if a direct edge is not existing    
    end
end
numberOfPaths=(length(module1)*length(module2));          % total number of paths between moduel1 and module2
finaledge=sum(sum(weightShortPath))/numberOfPaths;        % to find the final edge between module1 and module2, divides the sum of the weightShortPath by the total number of paths which composes the two modules. 
end