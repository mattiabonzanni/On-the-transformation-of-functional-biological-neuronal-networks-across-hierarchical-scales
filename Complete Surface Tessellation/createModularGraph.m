function [ModularG,ModularMatrix]=createModularGraph(distMatrix,numberModules,module)
%% To create the modular graph 
% Input:
% 1. distMatrix=Global matrix;
% 2. numberModules= number of ROI in which more than one node is identified;
% 3. module= structure containing the composition of each module.

% Output
% 1. ModularG= weighted unidirected modular graph;
% 2. ModularMatrix= modular matrix.

% written by Mattia Bonzanni
%%
count=1;
for m=1:numberModules
    module1=module{m,1};                                        % to extract module1
    for n=m+1:numberModules
        module2=module{n,1};                                    % to extract module2
        finaledge=(sum(full(distMatrix(module1, module2)),'all'))/(length(module1)*length(module2)); % to compute the edge weight between two modules based on the direct path. I just extract the submatrix from distMatrix with the rows equal to the element of Module1 and the columns equal to the element of Module2
        table(count,1)=m;                                       % s column for the graph
        table(count,2)=n;                                       % t column for the graph
        table(count,3)=finaledge;                               % weight edge column for the graph
        count=count+1;
    end
    fprintf('Iteration ModularG compuation: %d/%d\n',m,numberModules);
end
ModularG=graph(table(:,1),table(:,2));                          % the modular graph
ModularG.Edges.Weight=table(:,3);                               % the edges of the modular graph
ModularMatrix=full(adjacency(ModularG,'weight'));               % the modular matrix
end