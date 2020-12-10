clear
clc
RangeROIs=[8].^2;
myList=string(ls);
isCond= contains(myList,'NG');
NG=1000;
file=myList(isCond);
k=1;
for t=1:length(file)
    load(file(t))
    for i=1:length(RangeROIs)
        numberROIs=RangeROIs(i);
        [module,numberModules,xcentroid,ycentroid,zcentroid,HiddenNodes]=clusterNodesBOTTOM(numberROIs,Results.Input.xcoord,Results.Input.ycoord,Results.Input.zcoord,NG,0.2);
        [ModularRndG,ModularRandMatrix]=createModularGraph(Results.Output.GlobalGraph.RandomMatrix,numberModules,module);
        H=ModularRndG;
        T = minspantree(ModularRndG);
        SpanningTree(k).NG=Results.Input.NG;
        SpanningTree(k).Beta=Results.Input.betarewiring;
        H.Edges.Weight=1./ H.Edges.Weight;
        T1= minspantree(H);
        edge=T1.Edges.EndNodes;
        for p=1:length(edge)
            s=edge(p,1);
            t=edge(p,2);
            position(p)=findedge(ModularRndG,s,t);
            for f=1:length(position)
            if position(f)>size(H.Edges)
                position(f)=0;
            end
            end
            position(position==0)=[];
            weight=ModularRndG.Edges(position,2);
        end
        SpanningTree(k).valueSum2=sum(table2array(weight));
        k=k+1
    end
end

