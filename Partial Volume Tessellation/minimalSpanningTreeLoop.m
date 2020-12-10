clear
clc
myList=string(ls);
isCond= contains(myList,'NG500');
file=myList(isCond);
k=1;
for i=1:length(file)
    load(file(i))
    if Results.Input.edgeDensity==8
        for j=1:length(Results.Output.ModularGraph);
        if Results.Output.ModularGraph(j).numberROIs==64;
            H=Results.Output.ModularGraph(j).GraphVisible;
        T = minspantree(Results.Output.ModularGraph(j).GraphVisible);
        SpanningTree(k).NG=Results.Input.NG;
        SpanningTree(k).value=sum(T.Edges.Weight);
        SpanningTree(k).Beta=Results.Input.betarewiring;
        H.Edges.Weight=1./Results.Output.ModularGraph(j).GraphVisible.Edges.Weight
        T1= minspantree(H);      
%         x=Results.Output.ModularGraph(1).xcentroid;
%         y=Results.Output.ModularGraph(1).ycentroid;
%         z=Results.Output.ModularGraph(1).zcentroid;
%         h=plot(Results.Output.ModularGraph(j).GraphVisible,'XData',x,'YData',y,'ZData',z)
%         h.MarkerSize=8
%         h.NodeLabel=[];
%         h.LineWidth=2;
%         h.NodeColor=[0.5 0.5 0.5];
%         h.EdgeColor=[0.5 0.5 0.5];
% %         highlight(h,T1,'EdgeColor','r','LineWidth',7)
%         set(gca,'XTick',[], 'YTick', [],'ZTick', [])
% 
        edge=T1.Edges.EndNodes;
        for p=1:length(edge)
            s=edge(p,1);
            t=edge(p,2);
            position(p)=findedge(Results.Output.ModularGraph(j).GraphVisible,s,t);
            for f=1:length(position)
            if position(f)>size(H.Edges)
                position(f)=0;
            end
            end
            position(position==0)=[];
            weight=Results.Output.ModularGraph(j).GraphVisible.Edges(position,2);
        end
        SpanningTree(k).valueSum2=sum(table2array(weight));
        k=k+1
        end
        end
end
end
