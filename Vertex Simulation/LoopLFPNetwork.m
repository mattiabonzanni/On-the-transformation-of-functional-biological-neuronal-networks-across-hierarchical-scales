%% Loop Vertex Simulation
% Aim: To study the network properties of microscale vs mesoscale signals
% from the same artificial neuronal network. 
%% Nomenclature:
% - Global (or Whole) graph: to mimic microscale neuronal networks (computed from neuronal signals);
% - Modular Graph: to mimic mesoscale neuronal networks (computed from LFP signals);
% - Regular (or Lattice): regular null model (each node is attached to its neighbours);
% - Random: random null model (random wiring);
% - Net: the studied network.
%% Required codes
% 1. avg_clus_matrix.m;
    % 1.1 clustering_coef_matrix.m;
% 2. randomize_matrix.m;
% 3. avg_path_matrix.m;
% 4. tutorial5LFP
%% Inputs-Outputs
        %Inputs:
% 1. density= TissueParams.neuronDensity from tutorial5LFP ;
% 2. a= depth of MEA from tutorial5LFP;
% 3. b= depth of MEA from tutorial5LFP;
% 4. cutoff= arbitrary clustering coefficient cutoff to define a which neuron influence a given electrode; 
        %Outputs:
% A Network structure containing;
% Field1 - Matrix
% 1. NeuronsLFPCorr= correlation matrix between each neuron and electrodes;
% 2. WholeMatrix= correlation matrix of the single neurons;
% 3. LFPMatrix= correlation matrix of the single LFP signals;
% 4. WholeMatrixRescaled= rescaled WholeMatrix (to avoid negative values);
% 5. LFPMatrixRescaled= rescaled LFPMatrix

% Field2 - Input
% See Vertex instructions for the definition of the variables
% 1. spikes
% 2. LFP: electrode traces - mesoscale
% 3. v_m: voltage traces - microscale
% 4. l_syn
% 5. params
    % 5.1 TissueParams
    % 5.2 NeuronParams
    % 5.3 ConnectionParams
    % 5.4 RecordingSettings
    % 5.5 SimulationSettings

% Field 3 - SummaryParameters
% Table format
% 1. Neurons= number of neurons;
% 2. Neurons/module= average density of neuron per electrode;
% 3. Electrodes= number of electrodes;
% The next parameters are for both Whole and Modular networks
% 4-5. CC= average clustering coefficient;
% 6-7. PL= average path length;
% 8-9. SumWeight= the sum of the weighted edges;
% 10. SecondRatio= ratio between the second smallest eigenvalues of the Modular and Whole graphs
% 11. LargestRatio= ratio between the largest eigenvalues of the Modular and Whole graphs

% Field 4 - FileName
%% References codes
% 1) tutorial5LFP (http://vertexsimulator.org/tutorial-5-models-with-several-layers/)
% 1) avg_clus_matrix (written by Eric Bridgeford);
% 2) avg_path_matrix (written by Eric Bridgeford);
% 3) clustering_coef_matrix (code originally written by Mika Rubinov,UNSW, 2007-2010 and modified/written by Eric Bridgeford);
% 4) latmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);
% 5) randmio_und_connected (Mika Rubinov, UNSW; Jonathan Power, WUSTL and Olaf Sporns, IU);

% Written by Mattia Bonzanni 

clear
clc
densityRange=100;
aRange=[1000];
bRange=-100;
cutoff=0.6;
varNames={'Neurons','Neurons/module','Electrodes','CCWhole','CCLFP','PLWhole','PLLFP','SumWeightWhole','SumWeightLFP','SecondRatio','LargestRatio'};
for Q=1:length(densityRange)
    for E=1:length(bRange)
        Density=densityRange(Q);
        a=aRange;
        b=bRange(E);
        Results=tutorial5LFP(Density,a,b);
        %% correlation-based cutoff
        dimNeu=size(Results.v_m,1);                             % number of neurons
        dimLFP=Results.params.RecordingSettings.numElectrodes;  % number of electrodes
        sz=[dimNeu,dimLFP];
        NeuronsLFPcorr=zeros(sz);                               % preallocation matrix
        for i=1:dimNeu
            for j=1:dimLFP
                corroeficient=corrcoef(Results.v_m(i,:),Results.LFP(j,:));
                NeuronsLFPcorr(i,j)=abs(corroeficient(1,2));                % contribution of each neuron to the resulting LFP signal per each electrode; 
            end
        end
        for l=1:dimLFP
            trace=NeuronsLFPcorr(:,l);
            minT=min(trace);
            maxT=max(trace);
            tracerenorm=(trace-minT)/(maxT-minT);
            module(l)=sum(tracerenorm>cutoff);  % cell variable in which I store the composition of each module
            fprintf('column %d\n',l)
        end
        AverageNeurons=mean(module);
        %% Correlation Coefficient matrices
        WholeMatrix=zeros(dimNeu);  % preallocation matrix
        % load WholeMatrix
        for n=1:dimNeu
            trace1=Results.v_m(n,:);    %first trace
            for m=n+1:dimNeu
                cc=corrcoef(trace1,Results.v_m(m,:));
                WholeMatrix(n,m)=cc(1,2);
                WholeMatrix(m,n)=cc(1,2);
                fprintf('Whole Matrix n=%d - m=%d\n',n,m)
            end
        end
        for n=1:dimLFP-1
            for m=n+1:dimLFP
                cc=corrcoef(Results.LFP(n,:),Results.LFP(m,:));
                LFPMatrix(n,m)=cc(1,2);
                LFPMatrix(m,n)=cc(1,2);
            end
        end
        %% Rescale matrices
        minWM=min(min(WholeMatrix));
        maxWM=max(max(WholeMatrix));
        WholeMatrixRescaled=zeros(length(WholeMatrix));
        for i=1:length(WholeMatrix)
            for j=i+1:length(WholeMatrix)
                WholeMatrixRescaled(i,j)=(WholeMatrix(i,j)-minWM)/(maxWM-minWM);
                WholeMatrixRescaled(j,i)=WholeMatrix(i,j);
                fprintf('Whole Matrix Rescale i=%d - j=%d\n',i,j)
            end
        end
        minRM=min(min(LFPMatrix));
        maxRM=max(max(LFPMatrix));
        for i=1:length(LFPMatrix)
            for j=i+1:length(LFPMatrix)
                LFPMatrixRescaled(i,j)=(LFPMatrix(i,j)-minRM)/(maxRM-minRM);
                LFPMatrixRescaled(j,i)=LFPMatrixRescaled(i,j);
            end
        end
        %% Parameters
        net_path_LFP = avg_path_matrix(1./LFPMatrixRescaled);      % average path of the network
        net_clus_LFP = avg_clus_matrix(LFPMatrixRescaled,'O');     % average clustering coef. with the Onella method ('O')
        net_clus_Whole = avg_clus_matrix(WholeMatrixRescaled,'O');     % average clustering coef. with the Onella method ('O')
        net_path_Whole = avg_path_matrix(1./WholeMatrixRescaled);      % average path of the network
        SumWeightLFP=sum(sum(LFPMatrixRescaled));
        SumWeightWhole=sum(sum(WholeMatrixRescaled));
        dvectorWhole = sum(WholeMatrixRescaled,1);
        W = diag(dvectorWhole);
        %% Calculate Laplacian L
        LW = W - WholeMatrixRescaled;
        %% Compute the eigenspectrum of L(t)
        eW = eig(LW);
        %% Calculate synchronizabilty s
        % Sort smallest to largest eigenvalue
        eW = sort(eW);
        secondsmallestW=eW(2);
        largestW=eW(end);
        syncWhole=abs(secondsmallestW/largestW);
        dvectorLFP = sum(LFPMatrixRescaled,1);
        % convert it into a diagonal matrix D
        DLFP = diag(dvectorLFP);
        %% Calculate Laplacian L
        LLFP = DLFP - LFPMatrixRescaled;
        %% Compute the eigenspectrum of L(t)
        eLFP = eig(LLFP);
        %% Calculate synchronizabilty s
        % Sort smallest to largest eigenvalue
        eLFP = sort(eLFP);
        % Compute synchronizability: ratio of second smallest eigenvalue to largest
        % eigenvalue
        secondsmallestLFP=eLFP(2);
        largestLFP=eLFP(end);
        syncLFP = abs(eLFP(2)/eLFP(end));
        secondRatio=secondsmallestLFP/secondsmallestW;
        largestRatio=largestLFP/largestW;
        %% Organize Output
        Network.Matrix.NeuronsLFPcorr=NeuronsLFPcorr;
        Network.Matrix.WholeMatrix=WholeMatrix;
        Network.Matrix.LFPMatrix=LFPMatrix;
        Network.Matrix.WholeMatrixRescaled=WholeMatrixRescaled;
        Network.Matrix.LFPMatrixRescaled=LFPMatrixRescaled;
        Network.Input=Results;
        Network.SummaryParameters=table(dimNeu,AverageNeurons,dimLFP,net_clus_Whole,net_clus_LFP,net_path_Whole,net_path_LFP,SumWeightWhole,SumWeightLFP,secondRatio,largestRatio,'VariableNames',varNames);
        time=datestr(now,'mmddyy_HHMMSS');
        s=sprintf('NG%d_',dimNeu);
        filename=[s,time];
        Network.FileName=filename;
        save(filename,'Network') 
    end
end