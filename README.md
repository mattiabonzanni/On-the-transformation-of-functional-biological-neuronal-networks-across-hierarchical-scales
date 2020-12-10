# On-the-transformation-of-functional-biological-neuronal-networks-across-hierarchical-scales

WHY: 
The brain possesses structural and functional hierarchical architectures organized over multiple scales. Considering that functional recordings commonly focused on a single spatial level, and since multiple scales interact with one another, we should explore the behavior of neuronal networks across different scales. 

HOW: We can study the behavior of a given network across scales applying a renormalization approach. In the codes, the renormalization approach is based on the Euclidian distance between nodes.

WHAT: In the simulation, the following parameters are studied across scales: 
1. Average clustering coefficient
2. Average path length
3. Small-world propensity
4. Average edge weight
5. Modularity
6. Synchronizability 
7. Fraction of long-term connections

There are three folders:

Complete surface tessellation - 
Aim1: To generate a 2D network (Global Network) and, based on spatial arguments,completelly renormalize it (Modular Network). The complete renormalization implies a tessellation of the entire network.  
Aim2: To study the network properties between the Global vs Modular networks.

Partial volume tessellation - 
Aim1: To generate a 3D network (Global Network) and, based on spatial arguments,partially renormalize it (Modular Network). The partial renormalization implies a tessellation of a fraction of the entire network. 
Aim2: To study the network properties between the Global vs Modular networks.

Vertex Simulation - 
To study the network properties of microscale vs mesoscale signals from the same artificial neuronal network. Please refer for Vertex simulations to (developed by Tomsett RJ): http://vertexsimulator.org/

Each folder contains:

1. All the required codes
2. A code to run the simulation (either as single- to run a single set of inputs- or as a loop - to run several combinations of inputs-). Use the following codes for the simulations:
    2.1 Complete surface tessellation
        2.1.1 CompleteSurfaceTessellation.m
        2.1.2 LoopCompleteSurfaceTessellation.m
    2.2 Partial Volume Tessellation
        2.2.1 PartialVolumeTessellation.m
        2.2.2 Loop_PartialVolumeTessellation.m
    2.3 Vertex Simulation (unzip the folder)
        2.3.1 LoopLFPNetwork.m
All the simulations are automatically saved in the folder at the end with the data and time. 
