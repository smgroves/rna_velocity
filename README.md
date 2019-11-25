# Landscape Theory and RNA Velocity

This repository contains code used for interrogating the methods in RNA velocity analysis. We use simple network models to derive CME-based equations for each monomer in the model to use as a baseline for comparing the accuracy of RNA-velocity-derived metrics.

## Building the network model and generating fake data

We have implemented two models so far: 1 self-regulating gene and a 2 gene model with mutual inhibition and self-activation. Both of these models at the current parameters have 2 steady states. 

## Running RNA velocity

We use the original RNA velocity package from La Manno et al (2018) [`velocyto`] to calculate RNA velocity for each gene for each cell in the fake dataset. Other methods (such as scVelo) may be implemented in the future.

We initialize `N` cells randomly and then let the system evolve under the dynamics of the network model. With multiple fake timepoints, we can compare three things:
1. Extrapolated t1, based on RNA velocity at t0, *vs* t1, based on second time point in dataset
2. Extrapolated t1, based on network model dynamics, *vs* t1, based on second time point in dataset
3. Extrapolated t1, based on RNA velocity, *vs* extrapolated t1, based on network model dynamics.

Comparison # 1 is done in the original La Manno paper, and will be used to investigate the accuracy of RNA velocity based on the actual evolution of the system. Comparison # 2 should give us a very small error, and this will give us a reference for the magnitude of error to expect based solely on stochasticity in the model. Comparison # 3 is what we are actually interested in, and we will make this comparison by interrogating the transition matrices derived from the network model and from RNA velocity between all time points. This comparison thus requires that we build an RNA-velocity-derived transition matrix by taking the cosine correlation between sampled points and predicted velocity, thus converting a sparse velocity field into a Markov chain model. This step in and of itself may cause errors due to the limitations of projecting velocity vectors onto pre-defined discrete sampled points, which we may choose to investigate at a later time.

## Quantifying plasticity

One of the goals of this project is to provide a theoretical basis for quantifying plasticity from RNA velocity. Plasticity has not been well defined in the literature, and it is often used interchangeably with heterogeneity of a system, which alternatively does not actually characterize the dynamics of the system in the way plasticity should. We therefore want to quantify multiple components of plasticity:
1. Plasticity should describe how far a cell is likely to travel from it's current cell state. This should correspond to some notion of instability in dynamical systems theory, and we can compare these two metrics using the relative instability of states and the plasticity metric we derive from RNA velocity transitions. In this quantification, we may expect a transit-amplifying cell to look more plastic than a stem cell, depending on the local stability of each state in the landscape.
2. Plasticity should describe how variable a cell's fate is; i.e. it should correlate with the degree of multipotency of a cell type. In this case, a stem cell that has multiple possible fates should be considered more plastic than a differentiated cell, regardless (or in despite of) of proliferation rates. 
3. Plasticity most likely correlates with the openness of chromatin, which can be modeled in our system as a multiplicative factor describing the likelihood of TF binding (i.e. closed chromatin has lower likelihood of TF binding). By changing the parameters of the system, we should be able to increase or decrease plasticity. 
