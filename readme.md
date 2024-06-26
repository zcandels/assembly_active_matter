# Project Assembly Active Matter

The purpose of this project is to be able to extend assembly theory from
static objects (eg molecules) to dynamic structures like clusters of 
organisms or active particles. We wish to determine the circumstances under which
the assembly of a system can increase.

# Model

Here we use a simple Vicsek model for our active particles. In the Vicsek model, 
each particle tries to orient its direction of motion with those of the particles 
a certain distance away from it. This tendency, however, is frustrated by thermal
noise which pushes the particle in random directions. The governing equations 
(expressed in discrete time) can be written as 

$$ \mathbf{x}_i^{(n+1)} = \mathbf{x}_i^{(n)} + \mathbf{v}_i^{(n)}\delta t $$

where 

$\mathbf{v}_i^{(n)} = v\cdot \left[\cos \theta^{(n)} \atop 
\sin \theta^{(n)}\right$. Here, $v$ is a constant velocity and the direction of motion
is determined in each timestep by

$$ \theta^{(n+1)} = \langle \theta^{(n)} \rangle_r + \Delta \theta $$

where $\langle\theta^{(n)} \rangle_r$ is the average direction of the particles
within a distance $r$ of particle $i$ and $\Delta \theta$ represents 
delta-correlated Gaussian noise.
