import numpy as np
import matplotlib.pyplot as plt
from rdfpy import rdf

Np = 25
Nx = int(np.sqrt(Np))
Np = Nx**2
L_x = 2
x_vals = np.linspace(-L_x/2, L_x/2, Nx)
y_vals = np.linspace(-L_x/2, L_x/2, Nx)
rho = Np/L_x**2



Particles = {}
for i in range(Nx):
    for j in range(Nx):
        Particles[i + j*(Nx)] = np.array([ x_vals[i], y_vals[j] ])
        
dist_array = []
r_cutoff = 0.1
bins = np.arange(0, L_x, r_cutoff)
rdf = np.zeros(len(bins) -1)


for i in range(Np):
    dist_array = []
    p_pos = Particles[i]
    for j in range(Np):
        q_pos = Particles[j]
        if i == j:
            continue
        pq_dist = p_pos - q_pos
        dist_array.append(np.linalg.norm(p_pos - q_pos))
    rdf += np.histogram(dist_array, bins)[0]
    #rdf = np.divide(rdf, bins[1:]**2)
rdf = np.divide(rdf, bins[1:]**2)
rdf /= (Np)

rdf /= (4*np.pi*rho)

#%%
plt.figure()
plt.plot(bins[1:], rdf)
plt.xlabel(r"$r$", fontsize=20)
plt.ylabel(r"$g(r)$", fontsize=20)
    
        
        
        
        