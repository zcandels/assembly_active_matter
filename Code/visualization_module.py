import numpy as np
import matplotlib.pyplot as plt

def visualize_clusters(cluster_dict, step):
    
    import matplotlib.cm as cm
    N_clus = len(cluster_dict)
    start = 0.0
    stop = 1.0
    cm_subsection = np.linspace(start, stop, N_clus)
    colors = [ cm.jet(x) for x in cm_subsection ]
    
    
    plt.rcParams['text.usetex'] = True
    plt.figure()
    for n in range(N_clus):
        cluster = cluster_dict[n]
        clusterDataFrame = cluster.clusterDataFrame
        x = clusterDataFrame['x-position']
        y = clusterDataFrame['y-position']
        theta = cluster['angle']
        u = np.cos(theta)
        v = np.sin(theta)
        plt.quiver(x,y,u,v,color=colors[n])
        plt.xlabel(r'$x$', fontsize=20)
        plt.ylabel(r'$y$', fontsize=20)
    plt.show()
    
    figPath = "C:/Users/2941737C/Research/assembly_active_matter/figures/data/"
    plt.savefig(figPath + f"step{step}.png", bbox_inches='tight')
    
    return
    
def cluster_histogram(cluster_dict):
    particle_distribution = []
    for _, cluster in cluster_dict.items():
        numParticles = cluster.numParticles
        particle_distribution.append(numParticles)
    mean = np.mean(particle_distribution)
    variance = np.var(particle_distribution)
    
    plt.rcParams['text.usetex'] = True
    fig, ax = plt.subplots()
    ax.hist(particle_distribution, bins="auto", rwidth=0.8)
    ax.set_xlabel(r"Particles per Cluster")
    ax.set_ylabel(r"Count")
    ax.annotate(r'$\mu = %d$' %mean, xy=(0.95, 0.9), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    ax.annotate(r'$\sigma = %g$' %variance, xy=(0.95, 0.8), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    #ax.annotate(r"$\sigma = %g$" % variance,  xytext=(0.75, 0.7), # fraction, fraction
    #        textcoords='figure fraction',fontsize=15)

    return