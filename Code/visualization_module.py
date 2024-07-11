import numpy as np
import matplotlib.pyplot as plt
import os

def visualize_clusters(cluster_dict, step):
    
    
    import matplotlib.cm as cm
    N_clus = len(cluster_dict)
    start = 0.0
    stop = 1.0
    cm_subsection = np.linspace(start, stop, N_clus)
    colors = [ cm.jet(x) for x in cm_subsection ]
    
    
    plt.rcParams['text.usetex'] = True
    plt.figure()
    for n in range(1, N_clus):
        cluster = cluster_dict[n]
        clusterDataFrame = cluster.clusterDataFrame
        x = clusterDataFrame['x-position']
        y = clusterDataFrame['y-position']
        theta = clusterDataFrame['angle']
        u = np.cos(theta)
        v = np.sin(theta)
        plt.quiver(x,y,u,v,color=colors[n])
        plt.xlabel(r'$x$', fontsize=20)
        plt.ylabel(r'$y$', fontsize=20)
        plt.title('timestep = %d' %step)
    plt.show()
    
    path_to_code = os.getcwd()
    assembly_active_matter_path = os.path.dirname(path_to_code)
    figPath = assembly_active_matter_path + "/figures/data"
    plt.savefig(figPath + f"step{step}.png", bbox_inches='tight')
    
    return
    
def cluster_histogram(cluster_dict):
    particle_distribution = []
    for _, cluster in cluster_dict.items():
        numParticles = cluster.numParticles
        particle_distribution.append(numParticles)
    mean = np.mean(particle_distribution)
    std = np.std(particle_distribution)
    
    plt.rcParams['text.usetex'] = True
    fig, ax = plt.subplots()
    ax.hist(particle_distribution, bins="auto", rwidth=0.8)
    ax.set_xlabel(r"Particles per Cluster")
    ax.set_ylabel(r"Count")
    ax.annotate(r'$\mu = %d$' %mean, xy=(0.95, 0.9), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    ax.annotate(r'$\sigma = %g$' %std, xy=(0.95, 0.8), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    #ax.annotate(r"$\sigma = %g$" % variance,  xytext=(0.75, 0.7), # fraction, fraction
    #        textcoords='figure fraction',fontsize=15)

    return

def object_histogram(object_dict, sim_vicsek):
    eta = sim_vicsek.eta
    assembly_indices = []
    for _, obj in object_dict.items():
        assembly_indices.append(obj.assemblyIndex[-1])
    mean = np.mean(assembly_indices)
    std = np.std(assembly_indices)
    plt.rcParams['text.usetex'] = True
    fig, ax = plt.subplots()
    ax.hist(assembly_indices, bins="auto", rwidth=0.8)
    ax.set_xlabel(r"Assembly Indices")
    ax.set_ylabel(r"Count")
    ax.annotate(r'$\mu = %d$' %mean, xy=(0.95, 0.9), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    ax.annotate(r'$\sigma = %g$' %std, xy=(0.95, 0.8), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    
    
    path_to_code = os.getcwd()
    assembly_active_matter_path = os.path.dirname(path_to_code)
    figPath = assembly_active_matter_path + "/figures/data"
    plt.savefig(figPath + f"hist_eta_{eta}.png", bbox_inches='tight')
    
    
        

def mean_assembly_ind(assembly_mean_var, sim_vicsek):
    
    eta = sim_vicsek.eta
    rho = sim_vicsek.numParticles/sim_vicsek.L**2
    
    mean_ass_ind = assembly_mean_var[:, 0]
    var = assembly_mean_var[:, 1]
    time = np.asarray(range(len(mean_ass_ind)))
    plt.rcParams['text.usetex'] = True
    fig, ax = plt.subplots()
    ax.plot(mean_ass_ind, linewidth=2)
    #ax.fill_between(time, mean_ass_ind - var, mean_ass_ind + var,
    #                alpha=0.2, label='std dev')
    ax.set_xlabel(r"time $t$", fontsize=20)
    ax.set_ylabel(r"$\bar{a}$",fontsize=20)
    ax.annotate(r'$\eta = %g$' %eta, xy=(0.95, 0.7), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    ax.annotate(r'$\rho = %g$' %rho, xy=(0.95, 0.62), xycoords='axes fraction',
            fontsize=20, ha='right', va='top')
    ax.tick_params(direction='in')
    
    figPath = "C:/Users/2941737C/Research/"\
        + "assembly_active_matter/figures/data/"
        
    path_to_code = os.getcwd()
    assembly_active_matter_path = os.path.dirname(path_to_code)
    figPath = assembly_active_matter_path + "/figures/data/"
        
    plt.savefig(figPath + f"mean_assembly_eta_{eta}.png", bbox_inches='tight')