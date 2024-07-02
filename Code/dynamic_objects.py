import numpy as np
import matplotlib.pyplot as plt
import clustering_fns as clus
import vic as vm


class cluster():
    def __init__(self, cluster_df):
        self.clusterDataFrame = cluster_df
        self.numParticles = cluster_df.shape[0]
    
    def computeCentroid(self):
        
        positions = self.clusterDataFrame[ ['x-position', 'y-position'] ]
        x_pos, y_pos = positions.iloc[:,0], positions.iloc[:,1]
        x_com, y_com = 0., 0.
        for i in range(len(x_pos)):
            # x_com += positions.iloc[i,0]
            # y_com += positions.iloc[i,1]
            x_com += x_pos.iloc[i]
            y_com += y_pos.iloc[i]
        x_com = x_com/(len(x_pos))
        y_com = y_com/(len(x_pos))
        
        self.centroid = np.array([x_com, y_com])
        return

    
    def getAssemblyIndex(self):
        import subprocess
        import re
        clus.generateMolFile(self.clusterDataFrame)
        
        subprocess.run(["assemblyCpp_256.exe", "mol_file"])
        fName = "mol_fileOut"
        with open(fName, 'r') as file:
            content = file.read()
        numbers = re.findall(r'\d+', content)
        
        assembly_index = int(numbers[1])
        self.assemblyIndex = assembly_index


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

def persistence(centroid):
    return

def get_identifier(cluster):
    return


def main():
    from time import sleep
    plt.close('all')
    N = 300
    L = 7
    v = 0.03
    r = 1
    dt = 1
    steps = 40
    num_clusters = np.zeros(steps)
    objects = []

    eta = 0.1
    
    sim_vicsek = vm.VicsekModel(N, L, v, eta, r, dt)

    for n in range(steps):
            cluster_dict = {}

            sim_vicsek.update()
            pos = sim_vicsek.get_positions()
            angles = sim_vicsek.get_velocities()
            kinematic_df = clus.make_kinematic_df(pos, angles)
            
            clusterList = clus.makeClustersDbscan(kinematic_df)
            for ind in range(len(clusterList)):
                cluster_dict[ind] = cluster(clusterList[ind])
                cluster_dict[ind].computeCentroid()
                cluster_dict[ind].getAssemblyIndex()
            
            a = cluster_dict[0]
            a.computeCentroid()
            print(a.centroid)
            
            cluster_histogram(cluster_dict)
            
            
      #      print(num_clusters[n])
            
            # visualize_clusters(clusters, n)
       #     if n%10 == 0:
       #         cluster_histogram(clusters)
            
            # generate_mol_file(clusters)

            # assembly_indices_step_n = call_assembly_code(num_clusters)
            
            # diff = np.max(assembly_indices_step_n)\
            #     - np.min(assembly_indices_step_n)

            # plt.figure()
            # plt.hist(assembly_indices_step_n, diff, histtype='bar')
            # plt.xlabel(r"Assembly Index $a$")
            # plt.ylabel("Count")
            
            # assembly_indices_all_time.append(assembly_indices_step_n)
            
    
            
            
            

if __name__ == '__main__':
    main()

