import numpy as np
import matplotlib.pyplot as plt
import clustering_fns as clus
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import vic as vm


class cluster():
    def __init__(self, cluster_df):
        self.clusterDataFrame = cluster_df
        self.numParticles = cluster_df.shape[0]
        self.centroid = "default"
        self.assemblyIndex = "default"

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
        
        assembly_index = int(numbers[0])
        self.assemblyIndex = assembly_index


class dynamic_object():
    def __init__(self, centroidPosition, numParticles, assemblyIndex, area):
        self.centroidPosition = centroidPosition
        self.numParticles = numParticles
        self.assemblyIndex = assemblyIndex
        self.area = area
        self.lifeTime = 1
        self.data_dict = {"CoM": [centroidPosition],
                            "numParticles": [numParticles],
                            "assemblyIndex": [assemblyIndex],
                            "area": [area], "lifeTime": 1}
    
    def updateObject(self, centroidPosition,
                         numParticles, assemblyIndex, area):
        self.centroidPosition = centroidPosition
        self.numParticles = numParticles
        self.assemblyIndex = assemblyIndex
        self.area = area
        self.lifeTime += 1
        objDict = self.data_dict
        objDict["CoM"].append(centroidPosition)
        objDict["numParticles"].append(numParticles)
        objDict["assemblyIndex"].append(assemblyIndex)
        objDict["area"].append(area)
        objDict["lifeTime"] += 1
        return
    
def getObjectKey(object_dict, centroid, epsilon):
    for obj in object_dict:
        obj_centroid = object_dict[obj].data_dict["CoM"][-1]
        if np.linalg.norm(obj_centroid - centroid) < epsilon:
            key = obj
    return key


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
    plt.close('all')
    N = 300
    L = 3.1
    v = 0.03
    r = 1
    dt = 1
    steps = 5
    
    
    num_clusters = np.zeros(steps, dtype=int)
    object_dict = {}
    centroids_nmp1 = []
    epsilon = 3*v

    eta = 0.1
    
    sim_vicsek = vm.VicsekModel(N, L, v, eta, r, dt)
    
    newObjectCtr = 0
    for n in range(steps):
            cluster_dict = {}

            sim_vicsek.update()
            pos = sim_vicsek.get_positions()
            angles = sim_vicsek.get_velocities()
            kinematic_df = clus.make_kinematic_df(pos, angles)
            
            clusterList = clus.makeClustersDbscan(kinematic_df)
            #centroids = []
            
            for ind in range(1, len(clusterList)):
                # Change starting index of loop so we don't 
                # create a cluster object for straggler particles
                # i.e. particles for which their cluster label is -1.
                cluster_dict[ind] = cluster(clusterList[ind])
                cluster_dict[ind].computeCentroid()
                cluster_dict[ind].getAssemblyIndex()
                #print(cluster_dict[ind].assemblyIndex)
                
                if n == 2:
                    centroid = cluster_dict[ind].centroid
                    directed_dist_vec = centroids_nmp1 - centroid
                    distances = np.linalg.norm(directed_dist_vec, axis=1 )
                    for d in distances:
                        if d < epsilon:
                            newObjectCtr += 1
                            numParticles = cluster_dict[ind].numParticles
                            assemblyIndex = cluster_dict[ind].assemblyIndex
                            area = 1
                            # add some statement to see if the condition
                            # is satisfied for multiple clusters
                            key = f"obj{newObjectCtr}"
                            object_dict[key] = dynamic_object(
                                centroid, numParticles, assemblyIndex, area)
#                            print("n = ", n)
#                            print(object_dict[key].centroidPosition)
                elif n > 2:
                    centroid = cluster_dict[ind].centroid
                    directed_dist_vec = centroids_nmp1 - centroid
                    distances = np.linalg.norm( directed_dist_vec, axis=1 )
                    for d in distances:
                        if d < epsilon:
                            distances = []
                            for _, val in object_dict.items():
                                obj_CoM = val.centroidPosition
                                if np.linalg.norm(obj_CoM - centroid) < epsilon:
                                    aInd = cluster_dict[ind].assemblyIndex
                                    numParticles = cluster_dict[ind].numParticles
                                    area = 1. 

                                    key = getObjectKey(object_dict,
                                                       centroid,
                                                       epsilon)
                                    
                                    object_dict[str(key)].updateObject(centroid,
                                                             aInd, 
                                                             numParticles,
                                                             area)
            ##############################################################
                            newObjectCtr += 1 
                            numParticles = cluster_dict[ind].numParticles
                            assemblyIndex = cluster_dict[ind].assemblyIndex
                            area = 1
                            key = f"obj{newObjectCtr}"
                            object_dict[key] = dynamic_object(
                                centroid, numParticles, assemblyIndex, area)
                                    
                                
           ### I think this works but tomorrow check to make sure that when
           ### executing everything below line 213 that a new object is
           ### created and that you're not just reassigning new properties
           ### to a dictionary key that already exists.
            
            
            # visualize_clusters(clusters, n)
            if n%10 == 0:
                cluster_histogram(cluster_dict)

            centroids_nmp1 = []
            for _, clusterObject in cluster_dict.items():
                centroids_nmp1.append(clusterObject.centroid)
            centroids_nmp1 = np.asarray(centroids_nmp1)

            # diff = np.max(assembly_indices_step_n)\
            #     - np.min(assembly_indices_step_n)

            # plt.figure()
            # plt.hist(assembly_indices_step_n, diff, histtype='bar')
            # plt.xlabel(r"Assembly Index $a$")
            # plt.ylabel("Count")
            
            # assembly_indices_all_time.append(assembly_indices_step_n)
            

if __name__ == '__main__':
    main()

