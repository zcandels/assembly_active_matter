import numpy as np
import matplotlib.pyplot as plt
import clustering_fns as clus
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import vic
import visualization_module as vm


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
        self.centroidPosition = [centroidPosition]
        self.numParticles = [numParticles]
        self.assemblyIndex = [assemblyIndex]
        self.area = [area]
        self.lifeTime = 1
        self.DoA = "Alive"
        # self.data_dict = {"CoM": [centroidPosition],
        #                    "numParticles": [numParticles],
        #                    "assemblyIndex": [assemblyIndex],
        #                    "area": [area], "lifeTime": 1}
    
    def updateObject(self, centroidPosition,
                         numParticles, assemblyIndex, area):
        self.centroidPosition.append(centroidPosition)
        self.numParticles.append(numParticles)
        self.assemblyIndex.append(assemblyIndex)
        self.area.append(area)
        self.lifeTime += 1
        self.DoA = "Alive"
        #objDict = self.data_dict
        #objDict["CoM"].append(centroidPosition)
        #objDict["numParticles"].append(numParticles)
        #objDict["assemblyIndex"].append(assemblyIndex)
        #objDict["area"].append(area)
        #objDict["lifeTime"] += 1
        return
    
    def deathCertificate(self):
        self.DoA = "Dead"
    
    
def getObjectKey(object_dict, centroid, epsilon):
    for obj in object_dict:
        obj_centroid = object_dict[obj].centroidPosition[-1]
        if np.linalg.norm(obj_centroid - centroid) < epsilon:
            key = obj
    return key



def do_timesteps(steps, sim_vicsek, epsilon):
    dynObjectId = 0
    object_dict = {}
    DoA_dict = {}
            
    centroids_nmp1 = []
    cluster_dict = {}
    for n in range(steps):

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
                        dynObjectId += 1
                        numParticles = cluster_dict[ind].numParticles
                        assemblyIndex = cluster_dict[ind].assemblyIndex
                        area = 1
                        # add some statement to see if the condition
                        # is satisfied for multiple clusters
                        key = dynObjectId
                        object_dict[dynObjectId] = dynamic_object(
                            centroid, numParticles, assemblyIndex, area)
                        DoA_dict[dynObjectId] = 1
            elif n > 2:
                centroid = cluster_dict[ind].centroid
                directed_dist_vec = centroids_nmp1 - centroid
                distances = np.linalg.norm( directed_dist_vec, axis=1 )
                for d in distances:
                    if d < epsilon:
                        distances = []
                        # Loop through all extant objects to see if 
                        # the current cluster already has a 
                        # corresponding object of type dynamic_object.
                        for _, val in object_dict.items():
                            obj_CoM = val.centroidPosition
                            if np.linalg.norm(obj_CoM - centroid) < epsilon:
                                aInd = cluster_dict[ind].assemblyIndex
                                numParticles = cluster_dict[ind].numParticles
                                area = 1. 

                                dynObjectId = getObjectKey(object_dict,
                                                   centroid,
                                                   epsilon)
                                
                                object_dict[dynObjectId].updateObject(centroid,
                                                         aInd, 
                                                         numParticles,
                                                         area)
                                DoA_dict[dynObjectId] = 1
                                
                
        ##############################################################
                        # If the object is new, create a new 
                        # instance of the class dynamic_object.
                        dynObjectId += 1 
                        numParticles = cluster_dict[ind].numParticles
                        assemblyIndex = cluster_dict[ind].assemblyIndex
                        area = 1
                        object_dict[dynObjectId] = dynamic_object(
                            centroid, numParticles, assemblyIndex, area)
                        DoA_dict[dynObjectId] = 1
                        # Might not need the above line since 
                        # We only update the DoA dictionary
                        # when an object is updated to show that it is
                        # still alive. By default, when an object dyn_obj
                        # is created, dyn_obj.DoA = "Alive"
                        
                       # statements only executed if n > 2.
                for dynObjectId, DoA in DoA_dict.items():
                    if DoA == 0:
                        object_dict[dynObjectId].deathCertificate()
        
        
        #vm.visualize_clusters(clusters, n)
        if n%10 == 0:
            vm.cluster_histogram(cluster_dict)

        centroids_nmp1 = []
        for _, clusterObject in cluster_dict.items():
            centroids_nmp1.append(clusterObject.centroid)
        centroids_nmp1 = np.asarray(centroids_nmp1)
        
        if n >= 2:
            for key in object_dict:
                DoA_dict[key] = 0
    return object_dict
         


def main():
    plt.close('all')
    N = 300
    L = 3.1
    v = 0.03
    r = 1
    dt = 1
    steps = 5
    
    
    num_clusters = np.zeros(steps, dtype=int)
    epsilon = 3*v

    eta = 0.1
    
    sim_vicsek = vic.VicsekModel(N, L, v, eta, r, dt)
    
    object_dict = do_timesteps(steps, sim_vicsek, epsilon)
            
            

if __name__ == '__main__':
    main()

