import numpy as np
import matplotlib.pyplot as plt
import clustering_fns as clus
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import vic
import visualization_module as vm
import cluster_dyn_objs as cdo
import networkx as nx
import os
import sys    


def getObjectKey(object_dict, centroid, epsilon):
    for obj in object_dict:
        obj_centroid = object_dict[obj].centroidPosition[-1]
        if np.linalg.norm(obj_centroid - centroid) < epsilon:
            key = obj
    return key

def jas(obj1, obj2, OS):
    import subprocess
    import re
    
    graph1, graph2 = nx.Graph(obj1.adj_mat), nx.Graph(obj2.adj_mat)
    
    composite_graph = nx.compose(graph1, graph2)
    comp_adj_mat = nx.adjacency_matrix(composite_graph)
    edges, vertices = clus.adj_mat_to_edgeVert(comp_adj_mat)
    
    clus.generateMolFile(comp_adj_mat)
    
    pwd = os.getcwd()
    
    if OS == "win": 
        exe_path = pwd + "/assemblyCpp_256.exe"
        subprocess.run([exe_path, "mol_file"],
                          stdout = subprocess.DEVNULL)
    elif OS == "nix":
        exe_path = pwd + "/assemblyCpp"
        subprocess.run([exe_path, "mol_file"],
                          stdout = subprocess.DEVNULL)
        
    fName = "mol_fileOut"
    with open(fName, 'r') as file:
        content = file.read()
        numbers = re.findall(r'\d+', content)
    
    assembly_index = int(numbers[0])
    return assembly_index
    

def assign_copy_nums(object_dict, OS):
    unique_obj_iter = 0
    for _, o1  in object_dict.items():
        if o1.DoA == "Alive":
            if o1.label == "default":
                for _, o2 in object_dict.items():
                    if o1.assemblyIndex[-1] == o2.assemblyIndex[-1]:
                        mu = o1.assemblyIndex[-1] + 1
                        if mu == jas(o1, o2, OS):
                            print("multiple objects!")
                            o1.copyNum += 1
                            o2.copyNum += 1
                            o1.label = unique_obj_iter
                            o2.label = unique_obj_iter
                            

def computeFullAssembly(object_dict, fullAssembly, step):
    accounted_for = []
    for _, obj in object_dict.items():
        label = obj.label
        if label != "default" and label in accounted_for:
            continue
        accounted_for.append(label)
        assInd = obj.assemblyIndex[-1]
        copyNum = obj.copyNum
        fullAssembly[step] += np.exp(assInd)*(copyNum - 1)
    fullAssembly[step] /= len(object_dict)
        
        
                        


def do_timesteps(steps, sim_vicsek, epsilon, OS):
    dynObjectId = 0
    object_dict = {}
    DoA_dict = {}
    
    assemblyIndices_mean_std = np.zeros( (steps, 2) )
    assembly_indices_current_step = []
    
    fullAssembly = np.zeros(steps)
            
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
            cluster_dict[ind] = cdo.cluster(clusterList[ind], OS)
            cluster_dict[ind].computeCentroid()
            cluster_dict[ind].get_adj_mat()
            cluster_dict[ind].getAssemblyIndex()
            assembly_indices_current_step.append(cluster_dict[ind].assemblyIndex)
            #print(cluster_dict[ind].assemblyIndex)
            
            if n == 2:
                centroid = cluster_dict[ind].centroid
                directed_dist_vec = centroids_nmp1 - centroid
                distances = np.linalg.norm(directed_dist_vec, axis=1 )
                for d in distances:
                    if d < epsilon:
                        dynObjectId += 1
                        numParticles = cluster_dict[ind].numParticles
                        adj_mat = cluster_dict[ind].adj_mat
                        assemblyIndex = cluster_dict[ind].assemblyIndex
                        area = 1
                        # add some statement to see if the condition
                        # is satisfied for multiple clusters
                        key = dynObjectId
                        object_dict[dynObjectId] = cdo.dynamic_object(
                            centroid, adj_mat, numParticles,
                            assemblyIndex, area)
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
                                adj_mat = cluster_dict[ind].adj_mat
                                aInd = cluster_dict[ind].assemblyIndex
                                numParticles = cluster_dict[ind].numParticles
                                area = 1. 

                                dynObjectId = getObjectKey(object_dict,
                                                   centroid,
                                                   epsilon)
                                
                                object_dict[dynObjectId].updateObject(centroid,
                                                         adj_mat,
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
                        adj_mat = cluster_dict[ind].adj_mat
                        area = 1
                        object_dict[dynObjectId] = cdo.dynamic_object(
                            centroid, adj_mat,
                            numParticles, assemblyIndex, area)
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

        centroids_nmp1 = []
        for _, clusterObject in cluster_dict.items():
            centroids_nmp1.append(clusterObject.centroid)
        centroids_nmp1 = np.asarray(centroids_nmp1)
        
        if n >= 2:
            for key in object_dict:
                DoA_dict[key] = 0
            if n%10 == 0:  
                vm.object_histogram(object_dict, sim_vicsek)
                
#            assign_copy_nums(object_dict, OS)
#            computeFullAssembly(object_dict, fullAssembly, n)
            #vm.visualize_clusters(cluster_dict, n)
                
        assemblyIndices_mean_std[n,0] = np.mean(assembly_indices_current_step)
        assemblyIndices_mean_std[n,1] = np.std(assembly_indices_current_step)
        

        
        
    fName = "meanStdAssemblyIndexOverTime.dat"
    np.savetxt(fName, assemblyIndices_mean_std)
    
    f2Name = "fullAssemblyOverTime.dat"
    np.savetxt(f2Name, fullAssembly)
    
    return assemblyIndices_mean_std, fullAssembly
         


def main(rho, eta):
    import time
    
    tic = time.time()
    plt.close('all')
    #rho = float(sys.argv[1]) #97
    L = 10
    N = round(rho*L**2)
    v = 0.03
    r = 1
    dt = 1
    steps = 4
    
    epsilon = 3*v

    #eta = float(sys.argv[1]) # 0.001
    
    OS = "nix"
    sim_vicsek = vic.VicsekModel(N, L, v, eta, r, dt)
    assemblyIndices_mean_std, fullAssembly\
            = do_timesteps(steps, sim_vicsek, epsilon, OS)
    
        #vm.mean_assembly_ind(assemblyIndices_mean_std, sim_vicsek)
        #vm.full_assembly(fullAssembly, sim_vicsek)
    
    toc = time.time()
    
    print(toc - tic)
            
            

if __name__ == '__main__':
    rho = float(sys.argv[1])
    eta = float(sys.argv[2])
    main(rho, eta)

