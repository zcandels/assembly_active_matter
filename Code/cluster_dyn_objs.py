import numpy as np
import clustering_fns as clus

class cluster():
    def __init__(self, cluster_df):
        self.clusterDataFrame = cluster_df
        self.numParticles = cluster_df.shape[0]
        self.centroid = "default"
        self.assemblyIndex = "default"
        self.adj_mat = "default"

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
    
    
    def get_adj_mat(self):
        adj_mat = clus.make_adj_mat(self.clusterDataFrame)
        self.adj_mat = adj_mat
    
    def getAssemblyIndex(self):
        import subprocess
        import re
        clus.generateMolFile(self.adj_mat)
        
        subprocess.run(["assemblyCpp_256.exe", "mol_file"])
        fName = "mol_fileOut"
        with open(fName, 'r') as file:
            content = file.read()
        numbers = re.findall(r'\d+', content)
        
        assembly_index = int(numbers[0])
        self.assemblyIndex = assembly_index


class dynamic_object():
    def __init__(self, centroidPosition,
                 adj_mat, numParticles, assemblyIndex, area):
        self.centroidPosition = [centroidPosition]
        self.numParticles = [numParticles]
        self.adj_mat = adj_mat
        self.assemblyIndex = [assemblyIndex]
        self.area = [area]
        self.lifeTime = 1
        self.DoA = "Alive"
        # self.data_dict = {"CoM": [centroidPosition],
        #                    "numParticles": [numParticles],
        #                    "assemblyIndex": [assemblyIndex],
        #                    "area": [area], "lifeTime": 1}
    
    def updateObject(self, centroidPosition, adj_mat,
                         numParticles, assemblyIndex, area):
        self.centroidPosition.append(centroidPosition)
        self.numParticles.append(numParticles)
        self.assemblyIndex.append(assemblyIndex)
        self.area.append(area)
        self.adj_mat = adj_mat
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
    