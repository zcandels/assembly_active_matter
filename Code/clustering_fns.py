
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
from sklearn import preprocessing

def make_kinematic_df(pos, orientations):
    '''
    DESCRIPTION
        Creates Pandas dataframe from positions and velocity
        direction of particles (passed as numpy arrays)
    
    Parameters
    -----------
    pos: array-like of shape (N,2) containing the (x,y) positions
    of N particles
    
    angle: array-like of shape N containing the velocity
    direction of each of the particles
    
    Returns
    --------
    position_angle_df : pandas dataframe
    position and orientation data put into a pandas dataframe.
    '''
    kinematic_dict = {'x-position': pos[:,0], 'y-position': pos[:,1],
                'angle':orientations}
    position_angle_df = pd.DataFrame(data=kinematic_dict)
    
    dfCluster = position_angle_df.copy()
    scaler = preprocessing.StandardScaler()
    
    dfCluster[['x-position', 'y-position', 'angle']]\
        = scaler.fit_transform(dfCluster[['x-position', 'y-position', 'angle']])
    dfCluster = dfCluster[['x-position', 'y-position', 'angle']]
    
    return dfCluster


def makeClustersKmeans(dfCluster):
    '''
    Creates clusters of data based on the K-means clustering algorithm
    implemented in sklearn.
    
    Parameters
    ----------
    dfCluster: pandas dataframe whose columns contain the (x,y) 
    position of N particles and the velocity direction of the particles.
    
    Returns
    ---------
    data_1,..., data_n: pandas dataframes containing where the entries
    in dataframe data_i correspond to the particles in cluster i.
    '''
    wcss = []
    
    from numpy import argmax
    from numpy import arange

    sillhouette_scores = []
    indexes = list(range(2,20))
    for i in range(2,20):
        k_means = KMeans(n_clusters=i, init='k-means++',
                         random_state=42)
        y = k_means.fit_predict(dfCluster)
        wcss.append(k_means.inertia_)

        sillhouette_scores.append(silhouette_score(dfCluster, y))
    
    optimal_cluster_ind = argmax(sillhouette_scores)
    optimal_cluster_num = indexes[optimal_cluster_ind]
    k_means_optimum = KMeans(optimal_cluster_num, init='k-means++',
                     random_state=42)
    y = k_means_optimum.fit_predict(dfCluster)
        
    # plt.figure()
    # plt.plot(arange(2,20),wcss)
    # plt.xlabel('Clusters')
    # plt.ylabel('SSE')
    # plt.show()

    dfCluster['cluster'] = y
    clusters = []
    for i in range(optimal_cluster_num):
        clusters.append(dfCluster[dfCluster.cluster==i])
    
    return clusters

def makeClustersDbscan(dfCluster):
    """
    Creates clusters of data based on the DBSCAN clustering algorithm
    implemented in sklearn.

    Parameters
    ----------
    dfCluster: pandas dataframe whose columns contain the (x,y) 
    position of N particles and the velocity direction of the particles.

    Returns
    ---------
    clusters: List of pandas dataframes where the columns of clusters[i] 
    are 'x-positions', 'y-positions', and 'cluster', the latter 
    containing the cluster to label to which a given particle belongs.
    """
    from numpy import unique 
    clustering = DBSCAN(eps = 0.2, min_samples=5).fit(dfCluster)
    y = clustering.labels_
    dfCluster['cluster-label'] = y
    clusters =[]
    [distinct_vals, _] = unique(y, return_counts=True)
    for num in distinct_vals:
       clusters.append(dfCluster[dfCluster.loc[:,'cluster-label']==num])

    return clusters



def computeEdgesVertices(clusterDataFrame):
    cluster_positions = clusterDataFrame[ ['x-position', 'y-position'] ]
    from sklearn.neighbors import kneighbors_graph
    graph1 = kneighbors_graph(cluster_positions, 2, mode='connectivity',
                          include_self=True)


    adj_mat = graph1.toarray()
    vertices = list(range(len(adj_mat)))

    Edges = []
    for i in range(len(adj_mat)):
        for j in range(len(adj_mat)):
            if i == j:
                continue
            if adj_mat[i,j] != 0:
                Edges.append({i,j})

    edges = []
    for i in range(len(Edges)):
        if Edges[i] not in edges:
            edges.append(Edges[i])
        

    return edges, vertices

def misc_fn(clusterDataFrame):
    from rdkit import Chem
    from rdkit import RDLogger 
    from rdkit.Chem.rdchem import RWMol
    a = 1
    b = 2
    
    def transfrom_bond(bond):
        if bond == 1.0:
            return Chem.rdchem.BondType.SINGLE
        if bond == 2.0:
            return Chem.rdchem.BondType.DOUBLE
        if bond == 3.0:
            return Chem.rdchem.BondType.TRIPLE
        return "error"
    
    edges, vertices = computeEdgesVertices(clusterDataFrame)
    

def generateMolFile(clusterDataFrame):
    from rdkit import Chem
    from rdkit import RDLogger 
    from rdkit.Chem.rdchem import RWMol

    def transfrom_bond(bond):
        if bond == 1.0:
            return Chem.rdchem.BondType.SINGLE
        if bond == 2.0:
            return Chem.rdchem.BondType.DOUBLE
        if bond == 3.0:
            return Chem.rdchem.BondType.TRIPLE
        return "error"
    
    def transfrom_bond_float(bond):
        if bond == "single":
            return 1.0
        if bond == "double":
            return 2.0
        if bond == "triple":
            return 3.0
        return "error"
    
    def tables2mol(tables):
        atoms_info, bonds_info = tables
        emol = RWMol()
        for v in atoms_info:
            emol.AddAtom(Chem.Atom(v[1]))
        for e in bonds_info:
            emol.AddBond(e[0], e[1], transfrom_bond(e[2]))
        mol = emol.GetMol()
        return mol
    
    edges, vertices = computeEdgesVertices(clusterDataFrame)
    bonds_info \
        = [(list(bond)[0],list(bond)[1], 1.0) for j,bond in enumerate(edges)]
    atom_list = [ (j,"C") for j,i in enumerate(vertices)]
    mol= tables2mol((atom_list,bonds_info))
    fileName = f"mol_file"
    print(Chem.MolToMolBlock(mol),file=open(fileName +".mol",'w+'))
