
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
    from numpy import unique 
    clustering = DBSCAN(eps = 0.2, min_samples=2).fit(dfCluster)
    y = clustering.labels_
    dfCluster['cluster'] = y
    clusters =[]
    [distinct_vals, _] = unique(y, return_counts=True)
    counts = len(distinct_vals)
    for i in range(counts):
       clusters.append(dfCluster[dfCluster.cluster==i])

    return clusters

