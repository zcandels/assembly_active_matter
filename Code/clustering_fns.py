
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

def make_kinematic_df(pos, orientations):
    '''
    DESCRIPTION
        Creates Pandas dataframe from positions and velocity
        direction of particles (passed as numpy arrays)
    
    Parameters
    -----------
    pos: array-like of shape (N,2) containing the (x,y) positions
    of N particles
    
    orientations: array-like of shape N containing the velocity
    direction of each of the particles
    
    Returns
    --------
    position_orientation_df : pandas dataframe
    position and orientation data put into a pandas dataframe.
    '''
    kinematic_dict = {'x-position': pos[:,0], 'y-position': pos[:,1],
                'orientations':orientations}
    position_orientation_df = pd.DataFrame(data=kinematic_dict)
    
    return position_orientation_df


def make_clusters(kinematic_df):
    '''
    Creates clusters of data based on the K-means clustering algorithm
    implemented in sklearn.
    
    Parameters
    ----------
    kinematic_df: pandas dataframe whose columns contain the (x,y) 
    position of N particles and the velocity direction of the particles.
    
    Returns
    ---------
    data_1,..., data_n: pandas dataframes containing where the entries
    in dataframe data_i correspond to the particles in cluster i.
    '''
    wcss = []
    
    import matplotlib.pyplot as plt
    import numpy as np

    sillhouette_scores = []
    indexes = list(range(2,20))
    for i in range(2,20):
        k_means = KMeans(n_clusters=i, init='k-means++',
                         random_state=42)
        y = k_means.fit_predict(kinematic_df)
        wcss.append(k_means.inertia_)

        a = silhouette_score(kinematic_df, y)
        sillhouette_scores.append(silhouette_score(kinematic_df, y))
    
    optimal_cluster_ind = np.argmax(sillhouette_scores)
    optimal_cluster_num = indexes[optimal_cluster_ind]
    k_means_optimum = KMeans(optimal_cluster_num, init='k-means++',
                     random_state=42)
    y = k_means_optimum.fit_predict(kinematic_df)
        
    # plt.figure()
    # plt.plot(np.arange(2,20),wcss)
    # plt.xlabel('Clusters')
    # plt.ylabel('SSE')
    # plt.show()

    kinematic_df['cluster'] = y
    clusters = []
    for i in range(optimal_cluster_num):
        clusters.append(kinematic_df[kinematic_df.cluster==i])
    
    return clusters