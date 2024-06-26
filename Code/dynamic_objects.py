import numpy as np
import matplotlib.pyplot as plt
import clustering_fns as clus
import vic as vm


def get_edges_vertices(cluster):
    
    cluster_positions = cluster[ ['x-position', 'y-position'] ]
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


def generate_mol_file(clusters):
    
    
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
    
    for iter in range(len(clusters)):
        edges, vertices = get_edges_vertices(clusters[iter])
        bonds_info \
            = [(list(bond)[0],list(bond)[1], 1.0) for j,bond in enumerate(edges)]
        atom_list = [ (j,"C") for j,i in enumerate(vertices)]
        mol= tables2mol((atom_list,bonds_info))
        fileName = f"mol_file_cluster_{iter}"
        print(Chem.MolToMolBlock(mol),file=open(fileName +".mol",'w+'))

def call_assembly_code(num_clusters):
    import subprocess
    import re
    
    assembly_indices_single_timestep = []
    for n in range(num_clusters):
        subprocess.run(["assemblyCpp_256.exe", f"mol_file_cluster_{n}"])
        
        fName = f"mol_file_cluster_{n}Out"
        with open(fName, 'r') as file:
            content = file.read()
        numbers = re.findall(r'\d+', content)
        
        assembly_index = int(numbers[1])
        assembly_indices_single_timestep.append(assembly_index)
        print(assembly_index,",")
        
    
    return assembly_indices_single_timestep

def visualize_clusters(clusters, step):
    
    import matplotlib.cm as cm
    N_clus = len(clusters)
    start = 0.0
    stop = 1.0
    cm_subsection = np.linspace(start, stop, N_clus)
    colors = [ cm.jet(x) for x in cm_subsection ]
    
    
    plt.rcParams['text.usetex'] = True
    plt.figure()
    for n in range(N_clus):
        cluster = clusters[n]
        x = cluster['x-position']
        y = cluster['y-position']
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
    
def cluster_histogram(clusters):
    particle_distribution = []
    for cluster in clusters:
        num_particles = cluster.shape[0]
        particle_distribution.append(num_particles)
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


def main():
    from time import sleep
    plt.close('all')
    N = 300
    L = 7
    v = 0.03
    r = 1
    dt = 1
    steps = 40
    assembly_indices_all_time = []
    num_clusters = np.zeros(steps)

    eta = 0.1
    
    sim_vicsek = vm.VicsekModel(N, L, v, eta, r, dt)

    for n in range(steps):

            sim_vicsek.update()
            pos = sim_vicsek.get_positions()
            angles = sim_vicsek.get_velocities()
            kinematic_df = clus.make_kinematic_df(pos, angles)
            
            clusters = clus.makeClustersDbscan(kinematic_df)
            num_clusters[n] = len(clusters)
            
            
            print(num_clusters[n])
            
            # visualize_clusters(clusters, n)
            if n%10 == 0:
                cluster_histogram(clusters)
            
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

