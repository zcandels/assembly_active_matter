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


def visualize_results(clusters):
    
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
    
    return
    


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
        fileName = f"assembly_index_cluster_{iter}"
        print(Chem.MolToMolBlock(mol),file=open(fileName +".mol",'w+'))


def main():
    plt.close('all')
    N = 40
    L = 3.1
    v = 0.03
    r = 1
    dt = 1
    steps = 5

    eta = 0.1
    
    sim_vicsek = vm.VicsekModel(N, L, v, eta, r, dt)

    for n in range(steps):
            sim_vicsek.update()
            pos = sim_vicsek.get_positions()
            angles = sim_vicsek.get_velocities()
            kinematic_df = clus.make_kinematic_df(pos, angles)
            
            clusters = clus.make_clusters(kinematic_df)
            
            
            #vertices, edges = get_edges_vertices(clusters)
            
            generate_mol_file(clusters)
            
            visualize_results(clusters)
            

if __name__ == '__main__':
    main()

