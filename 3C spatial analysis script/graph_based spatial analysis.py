"""
Date: 18/11/2020
@author: Ali Hashmi
"""
from scipy.spatial import Delaunay
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pylab
from mayavi import mlab
from functools import reduce
options = {'node_color':'#1f78b4','with_labels': False,'arrows':True,'node_size': 10,'width': 0.5,'arrowstyle': '-|>','arrowsize': 12};

# class MeshGraph generates an undirected graph from the Delaunay mesh of xyz coordinates
class MeshGraph:
    def __init__(self,assoc): #fin
        self.assoc = assoc
        self.labels,self.pts = list(self.assoc.keys()),tuple(map(tuple,self.assoc.values()))
        self.revassoc = dict(zip(self.pts,self.labels))
        self.pts = np.asarray(self.pts)
        def BuildGraph():
            tetrahed = Delaunay(self.pts)
            arr=tetrahed.simplices
            ls = []
            # extracting all edge connections
            for pt_ind in range(len(self.pts)):
                log=tetrahed.simplices == pt_ind
                ind=np.any(log,axis=1)
                part=np.unique(arr[ind].flatten())
                tuples = [(x,pt_ind) for x in part if x != pt_ind]
                ls.append(tuples)
            # creating graph/network
            G = nx.Graph()
            for arr in ls:
                tk = np.take(self.pts,arr,axis=0)
                for i in tk:
                    x,y = tuple(map(tuple,i))
                    G.add_edge(x,y)
            return G

        self.graph = BuildGraph() # this gives the graph

    def to_undirected(self): #fin
        H = self.graph
        self.graph = H.to_undirected()

    def is_directed(self): #fin
        return not self.graph.is_directed()

    def is_empty(self): #fin
        return nx.is_empty(self.graph)

    def getNeighboursOfNode(self,node):
        nodesel = tuple(self.assoc[node])
        return list(map(lambda x: self.revassoc[x],self.graph.neighbors(nodesel)))

    def getAllNeighbours(self):
        return list( map(lambda x : list(map(lambda y: self.revassoc[y], self.graph.neighbors(tuple(self.assoc[x])))),np.sort(self.labels)) )

    def degree_dist(self): #fin
        l=nx.degree_histogram(self.graph)
        return list(zip(range(len(l)),l))

    def __PruneGraph__(self,criteria):
        # method to prune the graph edges based on euclidean distances between nuclei
        edges = list(self.graph.edges())
        print(len(edges))
        edgedist=np.array(list(map(lambda x : np.linalg.norm(np.array(x[0])-np.array(x[1])), edges)))
        bool = list(edgedist >= criteria)
        boolinds=[i for i, x in enumerate(bool) if x]
        print(type(edges))
        print("pruning", len(boolinds), "edges")
        for j in boolinds:
            e = edges[j]
            self.graph.remove_edge(*e)
        return None

    def __getAdjacencyMatrix__(self): #fin
        return nx.adjacency_matrix(self.graph).todense()

    def nodes(self): #fin
        return list(map(lambda x: self.revassoc[x],self.graph.nodes()))

    def edges(self): #fin
        #return list(self.graph.edges())
        return list(map(lambda x: (self.revassoc[x[0]],self.revassoc[x[1]]),self.graph.edges()))

    def info(self):  #fin
        print(nx.info(self.graph) + "\ndensity: ",nx.density(self.graph))

    @staticmethod
    def print3Dembeddingin2D(graph,options):
         nx.draw_networkx(graph,pos=nx.spring_layout(graph),**options)

    @staticmethod
    def print3dGraph(graph):
        H = nx.convert_node_labels_to_integers(graph)
        pos = nx.spring_layout(H, dim=3)
        xyz = np.array([pos[v] for v in sorted(H)])
        cols = np.array(list(H.nodes()))
        points = mlab.points3d(xyz[:, 0],xyz[:, 1],xyz[:, 2],cols,scale_factor=0.1,scale_mode="none",colormap="Blues",resolution=20);
        points.mlab_source.dataset.lines = np.array(list(H.edges()))
        tubes = mlab.pipeline.tube(points, tube_radius=0.01)
        mlab.pipeline.surface(tubes, color=(0.8, 0.8, 0.8))
        mlab.show()

    def prop(self):
        print({"num_of_nodes": self.graph.number_of_nodes(),"num_of_edges": self.graph.number_of_edges()})


########################################################################################################################################
# creating synthetic xyz coordinates
points=np.random.rand(500,3);
labels = range(len(points))
assoc = dict(zip(labels,points))

#buildgraph
graph = MeshGraph(assoc)
graph.to_undirected()
graph.is_empty()

# basic check, points[graph_nodes] and graph.assoc[graph_nodes] -- not relevant for users
#np.all(list(map(lambda x: np.all(points[x]==graph.assoc[x]), graph.nodes())))

#print functions
MeshGraph.print3Dembeddingin2D(graph.graph,options)
MeshGraph.print3dGraph(graph.graph)

# get adjacency matrix
graph.__getAdjacencyMatrix__()

# basic graph properties
graph.edges()
graph.nodes()
graph.prop()
graph.is_directed()
graph.info()
graph.prop()
graph.degree_dist()

# neighbour of a given node (0 ...n-1 )
graph.getNeighboursOfNode(0)
graph.getAllNeighbours()

# prune edges of graph if needed *(connections between nuclei that are more than the threshold distance apart are severed)
graph.__PruneGraph__(0.20)
graph.prop()

MeshGraph.print3Dembeddingin2D(graph.graph,options)
MeshGraph.print3dGraph(graph.graph)

# :-> similar methods can be copied to the class for adding and removing nodes/edges
#graph.add_edge()
#graph.add_edges_from()
#graph.remove_edge()
#graph.remove_edges_from()
#graph.add_node()
#graph.add_nodes_from()
#graph.remove_node()
#graph.remove_nodes_from()
