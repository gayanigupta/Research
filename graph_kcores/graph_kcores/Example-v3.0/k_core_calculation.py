import networkx as nx
import matplotlib.pyplot as plt
G = nx.Graph()
G.add_nodes_from([1,2,3,4,5,6,7,8])
G.add_edges_from([(0,1),(0,2),
                  (1,0),(1,2),
                  (2,0),(2,1),(2,5),(2,4),(2,3),(2,6),
                  (3,4),(3,2),(3,6),(3,7),
                  (4,2),(4,6),(4,3),(4,7),
                  (5,1),(5,2),(5,6),(5,8),
                  (6,2),(6,4),(6,3),(6,5),(6,7),(6,8),
                  (7,4),(7,3),(7,6),
                  (8,5),(8,6)
                ])
nx.draw(G, with_labels=True, pos=nx.spring_layout(G))
plt.show()
#print(G.nodes())
#print(G.edges())
#print(nx.k_core(G,k=3).edges())
#graphs = nx.connected_component_subgraphs(nx.k_core(G,k=2))
G = nx.k_core(G, k=3)
nx.draw(G, with_labels=True, pos=nx.spring_layout(G))
plt.show()