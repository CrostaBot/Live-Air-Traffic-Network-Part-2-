import community
import networkx as nx
import matplotlib.pyplot as plt
import pprint as pp

G = nx.read_edgelist("dataset_3.txt",comments='#',create_using=nx.Graph())

partition = community.best_partition(G)
modularity = community.modularity(partition,G)
print(modularity)

size = len(set(partition.values()))
community = dict()

count = 0
for i in range(size):
    community[i] = count
    count = 0
    for nodes, comm in partition.items():
        if (comm == i):
            count += 1

pp.pprint(community)

values = [partition.get(node) for node in G.nodes()]
nx.draw_spring(G, cmap = plt.get_cmap('jet'), node_color = values, node_size=30, with_labels=False)
plt.show()
