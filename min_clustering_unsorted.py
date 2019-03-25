from typing import List
import distance
import numpy as np
from collections import defaultdict
from itertools import combinations,chain,combinations_with_replacement
import scipy
from scipy.cluster import hierarchy
from scipy.cluster import setup
import plotly.plotly as py
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
import plotly
import scipy.cluster.hierarchy as sch

"""function to find the minimum element indexes 
among all the elements in the proximity matrix"""
def find_minimum(mat,n):
    mini = 100000
    h1 = 0
    h2 = 0
    for i in range(0,n):
        for j in range(0,n):
            if mat[i][j] != 9999:
                if mat[i][j] < mini and i < j:
                    mini = mat[i][j]
                    h1 = i
                    h2 = j
    return h1, h2, mini

"""function to find the find the minimum distance between two data points,(they could be in
combined form or in singular form)"""

def dist_calculate(y, z, matrix):
    minim = 100000
    for i in range(0,len(y)):
        for j in range(0,len(z)):
            if matrix[y[i]][z[j]] < minim and y[i] != z[j]:
                minim = matrix[y[i]][z[j]]
    return minim


"""function to draw the dendrogram"""

def augmented_dendrogram(*args, **kwargs):
        data = scipy.cluster.hierarchy.dendrogram(*args, **kwargs)
        if not kwargs.get('no_plot', False):
            for i, d in zip(data['icoord'], data['dcoord']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                plt.plot(x, y, 'ro')
                plt.annotate("%.3g" % y, (x, y), xytext=(0,12),textcoords='offset points',va='top', ha='center')
        return data

""" code to open and write data to a np array initialised with zeroes """

s = input()
data = list()
d = dict()
with open(s,'r') as f:
    s = ""
    head = ""
    count = 0
    for row in f:
        string = row[:-1]
        #print(row)
        if(string.startswith('>')):
            if(head != ""):
                data.append((head, s, count))
                d[head] = s
            head = string
            s = ""
            count +=  1
        else:
            s = s+string[:-1]
    data.append((head,s,count))
    d[head] = s

dist1 = np.zeros(shape=(len(data), len(data)))
finalArray = np.zeros(shape=(len(data), len(data)))

"""code to fill the abovementioned array with levenstein distance values
   and hence making our distance matrix"""
for i in range(0, len(data)):
    for j in range(0, len(data)):
        if i != j and i < j:
            dist1[i][j] = dist1[j][i] = distance.levenshtein(data[i][1], data[j][1])
            finalArray[i][j] = finalArray[j][i] = distance.levenshtein(data[i][1], data[j][1])

groups = {}
r = dict()
global_cluster = list()
tmpZ = np.zeros(shape=(len(data)-1, 4))
g1 = len(data)
print(g1, "length of g")
ptr = 0
list1 = []
list2 = []
list_val = []
for i in range(0, len(data)):
    cluster_id = list()
    cluster_id.append(i)
    global_cluster.append(i)
    groups[i] = cluster_id
#print(dist1)
for i in range(0,len(data)-1):
    for h in range(0,len(data)):
        for g in range(0,len(data)):
            if dist1[h][g] != 9999 :
                flag = 1
                break
    if flag == 1:
        a, b, min_value = find_minimum(dist1, len(data))
        print("a:", a, "b:", b)
        print("groups[a]:", groups[a], "groups[b]:", groups[b])
        """filling the linkage-matrix tmpZ"""

        if len(groups[a]) == 1 and len(groups[b]) == 1:
            tmpZ[ptr][0] = a
            tmpZ[ptr][1] = b
            tmpZ[ptr][2] = min_value
            tmpZ[ptr][3] = len(groups[a])+len(groups[b])
            ptr = ptr + 1
        elif len(groups[a]) != 1 and len(groups[b]) == 1:
            tmpZ[ptr][0] = g1
            tmpZ[ptr][1] = b
            tmpZ[ptr][2] = min_value
            tmpZ[ptr][3] = len(groups[a])+len(groups[b])
            ptr = ptr + 1
            g1 = g1 + 1
        elif len(groups[a]) == 1 and len(groups[b]) != 1:
            print("value g", g1)
            tmpZ[ptr][0] = a
            tmpZ[ptr][1] = g1
            tmpZ[ptr][2] = min_value
            tmpZ[ptr][3] = len(groups[a])+len(groups[b])
            ptr = ptr + 1
            g1 = g1 + 1
        else:
            tmpZ[ptr][0] = g1
            tmpZ[ptr][1] = g1 + 1
            tmpZ[ptr][2] = min_value
            tmpZ[ptr][3] = len(groups[a])+len(groups[b])
            ptr = ptr + 1
            g1 = g1 + 2
        clusters = list()
        clusters = list(set(groups[a]) | set(groups[b]))
        groups[a] = clusters
        print(groups[a], "printing merging sets", groups[b])
        print(groups)

        """updating the proximity matrix"""

        for j in range(0, len(data)):
            dist1[j][groups[b][groups[b].index(min(groups[b]))]] = 9999
            dist1[groups[b][groups[b].index(min(groups[b]))]][j] = 9999
        for j in range(0, len(data)):
            for k in range(0, len(data)):
                if dist1[j][k] != 9999 and j != k:
                    dist1[j][k] = dist_calculate(groups[j], groups[k], finalArray)
        flag = 0
        print(tmpZ)
        print(dist1, "distane matrix")
    else:
        break


# Plot dendrogram
names = [data[i][0] for i in range(0, len(data))]
plt.figure(figsize=(25, 25))
plt.title('Hierarchical Clustering Dendrogram (Agglomerative)')
plt.xlabel('Sequence No.')
plt.ylabel('Distance')
augmented_dendrogram(tmpZ, labels=names, show_leaf_counts=True, p=200, truncate_mode='lastp')
plt.show()