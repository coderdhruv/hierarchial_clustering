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

def find_maximum(mat, n):
    maxi = 0
    h1 = 0
    h2 = 0
    for i in range(0, n):
        for j in range(0, n):
            if mat[i][j] != 9999:
                if mat[i][j] > maxi and i < j:
                    if set(groups[i]).issubset(set(groups[j])) == 0:
                        maxi = mat[i][j]
                        h1 = i
                        h2 = j

    return h1, h2, maxi


def max_dist_calculate(y, z, matrix):
    maxim = 0
    for i in range(0, len(y)):
        for j in range(0, len(z)):
            if matrix[y[i]][z[j]] > maxim and y[i] != z[j]:
                maxim = matrix[y[i]][z[j]]
    return maxim



def augmented_dendrogram(*args, **kwargs):
        data = scipy.cluster.hierarchy.dendrogram(*args, **kwargs)
        if not kwargs.get('no_plot', False):
            for i, d in zip(data['icoord'], data['dcoord']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                plt.plot(x, y, 'ro')
                plt.annotate("%.3g" % y, (x, y), xytext=(0,12),textcoords='offset points',va='top', ha='center')
        return data


s = input()
data = list()
d = dict()
with open(s, 'r') as f:
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
            s = s+string
    data.append((head,s,count))
    d[head] = s

dist1 = np.zeros(shape=(len(data), len(data)))
finalArray = np.zeros(shape=(len(data), len(data)))

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
        a, b, min_value = find_maximum(dist1, len(data))
        print("a:", a, "b:", b)
        list1.append(groups[a]);
        list2.append(groups[b])
        list_val.append(min_value)

        """if len(groups[a]) == 1 and len(groups[b]) == 1:
            tmpZ[ptr][0] = a
            tmpZ[ptr][1] = b
            tmpZ[ptr][2] = min_value
            ptr = ptr + 1
        elif len(groups[a]) != 1 and len(groups[b]) == 1:
            tmpZ[ptr][0] = g1
            tmpZ[ptr][1] = b
            tmpZ[ptr][2] = min_value
            ptr = ptr + 1
            g1 = g1 + 1
        elif len(groups[a]) == 1 and len(groups[b]) != 1:
            print("value g", g1)
            tmpZ[ptr][0] = a
            tmpZ[ptr][1] = g1
            tmpZ[ptr][2] = min_value
            ptr = ptr + 1
            g1 = g1 + 1
        else:
            tmpZ[ptr][0] = g1
            tmpZ[ptr][1] = g1+1
            tmpZ[ptr][2] = min_value
            ptr = ptr + 1
            g1 = g1 + 2"""
        clusters = list()
        clusters = list(set(groups[a]) | set(groups[b]))
        groups[a] = clusters
        print(groups[a], "printing merging sets", groups[b])
        print(groups)
        for j in range(0, len(data)):
            dist1[j][groups[b][groups[b].index(min(groups[b]))]] = 9999
            dist1[groups[b][groups[b].index(min(groups[b]))]][j] = 9999
        for j in range(0, len(data)):
            for k in range(0, len(data)):
                if dist1[j][k] != 9999 and j != k:
                    dist1[j][k] = max_dist_calculate(groups[j], groups[k], finalArray)
        flag = 0
        print(tmpZ)
        print(dist1, "distane matrix")
    else:
        break
n = len(data)
print(list1,"printing 1")
print(list2,"printing 2")
print(list_val,"printing distance")
for i in range(0, len(list1)):
    for j in range(i+1, len(list1)):
        sum_i = len(list1[i]) + len(list2[i])
        sum_j = len(list1[j]) + len(list2[j])
        if sum_j < sum_i:
            list1[i], list2[i], list1[j], list2[j], list_val[i], list_val[j] = list1[j], list2[j], list1[i], list2[i], list_val[j], list_val[i]


print("list1=\n", list1, "list2=\n", list2, "list_val=\n", list_val)

final = []
for i in range(0, len(list1)):
    if len(list1[i]) == 1 and len(list2[i]) == 1:
        t = []
        t.append(list1[i][0])
        t.append(list2[i][0])
        t.append(list_val[i])
        t.append(2)
        final.append(t)
    elif len(list1[i]) == 1 and len(list2[i]) != 1:
        t = []
        t.append(list1[i][0])
        list2[i].sort()
        for j in range(0, i):
            if len(list1[j]) + len(list2[j]) == len(list2[i]):
                mixed = []
                for k in range(0, len(list1[j])):
                    mixed.append(list1[j][k])
                for k in range(0, len(list2[j])):
                    mixed.append(list2[j][k])
                mixed.sort()
                if mixed == list2[i]:
                    t.append(n+j)
                    t.append(list_val[i])
                    t.append(len(list1[i]) + len(list2[i]))
                    break
        final.append(t)
    elif len(list1[i]) != 1 and len(list2[i]) == 1:
        t = []
        t.append(list2[i][0])
        list1[i].sort()
        for j in range(0, i):
            if len(list1[j]) + len(list2[j]) == len(list1[i]):
                mixed = []
                for k in range(0, len(list1[j])):
                    mixed.append(list1[j][k])
                for k in range(0, len(list2[j])):
                    mixed.append(list2[j][k])
                mixed.sort()
                if mixed == list1[i]:
                    t.append(n + j)
                    t.append(list_val[i])
                    t.append(len(list1[i]) + len(list2[i]))
                    break
        final.append(t)
    else:
        t = []
        list1[i].sort()
        list2[i].sort()
        for j in range(0, i):
            if len(list1[j]) + len(list2[j]) == len(list1[i]):
                mixed = []
                for k in range(0, len(list1[j])):
                    mixed.append(list1[j][k])
                for k in range(0, len(list2[j])):
                    mixed.append(list2[j][k])
                mixed.sort()
                if mixed == list1[i]:
                    t.append(n + j)

        for j in range(0, i):
            if len(list1[j]) + len(list2[j]) == len(list2[i]):
                mixed = []
                for k in range(0, len(list1[j])):
                    mixed.append(list1[j][k])
                for k in range(0, len(list2[j])):
                    mixed.append(list2[j][k])
                mixed.sort()
                if mixed == list2[i]:
                    t.append(n + j)
        t.append(list_val[i])
        t.append(len(list1[i]) + len(list2[i]))
        final.append(t)
print("printing final=",final)
final_numpy = np.zeros(shape=(len(final), 4))
for i in range(0, len(final)):
    for j in range(0, 4):
        final_numpy[i][j] = final[i][j]
        print("i=", i, "j=", j,"val=", final[i][j])
print(list1, "printing 1 again")
print(list2, "printing 2 again")
print(list_val, "printing distance again")

# Plot dendrogram
names = [i for i in range(0, len(data))]
plt.figure(figsize=(25, 25))
plt.title('Hierarchical Clustering Dendrogram (Agglomerative)')
plt.xlabel('Sequence No.')
plt.ylabel('Distance')
augmented_dendrogram(final_numpy, labels=names, show_leaf_counts=True, p=70, truncate_mode = 'lastp')
plt.show()