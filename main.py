import random
import scipy
from scipy.optimize import Bounds
from scipy.optimize import minimize
import networkx as nx
import joblib
import networkx as nx
from networkx.algorithms import approximation as approx
import numpy as np
import math
import matplotlib.pyplot as plt


class positionObject:
    def __init__(self, position):
        self.pos = position


def weighti(dij, intNode, nodes, neighbors):
    sum = 0
    for i in neighbors:
        #print(i)
        if (dij[intNode, nodes.index(i)]) == 0 :
            print(i)
            print(nodes[intNode])
        sum = sum + 1/(dij[intNode, nodes.index(i)])
    return sum

def fullWeights(positionz, dij, nodes):
    sum = 0
    for i in nodes:
        for j in nodes:
            if not(i == j):
                if dij[nodes.index(j), nodes.index(i)] == 0:
                    continue
                distance = math.sqrt((positionz[i][0] - positionz[j][0])**2 + (positionz[i][1] - positionz[j][1])**2)
                weightedDistance = distance/dij[nodes.index(j), nodes.index(i)]
                weightedposition = math.log(1/weightedDistance, 10)
                sum = sum + (weightedposition)**2
    return sum


def getWeightSum(dij, intNode, nodes, weights):
    sum = 0
    weightCurrent = weights[intNode]
    for c in range(0,len(nodes)):
        if not(intNode == c):
            sum = sum + weightCurrent*weights[c]
    return sum


def fitnessEval(position,dij, positions, nodeInt, nodes, weights):
    #weightCurrent = weighti(dij, nodeInt, nodes)
    weightSum = getWeightSum(dij, nodeInt, nodes, weights)
    chisquared = 0
    c = 0
    sum = 0
    for c in range(0,len(nodes)):
        if not(c == nodeInt):
            distance = math.sqrt((positions[c][0] - position[0])**2 + (positions[c][1]- position[1])**2) * dij[nodeInt, c]
            sum = sum + ((weights[c]*weights[nodeInt])**(-0.5))*(1-2*(distance)**(0.5))
    chisquared = math.sqrt((sum/weightSum)**2)
    return chisquared

namelist=joblib.load('namelist')
dij=joblib.load('dijmatrix')
g=joblib.load('graph')
nodes = list(g.nodes())


######Add arbitrary values to distance matrix############
dijfull = np.zeros((1122,1122))
for rn in range(0,len(dij)):
    for cn in range(0,len(dij[0])):
       dijfull[rn,cn] = dij[rn,cn]


for n in nodes:
    sum = 0
    for i in dij[nodes.index(n)]:
        sum = sum + i
    average = sum/562
    for j in range(562,1122):
        dijfull[nodes.index(n), j] = average

#############################################


pos = nx.spring_layout(g)
protectedPos = positionObject(pos)
positions = {}


### Make a matrix containing the neighbors of each nodes ######
nodesNeighbors = {}
for j in nodes:
    nodesNeighbors[j] = []
    for i in g.neighbors(j):
        nodesNeighbors[j].append(i)
print(len(nodesNeighbors))
################################################


### Assign a weight to each node based on its connectivity #######
weights = {}
c = 0
for j in nodes:
    positions[c] = protectedPos.pos[nodes[c]]
    weights[c] = weighti(dijfull, c, nodes, nodesNeighbors)
    c = c + 1
print("Done with weight matrix")
###################################################


finalGraph = nx.Graph()
finalGraph = g

newPos = {}
newPos = pos


for z in range(0,5):
    for i in nodes:
        x0 = newPos[i]
        res = minimize(fitnessEval, x0, args=(dijfull, positions, nodes.index(i), nodes, weights),  method="Powell")
        print(i)
        print(res)
        print(" ")
        newPos[i] = res.x
        print(str(newPos[i]) + str(x0))
        positions[nodes.index(i)] = res.x


finalWeights = fullWeights(newPos, dijfull, nodes)
initalWeights = fullWeights(protectedPos.pos, dijfull, nodes)

diffpos = {}
for i in nodes:
    diffpos[i] = newPos[i] - protectedPos.pos[i]
    print(diffpos[i])


#nx.draw_networkx_nodes(finalGraph, pos = newPos, node_size=0.5)
nx.draw_networkx(finalGraph, pos = newPos, node_size=1, with_labels=False)
plt.show()
exit()



