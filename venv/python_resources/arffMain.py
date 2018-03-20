# This file uses the liac-arff package authored by Renato de Pontes Pereira and Matthias Feurer
# Authors: Conor Lorsung and Kyle Becker

import arff
import math
import random

# Choose the distance measure to run with:
#    man = manhattan distance
#    euc = euclidean distance
distMeasure = 'euc'

# Choose the K value to run with:
k = 4

# Choose which file to run on:
#    file = 0 -> SigGene.arff
#    file = 1 -> AllGenes.arff
file = 0

dataset = []

if file == 0:
    dataset = arff.load(open("../data_sources/SigGene.arff", "r"))
elif file == 1:
    dataset = arff.load(open("../data_sources/AllGenes.arff", "r"))

datalen = len(dataset['data'])-1
attrlen = len(dataset['attributes'])-1

classes = []

recentroids = []

finalClusters = []

testCluster = [(0, 34, 10872.81252758457), (0, 35, 15032.566296544312), (0, 36, 9734.507098975273), (0, 59, 18802.87375004151), (0, 60, 16246.020191419188), (0, 63, 16877.212040500053), (0, 64, 15164.62049838373), (0, 65, 13485.726434271162), (0, 66, 15712.412579231745), (0, 72, 15978.173057643358)]

# Sample calls
# To get data from first attribute
#   data['data'][0]


# Manhattan distance function between two vectors
def manhattan(vec1, vec2):
    dist = 0
    for x in range(0, len(vec1)-1):
        dist += abs(vec1[x] - vec2[x])
    return dist


# Euclidean distance function between two vectors
def euclidean(vec1, vec2):
    dist = 0
    for x in range(0, len(vec1)-1):
        dist += (vec1[x] - vec2[x])**2
    return math.sqrt(dist)


# K-means function for given data-set
def kmeans(data, k):

    centroids = []

    for i in range(0, 100):

        mandists = []
        eucdists = []
        mclusters = []
        eclusters = []

        for x in range(0, k):
            mclusters.append([])
            eclusters.append([])
        if len(centroids) == 0:
            for x in range(0, k):
                index = random.randint(0, len(data)-1)
                centroids.append(data[index])
            print('\n')
            for x in centroids:
               print('CENTROID: ',data.index(x)+1)
            i = 0
            for x in data:
                i += 1
                for y in centroids:
                    mandists.append((centroids.index(y), i, manhattan(y, x)))
                    eucdists.append((centroids.index(y), i, euclidean(y, x)))
        elif recentroids == centroids:
            print('\nCOMPLETED AFTER ', i, 'ITERATIONS\n')
            for x in recentroids:
                print('FINAL CENTROID: ', x)
            print('\n')
            for x in finalClusters:
                print('FINAL CLUSTERS: Length = ', len(x),' : ', x)
            return
        else:
            centroids = recentroids
            for x in centroids:
                print('NEW CEN: ', x)
            i = 0
            for x in data:
                i += 1
                for y in centroids:
                    mandists.append((centroids.index(y), i, manhattan(y, x)))
                    eucdists.append((centroids.index(y), i, euclidean(y, x)))

        # Determine which distance measure the user wants to run
        if distMeasure == 'man':
            # Separate data into K clusters
            for x in range(0, len(mandists),k):
                dists = []
                for i in range(0, k):
                    dists.append(mandists[x+i])
                val = min(dists, key=lambda t: t[2])
                mclusters[val[0]].append(val)

            # Sort the data within the clusters by gene number
            for x in mclusters:
                x.sort(key=lambda t: t[1])

            for x in mclusters:
                print('CLUSTER: ', x)

            recentroids = recompute(mclusters)
            finalClusters = mclusters
        elif distMeasure == 'euc':
            # Separate data into K clusters
            for x in range(0, len(eucdists), k):
                dists = []
                for i in range(0, k):
                    dists.append(eucdists[x + i])
                val = min(dists, key=lambda t: t[2])
                eclusters[val[0]].append(val)

            # Sort the data within the clusters by gene number
            for x in eclusters:
                x.sort(key=lambda t: t[1])

            for x in eclusters:
                print('CLUSTER: ', x)

            recentroids = recompute(eclusters)
            finalClusters = eclusters


# Recompute the centroids for the given data-set
def recompute(centroids):
    print('LEN G1: ', len(centroids[0]), '\nLEN G2: ', len(centroids[1]))
    newcentroids = []
    for x in range(0, len(centroids)):
        newcentroids.append(dataset['data'][centroids[x][0][1]-1][:-1])
        for y in range(1, len(centroids[x])):
            opvec = dataset['data'][centroids[x][y][1]-1]
            for n in range(0, len(opvec)-1):
                newcentroids[x][n] += opvec[n]
    for x in newcentroids:
        print(x)

    for x in range(0, len(newcentroids)):
        for y in range(0, len(newcentroids[x])):
            newcentroids[x][y] = newcentroids[x][y] / len(centroids[x])

    for x in newcentroids:
        print(x)

    return newcentroids



# Get the class values for the given data-set
def classvals(data):
    classlabels = []
    classets = data['attributes'][attrlen][1]
    for x in classets:
        classlabels.append(x)
    return classlabels


# Entropy purity measure
def entropy(cluster):
    ent = 0
    totsize = 0
    ccounts = [0]*len(classes)
    for x in cluster:
        totsize += 1
        ccounts[classes.index(getclass(x))] += 1
    for x in ccounts:
        ent += -(x / totsize) * math.log((x / totsize), 2)
    return ent


# Retrieve the class label for the specified gene
def getclass(gene):
    return dataset['data'][gene[1]-1][attrlen]


def returnVals():
    return dataset


print(manhattan(dataset['data'][0], dataset['data'][1]))

print(euclidean(dataset['data'][0], dataset['data'][1]))

classes = classvals(dataset)

print(len(dataset['data']))

kmeans(dataset['data'], k)

print(entropy(testCluster))