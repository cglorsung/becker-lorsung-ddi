# Get delimiting character from the user.
delimiter = raw_input("What is the delimiting value in this data?")

# This list contains all of the numeric gene values along
# with their class label in the last index.
geneData = []

# This list contains all of the distances between one vector and the rest
distData = []

# No use yet.
dataCollection = []

# Open the Significant Gene ARFF file and read all gene data into lists.
with open("../data_sources/SigGene.arff", "r") as dataFile:
    for x in dataFile.readlines():
        if x[0].isdigit():
            x = x.split(delimiter)
            for j in x:
                if j[0].isdigit():
                    j = float(j)
            geneData.append(x)
dataFile.close()

# Print all values that have been read and parsed.
#for x in geneData:
#    print x


def manhattan(dist1, dist2):
    totDist = 0
    for i in range(0, len(dist1)-1):
        totDist += abs(float(dist1[i])-float(dist2[i]))
    return totDist


for i in range(1, len(geneData)-1):
    distData.append(manhattan(geneData[0], geneData[i]))
    #print manhattan(geneData[0], geneData[i])

print distData

