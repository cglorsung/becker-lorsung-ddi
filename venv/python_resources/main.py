# Get delimiting character from the user.
delimiter = raw_input("What is the delimiting value in this data?")

# This list contains all of the numeric gene values along
# with their class label in the last index.
geneData = []

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
for x in geneData:
     print x