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
            geneData.append(x.split(delimiter))
dataFile.close()

# Parse all numbers from strings to integers.
for x in geneData:
    for i in range(0, len(x)-1):
        if x[i][0].isdigit():
            x[i] = float(x[i])

# Print all values that have been read and parsed.
for x in geneData:
    print x