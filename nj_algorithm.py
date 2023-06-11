from tree import Tree
import sys

# I ran the program on commandline with 'nj_algorithm.py [inputfile]'

# I had trouble getting the scikit-bio to work so I just made a basic Tree class
# a node in Tree has an id,parent,distance to parent, and a list of children
# the root has null value for parent, 0.0 for distance, and every node defaults to empty children list.

def systemInput():#systemInput checks if the function was called on commandline with a input file, if not it defaults to pre-selected file
    defaultFile = "example2.fna"
    try:
        if sys.argv[1] is None:
            return defaultFile
        else:
            return sys.argv[1]
    except IndexError as e:
        return defaultFile

#function that calculates difference between two sequences.
#returns the mismatches/(length of sequence) to get the percentage
def calculateDiff(seq1,seq2):
    mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
    #removes trailing 0 so that format matches example
    if mismatch == 0:
        return 0
    return mismatch / float(len(seq1))


# creates a matrix of differences using calculateDiff
def differenceMatrix(ids, idSeq):
    length = len(ids)
    matrix = [[0] * length for _ in range(length)]
    #fill matrix
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = calculateDiff(idSeq[id1], idSeq[id2])
    return matrix


# writes the difference matrix to file
def CreateDistanceFile(ids, matrix):
    with open('genetic_distances.txt', 'w') as f:
        #write the headers which are the ids
        f.write('\t' + '\t'.join(ids) + '\n')
        #fills in the first column as id and then the rest of the matrix
        for i in range(len(matrix)):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')

# reads the fna file and returns the ids and a dictionary of relationships
# for ids to sequence for fast reference
def readFile(filename):
    #the ids and their sequence
    dictionary = {}
    #the names of the sequence
    ids = []
    curr = None
    with open(filename) as f:
        content = f.read().splitlines()
    for i, line in enumerate(content):
        #sequences are found only on even lines
        if i % 2 == 0:
            curr = line[1:]
            ids.append(curr)
        else:
            dictionary[curr] = line
##    print("ids", ids)
##    print("dictionary", dictionary)
    return ids, dictionary

def neighbor_joining(inputIds, inputDistMatrix):
    # global root value and nodeId for use in recursion
    global root
    global nodeId
    
    # nodeId is based on the fna file
    nodeId = 120
    root = None
    #id name and its order of appearence
    idDictionary = {}
    #childParentDict stores the tree relationships as a dictionary
    #it is used to construct the tree later
    childParentDict = {}
    for i, id in enumerate(inputIds):
        idDictionary[id] = str(i + 1)
    print("idDictionary", idDictionary)

    #recursive function that calculates the Qmatrix while shrinking the distance matrix
    def recurse(ids, distMatrix):
        global root
        global nodeId
        N  = len(distMatrix)

        # base case if distance matrix is 2x2
        if N <= 2:
            
            # add the relationship into our dictionary and also set the root
            node1 = idDictionary[ids[0]] if ids[0] in idDictionary else ids[0]
            node2 = idDictionary[ids[1]] if ids[1] in idDictionary else ids[1]
            childParentDict[node1] = (node2, distMatrix[0][1])
            root = node2
            print("root: ", root)
            return

        # default q initialized to 0
        Q_matrix = [[0] * N for _ in range(N)]
        
        # variables to keep track of the minimum indices in the q matrix
        min_i = 0
        min_j = 0
        #set to inf so its larger than everything
        minVal = float('inf')

        #calculating the Q matrix
        #step 1
        for i in range(N):
            for j in range(N):
                if i == j:
                    continue
                # neighborjoining matrix definition.            sum of i row and j column
                Q_matrix[i][j] = (N - 2) * distMatrix[i][j] - sum(distMatrix[i]) - sum(distMatrix[j])
                # track the minimum value and indices in the Q matrix
                #step 2 find min val
                if Q_matrix[i][j] < minVal:
                    min_i = i
                    min_j = j
                    minVal = Q_matrix[i][j]

        # distance from min_i
        # given equation
        edge_i = 1 / 2.0 * distMatrix[min_i][min_j] + 1 / (2.0 * (N - 2)) * (sum(distMatrix[min_i]) - sum(distMatrix[min_j]))
        # distance from min_j
        edge_j = distMatrix[min_i][min_j] - edge_i

        # save the childParentDict relationship for reconstruction
        child_i = idDictionary[ids[min_i]] if ids[min_i] in idDictionary else ids[min_i]
        child_j = idDictionary[ids[min_j]] if ids[min_j] in idDictionary else ids[min_j]
        #add to dictionary
        childParentDict[child_i] = (str(nodeId), edge_i)
        childParentDict[child_j] = (str(nodeId), edge_j)
        # remove used ids from the idlist and add the new nodeid
        ids = [i for i in ids if i not in (ids[min_i], ids[min_j])] + [str(nodeId)]
        #decrement the nodeId
        nodeId -= 1

        #creates new distance matrix for next iteration
        new_distMatrix = get_new_distMatrix(distMatrix, min_i, min_j, N)

        # recursively call with the new distance matrix
        recurse(ids, new_distMatrix)

    #first call
    recurse(inputIds, inputDistMatrix)

    # root of tree from the childParentDict relations
    treeRoot = Tree(root)
    # queue for reconstruction
    queue = [treeRoot]
    # while loop to reconstruct tree
##    print("childParentDict: ",childParentDict)
##    for child, (parent, distance) in childParentDict.items():
##        print("parent: ", parent, "child: ", child)
    while len(queue) > 0:
        temp = []
        for node in queue:
            # if parent node matches current node, 
            # create a new child node and add its parent and distance
            for child, (parent, distance) in childParentDict.items():
                if parent == node.id:
                    childNode = Tree(child,parent,distance)
                    node.add_child(childNode, distance)
                    temp.append(childNode)
        # delete the node in the dict as its been used
        for node in temp:
            del childParentDict[node.id]
        # update queue for next iteration
        queue = temp

    # return the tree root
    return treeRoot

    #creates a new distMatrix to be used in the next recursion
    #decreases matrix size by 1
def get_new_distMatrix(distMatrix, min_i, min_j, N):
    updatedDist = [[0] * (N + 1) for _ in range(N + 1)]
    for i in range(N):
        for j in range(N):
            updatedDist[i][j] = distMatrix[i][j]
    # update distances to the new node
    for k in range(N):
        #apply the formula to get new distances
        updatedDist[N][k] = (0.5) * (distMatrix[min_i][k] + distMatrix[min_j][k] - distMatrix[min_i][min_j])
        updatedDist[k][N] = updatedDist[N][k]

    # Create a new distance matrix 
    new_distMatrix = [[0] * (N - 1) for _ in range(N - 1)]
    temp_i = temp_j = 0
    for i in range(N + 1):
        # continue if i or j is a min value
        if i == min_i or i == min_j:
            continue
        temp_j = 0
        for j in range(N + 1):
            # Replace these two with the new node
            if j == min_i or j == min_j:
                continue
            new_distMatrix[temp_i][temp_j] = updatedDist[i][j]
            temp_j += 1
        temp_i += 1

    return new_distMatrix

# uses preorder traversal generate the edges file
def createEdgeFile(root):
    traversal = []
    # recursive preorder traversal
    def preorder_traversal(node):
        #base case
        if root is None:
            return
        #iterate through child dictionary
        for child in node.childrenLst:
            # traverse root first
            traversal.append((node.id, child.id, child.distance))
            # recursively traverse children
            preorder_traversal(child)
    #first call
    preorder_traversal(root)
    
    print("traversal: ", traversal)
    print("end of traversal")
    #write to file
    with open('edges.txt', 'w') as f:
        for i,(parent,child,distance) in enumerate(traversal):
            traversal[i] = (int(parent),int(child),distance)
##        traversal.sort()
##        traversal.reverse()
        for i,(parent,child,distance) in enumerate(traversal):
            traversal[i] = (str(parent),str(child),distance)
        for (parent, child, distance) in traversal:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')

def CreateNewick(root, file):
    #if no children then print taxonId and distance
    if(len(root.childrenLst) == 0):
        file.write("Taxon" + str(root.id) + ':' + str(root.distance))
    #else start a new parenthesis for the next group of children
    else:
        file.write('(')
        for i in range(0, len(root.childrenLst)):
            #call again for each child
            CreateNewick(root.childrenLst[i], file)
            if(i != len(root.childrenLst) - 1):
                file.write(',')
        file.write(')')
        #0.0 is the default distance for the root
        if(root.distance != 0.0):
            file.write(':' + str(root.distance))

##MAIN##
if __name__ == '__main__':
    print("main")
    #calls main with systemInput, default file is 'example2.fna'
    #main will take 1 argument on commandline which is the input file
    filename = systemInput()
    # read fna file and return the ids and idSequences
    ids, idSequences = readFile(filename)

    # creating difference matrix
    distanceMatrix = differenceMatrix(ids, idSequences)

    # writting difference matrix to file
    CreateDistanceFile(ids, distanceMatrix)
    #making tree root using neighbor_joining
    root = neighbor_joining(ids, distanceMatrix)
    #writing files
    createEdgeFile(root)
    treeFile = open("tree.tre", "w+")
    CreateNewick(root,treeFile)
    #adds the semicolon to the end of the Newick file
    treeFile.write(";")
    print("done")

##END OF MAIN##
