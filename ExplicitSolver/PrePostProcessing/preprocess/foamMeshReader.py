def requireMeshPath():
    path = input(r"Enter path: ")
    return path

def meshFileImport(path):
    f = open(path, "r")
    data = f.read()
    content = data.split('\n')
    return content

def checkDataType(content):
    for line in content:
        if ("points;" in line):
            dataType = 0
        if ("faces;" in line):
            dataType = 1
        if ("owner;" in line):
            dataType = 2
        if ("neighbour;" in line):
            dataType = 3
        if ("boundary;" in line):
            dataType = 4
    return dataType

def parseData(content, dataType):
    parsedData = []
    boundaryNames = []
    boundaryTypes = []
    boundaryNFaces = []
    boundaryStartFaces = []
    singleDataID = 0
    totalNumber = 0
    for line in content:
        singleData = []
        if ((dataType == 0 or dataType == 1)):
            if (len(line.split('(')) != 1):
                if (line.split('(')[1] != ''):
                    singleData.append(singleDataID)
                    coordinates = line.split('(')[1].split(')')[0].split(' ')
                    for i in range(len(coordinates)):
                        singleData.append(coordinates[i])
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1               
                else:
                    totalNumber = int(content[content.index(line) - 1]);
        elif ((dataType == 2 or dataType == 3)):
            if (line.isnumeric()):
                if ('(' in content[content.index(line) + 1]):
                    totalNumber = int(line)
                else:
                    singleData.append(singleDataID)
                    singleData.append(int(line))
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1
        else:
            if ('(' in line):
                if ('inGroups' not in line):
                    totalNumber = int(content[content.index(line) - 1])
            if ('type' in line):
                boundaryTypes.append(line.split(' ')[-1].split(';')[0])
            if ('nFaces' in line):
                boundaryNFaces.append(line.split(' ')[-1].split(';')[0])
            if ('startFace' in line):
                boundaryStartFaces.append(line.split(' ')[-1].split(';')[0])
            if (line.split(' ')[-1].isalpha() and line != 'FoamFile'):
                boundaryNames.append(line.split(' ')[-1])
    for i in range(len(boundaryTypes)):
        singleData.append(boundaryNames[i])
        singleData.append(boundaryTypes[i])
        singleData.append(boundaryNFaces[i])
        singleData.append(boundaryStartFaces[i])
        parsedData.append(singleData)
        singleData = []
    return parsedData, totalNumber

path = requireMeshPath()

points = meshFileImport(path + "/polyMesh/points")
faces = meshFileImport(path + "/polyMesh/faces")
owners = meshFileImport(path +"/polyMesh/owner")
nbours = meshFileImport(path + "/polyMesh/neighbour")
boundary = meshFileImport(path + "/polyMesh/boundary")

pointDataType = checkDataType(points)
faceDataType = checkDataType(faces)
ownerDataType = checkDataType(owners)
nboursDataType = checkDataType(nbours)
boundaryDataType = checkDataType(boundary)

parsedPoints, totalPoints = parseData(points, pointDataType)
parsedFaces, totalFaces = parseData(faces, faceDataType)
parsedOwners, totalOwners = parseData(owners, ownerDataType)
parsedNbours, totalNbours = parseData(nbours, nboursDataType)
parsedBoundary, totalBoundary = parseData(boundary, boundaryDataType)

def findTotalCellNumbers(content):
    for line in content:
        if ('note' in line):
            totalCellNumber = int(line.split('nCells:')[1].split(' ')[0])
    return totalCellNumber

def readVertexCoordinates(index):
    xCoordinate = parsedPoints[index][1]
    yCoordinate = parsedPoints[index][2]
    zCoordinate = parsedPoints[index][3]
    coordinate = [xCoordinate, yCoordinate, zCoordinate]
    return coordinate

def findFacesOfCell(cellIndex):
     facesOfElement = []
     for i in range(len(parsedOwners)):
         if (parsedOwners[i][1] == cellIndex):
             facesOfElement.append(parsedFaces[i])
     for i in range(len(parsedNbours)):
         if (parsedNbours[i][1] == cellIndex):
             facesOfElement.append(parsedFaces[i])
     return facesOfElement

def determineLimitsOfBoundary(boundary):
    minFaceID = int(boundary[3])
    numberOfFaces = int(boundary[2])
    maxFaceID = minFaceID + numberOfFaces - 1
    return minFaceID, maxFaceID

def findLeftCells(faceID):
    leftCellID = parsedOwners[faceID][1]
    return leftCellID

def findRightCells(faceID):
    rightCellID = parsedNbours[faceID][1]
    return rightCellID

def findEachFaceCoordinates():
    allFaces = []
    singleFace = []
    for face in parsedFaces:
        for j in range(1,5):
            singleFace.append(readVertexCoordinates(int(face[j])))
        allFaces.append(singleFace)
        singleFace = []
    return allFaces

facesWCoords = findEachFaceCoordinates()