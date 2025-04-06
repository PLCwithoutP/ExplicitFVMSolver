def meshFileImport(path):
    ''' This function imports polyMesh file where is given in path'''
    f = open(path, "r")
    data = f.read()
    content = data.split('\n')
    return content

def checkDataType(content):
    ''' This function checks whether read data is point, face, or cell'''
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
    ''' This function processes given polyMesh file into mesh data'''
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
        elif (dataType == 2):
            if (line.isnumeric()):
                if ('(' in content[content.index(line) + 1]):
                    totalNumber = int(line)
                else:
                    singleData.append(singleDataID)
                    singleData.append(int(line))
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1
        elif (dataType == 3):
            if ('(' and ')' in line):
                totalNumber = int(line.split("(")[0])
                data = line.split("(")[1].split(")")[0]
                for i in range(len(data.split(" "))):
                    singleData.append(singleDataID)
                    singleData.append(int(data.split(" ")[i]))
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1
                    singleData = []
            elif (line.isnumeric()):
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

def findTotalCellNumbers(content):
    ''' Find the total number of cells '''
    for line in content:
        if ('note' in line):
            totalCellNumber = int(line.split('nCells:')[1].split('nFaces:')[0])
    return totalCellNumber

def readVertexCoordinates(pointIndex, parsedPoints):
    ''' This function finds the coordinates of a given point'''
    xCoordinate = parsedPoints[pointIndex][1]
    yCoordinate = parsedPoints[pointIndex][2]
    zCoordinate = parsedPoints[pointIndex][3]
    coordinate = [xCoordinate, yCoordinate, zCoordinate]
    return coordinate

def findFacesOfCell(cellIndex, parsedOwners, parsedNbours, parsedFaces):
    ''' This function finds the faces of a given cell''' 
    facesOfElement = []
    for i in range(len(parsedOwners)):
        if (parsedOwners[i][1] == cellIndex):
            facesOfElement.append(parsedFaces[i])
    for i in range(len(parsedNbours)):
        if (parsedNbours[i][1] == cellIndex):
            facesOfElement.append(parsedFaces[i])
    return facesOfElement

def determineLimitsOfBoundary(boundary):
    ''' This function parses boundary limits, and finds number of faces'''
    minFaceID = int(boundary[3])
    numberOfFaces = int(boundary[2])
    maxFaceID = minFaceID + numberOfFaces - 1
    return minFaceID, maxFaceID

def findLeftRightCells(faceID, cellID, parsedOwners, parsedNbours):
    ''' This function finds left and right cells of a face'''
    leftCellID = cellID
    for faces in parsedNbours:
        checkFaceID = faces[0]
        checkCellID = faces[1]
        if (faceID == checkFaceID and cellID != checkCellID):
            rightCellID = checkCellID
    for faces in parsedOwners:
        checkFaceID = faces[0]
        checkCellID = faces[1]
        if (faceID == checkFaceID and cellID != checkCellID):
            rightCellID = checkCellID
    return leftCellID, rightCellID


def findEachFaceCoordinates(parsedFaces, parsedPoints):
    ''' This function finds face coordinates from parsed points'''
    allFaces = []
    singleFace = []
    for face in parsedFaces:
        for j in range(1,5):
            singleFace.append(readVertexCoordinates(int(face[j]), parsedPoints))
        allFaces.append(singleFace)
        singleFace = []
    return allFaces