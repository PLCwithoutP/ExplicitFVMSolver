import foamMeshReader as fmr
import numpy as np

def calculateFaceArea(faceIndex):
    face = fmr.facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    vector1 = [0, 0, 0]
    vector2 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
        vector1[i] = point1[i] - point2[i]
        vector2[i] = point1[i] - point4[i]
    area = (np.cross(vector2, vector1))
    area = np.linalg.norm(area)
    return area

def calculateFaceNormalVector(faceIndex, cellIndex):
    face = fmr.facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    pointCenter = calculateFaceCenterPoint(faceIndex)
    vector1 = [0, 0, 0]
    vector2 = [0, 0, 0]
    vector3 = [0, 0, 0]
    vector4 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
        vector1[i] = point1[i] - pointCenter[i]
        vector2[i] = point2[i] - pointCenter[i]
        vector3[i] = point3[i] - pointCenter[i]
        vector4[i] = point4[i] - pointCenter[i]
    faceOrthoVector = 0.5*(np.cross(vector1, vector2) + np.cross(vector2, vector3) 
                + np.cross(vector3, vector4) + np.cross(vector4, vector1))
    faceNormalVector = faceOrthoVector/np.linalg.norm(faceOrthoVector)
    if not (cellIndex == fmr.parsedOwners[faceIndex][1]):
        for i in range(len(faceNormalVector)):
            faceNormalVector[i] = (-1)*faceNormalVector[i]
    return faceNormalVector
    
def calculateEdgeVector(faceIndex):
    face = fmr.facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    vector1 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        vector1[i] = point2[i] - point1[i]
    return vector1

def calculateFaceCenterPoint(faceIndex):
    face = fmr.facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
    xMax = max(point1[0], point2[0], point3[0], point4[0])
    yMax = max(point1[1], point2[1], point3[1], point4[1])
    zMax = max(point1[2], point2[2], point3[2], point4[2])
    xMin = min(point1[0], point2[0], point3[0], point4[0])
    yMin = min(point1[1], point2[1], point3[1], point4[1])
    zMin = min(point1[2], point2[2], point3[2], point4[2])
    xCoord = (xMax + xMin)/2
    yCoord = (yMax + yMin)/2
    zCoord = (zMax + zMin)/2
    centerPoint = [xCoord, yCoord, zCoord]
    return centerPoint

def calculateElementCenterPoint(cellIndex):
    facesOfElement = []
    xCenter = 0
    yCenter = 0
    zCenter = 0
    facesOfElement = fmr.findFacesOfCell(cellIndex)
    for i in range(len(facesOfElement)):
        faceCenter = calculateFaceCenterPoint(int(facesOfElement[i][0]))
        xCenter = xCenter + faceCenter[0]
        yCenter = yCenter + faceCenter[1]
        zCenter = zCenter + faceCenter[2]
    xCenter = xCenter/6
    yCenter = yCenter/6
    zCenter = zCenter/6
    centerPoint = [xCenter, yCenter, zCenter]
    return centerPoint

def calculateElementVolumePoly(cellIndex):
    facesOfElement = []
    volume = 0
    facesOfElement = fmr.findFacesOfCell(cellIndex)
    for face in facesOfElement:
        faceArea = calculateFaceArea(face[0])
        faceCenter = calculateFaceCenterPoint(face[0])
        faceNormal = calculateFaceNormalVector(face[0], cellIndex)
        volume = volume + (faceArea*(np.dot(faceCenter,faceNormal)))
    return volume/3
    

def calculateInterpolationVector(cellIndex, faceIndex):
    cellCenter = calculateElementCenterPoint(cellIndex)
    faceCenter = calculateFaceCenterPoint(faceIndex)
    vector = [0,0,0]
    for i in range(len(cellCenter)):
        vector[i] = cellCenter[i] - faceCenter[i]
    return vector