import sys

sys.path.append("./PrePostProcessing/preprocess")
sys.path.append("./PrePostProcessing/postprocess")
sys.path.append("./FluxCalculators")
sys.path.append("./InterpolationSchemes")

import foamMeshReader as fmr
import geometricCalculator as geoc
from ausmPlusUp import ausmPlusUp as ausmPU
import GreenGauss as gg
import numpy as np

### ------------------------  Prepare Mesh  ------------------------ ###
path = 'C:/Users/mehme/Desktop/Studies/ExplicitSolver/TestMeshes/simpleThreeCells/constant'

points = fmr.meshFileImport(path + "/polyMesh/points")
faces = fmr.meshFileImport(path + "/polyMesh/faces")
owners = fmr.meshFileImport(path +"/polyMesh/owner")
nbours = fmr.meshFileImport(path + "/polyMesh/neighbour")
boundary = fmr.meshFileImport(path + "/polyMesh/boundary")

pointDataType = fmr.checkDataType(points)
faceDataType = fmr.checkDataType(faces)
ownerDataType = fmr.checkDataType(owners)
nboursDataType = fmr.checkDataType(nbours)
boundaryDataType = fmr.checkDataType(boundary)

parsedPoints, totalPoints = fmr.parseData(points, pointDataType)
parsedFaces, totalFaces = fmr.parseData(faces, faceDataType)
parsedOwners, totalOwners = fmr.parseData(owners, ownerDataType)
parsedNbours, totalNbours = fmr.parseData(nbours, nboursDataType)
parsedBoundary, totalBoundary = fmr.parseData(boundary, boundaryDataType)

parsedMinBoundaryLimits = []
parsedMaxBoundaryLimits = []

for boundary in parsedBoundary:
    minFaceID, maxFaceID = fmr.determineLimitsOfBoundary(boundary)
    parsedMaxBoundaryLimits.append(maxFaceID)
    parsedMinBoundaryLimits.append(minFaceID)
    
facesWCoords = fmr.findEachFaceCoordinates(parsedFaces, parsedPoints)

totalCells = fmr.findTotalCellNumbers(owners)

### ------------------------ Initializing fields ------------------------ ###

Q_init = [1,1,0,0,1] # ρ, ρu, ρv, ρw, ρE
Q_field_current = [] # Conserved quantities in every cell center for now
Q_field_later = [] # Conserved quantities in every cell center for later
fluxSolver = "AUSM+Up"

for cell in range(totalCells):
    Q_field_current.append(Q_init)
    Q_field_later.append(Q_init)
    maxEdge, minEdge = geoc.calculateExtremumEdges(cell, parsedOwners, parsedNbours, parsedFaces, facesWCoords)

### ------------------------ Starting time loop ------------------------ ###
t_end = 1
x_min = minEdge
deltaT = 0.5
n_iter = int(t_end / deltaT)
currentTime = 0

if (fluxSolver == "AUSM+Up"):
    for i in range(0, n_iter):
        print("Time is: ", currentTime)
        currentTime = currentTime + deltaT
        for cell in range(totalCells):
            totalCellSum = 0
            cellSum = 0
            print("Cell number:", cell)
            faces = fmr.findFacesOfCell(cell, parsedOwners, parsedNbours, parsedFaces)
            volume = geoc.calculateElementVolume(cell, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
            for face in faces:
                faceSum = 0
                faceID = face[0]
                faceNormal = geoc.calculateFaceNormalVector(faceID, cell, parsedOwners, facesWCoords)
                faceArea = geoc.calculateFaceArea(faceID, facesWCoords)
                # Internal Faces where flux calculations happen
                if (faceID < min(parsedMinBoundaryLimits)):
                    print("Face number: ", faceID)
                    #print("Face normal: ", faceNormal)
                    leftCellID, rightCellID = fmr.findLeftRightCells(faceID, cell, parsedOwners, parsedNbours)
                    Q_cell_L = Q_field_current[leftCellID]
                    Q_cell_R = Q_field_current[rightCellID]
                    print("Left and right cell vectors: ", Q_cell_L, Q_cell_R)
                    Q_face_L, Q_face_R = gg.calculateFaceLeftRightStateVector(Q_cell_L, 
                                                                              Q_cell_R,
                                                                              leftCellID, 
                                                                              rightCellID, 
                                                                              faceID, 
                                                                              parsedOwners, 
                                                                              parsedNbours, 
                                                                              parsedFaces, 
                                                                              facesWCoords)
                    print("Left and right face vectors: ", Q_face_L, Q_face_R)
                    faceInviscidFlux = ausmPU.invokeAUSMPlusUp(Q_face_L, Q_face_R, faceNormal, 0.8)
                    print("Face flux vector: ", faceInviscidFlux)
                    faceViscousFlux = [0,0,0,0,0]
                    faceNetFlux = np.array(faceInviscidFlux, dtype = float) - np.array(faceViscousFlux, dtype = float)
                    faceSum = faceNetFlux * faceArea                
                    cellSum = faceSum + cellSum   
                    print("Face sum: ", faceSum)
            totalCellSum = cellSum * (-1/volume) * deltaT
            print("Cell sum: ", totalCellSum)
            Q_field_later[cell] = totalCellSum
            print("For now: ", Q_field_current)
            print("For later: ", Q_field_later)
        Q_field_current = []
        for i in range(len(Q_field_later)):
            Q_field_current.append(Q_field_later[i])