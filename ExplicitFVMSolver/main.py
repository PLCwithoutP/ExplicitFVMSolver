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
import postProcessor as pp
### ------------------------  Prepare Mesh  ------------------------ ###
path = './TestMeshes/shockTube/constant'

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

F_inviscid = [] # inviscid flux data for each face
F_inviscidOther = [] # inviscid flux data for each face
F_viscous = [] # viscous flux data for each face

for i in range(totalFaces): # For each face
    F_inviscid.append(np.zeros(5))
    F_inviscidOther.append(np.zeros(5))
    F_viscous.append(np.zeros(5))

Q_init_left = [1,0,0,0,25] # ρ, ρu, ρv, ρw, ρE
#Q_init_right = Q_init_left
Q_init_right = [0.125,0,0,0,2.5]
Q_field_current = [] # Conserved quantities in every cell center for now
Q_field_later = [] # Conserved quantities in every cell center for later
fluxSolver = "AUSM+Up"

for cell in range(0,int(totalCells/2)):
    Q_field_current.append(Q_init_left.copy())
    Q_field_later.append(np.zeros(5))
    maxEdge, minEdge = geoc.calculateExtremumEdges(cell, parsedOwners, parsedNbours, parsedFaces, facesWCoords)

for cell in range(int(totalCells/2), int(totalCells)):
    Q_field_current.append(Q_init_right.copy())
    Q_field_later.append(np.zeros(5))
    maxEdge, minEdge = geoc.calculateExtremumEdges(cell, parsedOwners, parsedNbours, parsedFaces, facesWCoords)

### ------------------------ Starting time loop ------------------------ ###
t_end = 1
x_min = minEdge
deltaT = 0.001
n_iter = int(t_end / deltaT)
currentTime = 0

pp.createLogFile()

if (fluxSolver == "AUSM+Up"):
    for i in range(0, n_iter):
        print("Time is: ", currentTime)
        pp.logData(f"Time = {currentTime}")
        print("********************************************************************************")
        pp.logData("********************************************************************************")
        print("--------------------------------------------------------------------------------")
        pp.logData("--------------------------------------------------------------------------------")
        print("Calculation phase is started!")
        pp.logData("Calculation phase is started!")
        print("--------------------------------------------------------------------------------")
        pp.logData("--------------------------------------------------------------------------------")
        pp.plot_shock_tube(Q_field_current, parsedPoints, totalCells, currentTime)
        currentTime = currentTime + deltaT
        for face in parsedFaces:
            faceID = face[0]
            # Orientation is not taken into the account
            agLeftCellID, agRightCellID = fmr.findCellsOfFace(faceID, parsedOwners, parsedNbours)
            faceNormal = geoc.calculateFaceNormalVector(faceID, agLeftCellID, parsedOwners, facesWCoords)
            if (agRightCellID != -1):
                inverseFaceNormal = geoc.calculateFaceNormalVector(faceID, agRightCellID, parsedOwners, facesWCoords)
            faceArea = geoc.calculateFaceArea(faceID, facesWCoords)
            # Internal Faces where flux calculations happen
            if (faceID < min(parsedMinBoundaryLimits)):
                #print("Face number: ", faceID)
                pp.logData(f"Internal face number: {faceID}")
                #print("Face normal: ", faceNormal)
                pp.logData(f"Face normal: {faceNormal}")
                #print("Inverse face normal: ", inverseFaceNormal)
                pp.logData(f"Inverse face normal: {inverseFaceNormal}")
                Q_cell_L = Q_field_current[agLeftCellID]
                Q_cell_R = Q_field_current[agRightCellID]
                #print("Left and right cell vectors: ", Q_cell_L, Q_cell_R)
                pp.logData(f"Left and right cell vectors: {Q_cell_L} {Q_cell_R}")
                Q_face_L, Q_face_R = Q_cell_L, Q_cell_R
# =============================================================================
#                 Q_face_L, Q_face_R = gg.calculateFaceLeftRightStateVector(Q_cell_L, 
#                                                                           Q_cell_R,
#                                                                           agLeftCellID, 
#                                                                           agRightCellID,
#                                                                           Q_field_current,
#                                                                           faceID, 
#                                                                           parsedOwners, 
#                                                                           parsedNbours, 
#                                                                           parsedFaces, 
#                                                                           facesWCoords)
# =============================================================================
                #print("Left and right face vectors: ", Q_face_L, Q_face_R)
                pp.logData(f"Left and right face vectors: {Q_face_L} {Q_face_R}")
                faceInviscidFlux = ausmPU.invokeAUSMPlusUp(Q_face_L, Q_face_R, faceNormal, 0.01)
                pp.logData(f"Left face normal flux: {faceInviscidFlux}")
                #print("For face",faceID,"Face inv. flux value:", faceInviscidFlux)
                faceInviscidFluxInverse = ausmPU.invokeAUSMPlusUp(Q_face_L, Q_face_R, inverseFaceNormal, 0.01)
                pp.logData(f"Right face normal flux: {faceInviscidFluxInverse}")
                #print("Inverse direction for face",faceID,"Face inv. flux value:", faceInviscidFluxInverse)
                for i in range(5):
                    F_inviscid[faceID][i] = faceInviscidFlux[i].copy()
                    F_inviscidOther[faceID][i] = faceInviscidFluxInverse[i].copy()
                #print("Total face inviscid flux: ", F_inviscid[faceID])
            elif (faceID >= min(fmr.determineLimitsOfBoundary(parsedBoundary[0])) 
                  and faceID < max(fmr.determineLimitsOfBoundary(parsedBoundary[0]))):
                #print("This is sides boundary!")
                #pp.logData(f"Boundary face number: {faceID}")
                Q_ghost_cell = Q_field_current[agLeftCellID].copy()
                Q_cell_L = Q_field_current[agLeftCellID].copy()
                Q_face_L, Q_face_R = Q_cell_L, Q_ghost_cell
                #print("Left and right face vectors: ", Q_face_L, Q_face_R)
                faceInviscidFlux = ausmPU.invokeAUSMPlusUp(Q_face_L, Q_face_R, faceNormal, 0.01)
                #print("For face",faceID,"Face inv. flux value:", faceInviscidFlux)
                for i in range(5):
                    F_inviscid[faceID][i] = faceInviscidFlux[i].copy()
            elif (faceID > min(fmr.determineLimitsOfBoundary(parsedBoundary[0])) 
                  and faceID <= max(fmr.determineLimitsOfBoundary(parsedBoundary[0]))):
                #print("This is sides boundary!")
                #pp.logData(f"Boundary face number: {faceID}")
                Q_ghost_cell = Q_field_current[agLeftCellID].copy()
                Q_cell_L = Q_field_current[agLeftCellID].copy()
                Q_face_L, Q_face_R = Q_cell_L, Q_ghost_cell
                #print("Left and right face vectors: ", Q_face_L, Q_face_R)
                faceInviscidFlux = ausmPU.invokeAUSMPlusUp(Q_face_L, Q_face_R, faceNormal, 0.01)
                #print("For face",faceID,"Face inv. flux value:", faceInviscidFlux)
                for i in range(5):
                    F_inviscid[faceID][i] = faceInviscidFlux[i].copy()
            else:
                F_inviscid[faceID] = np.array([0,0,0,0,0])
                F_inviscidOther[faceID] = np.array([0,0,0,0,0])
        print("--------------------------------------------------------------------------------")
        pp.logData("--------------------------------------------------------------------------------")
        print("Evaluation phase is started!")
        pp.logData("Evaluation phase is started!")
        print("--------------------------------------------------------------------------------")
        pp.logData("--------------------------------------------------------------------------------")
        for cell in range(totalCells):
            #print("Cell number: ", cell)
            pp.logData(f"Cell number: {cell}")
            faces = fmr.findFacesOfCell(cell, parsedOwners, parsedNbours, parsedFaces)
            volume = geoc.calculateElementVolume(cell, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
            totalCellSum = [0,0,0,0,0]
            cellSum = [0,0,0,0,0]
            for face in faces:
                faceSum = [0,0,0,0,0]
                flux = [0,0,0,0,0]
                faceIDEval = face[0]
                faceArea = geoc.calculateFaceArea(faceIDEval, facesWCoords)
                if (faceIDEval < min(parsedMinBoundaryLimits)):
                    #print("Face number: ", faceIDEval)
                    pp.logData(f"Face number: {faceIDEval}")
                    leftCellIDEval, rightCellIDEval = fmr.findLeftRightCells(faceIDEval, cell, parsedOwners, parsedNbours)
                    agLeftCellIDEval, agRightCellIDEval = fmr.findCellsOfFace(faceIDEval, parsedOwners, parsedNbours)   
                    if (leftCellIDEval == agLeftCellIDEval):
                        #print("Exiting to cell!")
                        flux = F_inviscid[faceIDEval].copy()                    
                        #print(flux)
                        pp.logData(f"Flux is : {flux}")
                    elif (leftCellIDEval == agRightCellIDEval):
                        #print("Entering to cell!")
                        flux = F_inviscidOther[faceIDEval].copy()
                        #print(flux)
                        pp.logData(f"Flux is : {flux}")
                    for i in range(len(flux)):
                        faceSum[i] = flux[i] * faceArea                
                        cellSum[i] = faceSum[i] + cellSum[i]   
                    #print("Face sum: ", faceSum)
                else:
                    flux = F_inviscid[faceIDEval].copy()
                    F_inviscid[faceIDEval] = np.array([0,0,0,0,0])
                    for i in range(len(flux)):
                        faceSum[i] = flux[i] * faceArea                
                        cellSum[i] = faceSum[i] + cellSum[i]                       
            for i in range(len(Q_field_later[cell])):
                totalCellSum[i] = cellSum[i] * (-1/volume) * deltaT
                Q_field_later[cell][i] = Q_field_current[cell][i] + totalCellSum[i]
            pp.logData(f"Cell center value: {totalCellSum}")
        for i in range(len(Q_field_later)):
            Q_field_current[i] = Q_field_later[i].copy()
        #print("---------------------------------------------------------------------------")
        #print(F_inviscid)
        #print("---------------------------------------------------------------------------")
        #print(F_inviscidOther)
        #print("---------------------------------------------------------------------------")
        for face in range(totalFaces):
            F_inviscid[face] = np.array([0,0,0,0,0])
            F_inviscidOther[face] = np.array([0,0,0,0,0])
        #print("***************************************************************************")
        #print(F_inviscid)
        #print("***************************************************************************")
        #print(F_inviscidOther)
        #print("***************************************************************************")
        #print(Q_field_current)
print("...................................................................................")
pp.logData("...................................................................................")
print("Calculation is completed.")
pp.logData("Calculation is completed.")


