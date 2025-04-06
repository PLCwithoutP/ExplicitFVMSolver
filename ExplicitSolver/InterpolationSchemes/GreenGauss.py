import numpy as np
import sys

sys.path.append("./PrePostProcessing/preprocess")
sys.path.append("./PrePostProcessing/postprocess")
sys.path.append("./FluxCalculators")

import foamMeshReader as fmr
import geometricCalculator as geoc

gamma = 1.4

def gradValue(faceCenterValue, faceIndex, cellIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords):
    cellVolume = geoc.calculateElementVolume(cellIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
    testElementFaces = []
    gradQ_x = [0,0,0,0,0]
    gradQ_y = [0,0,0,0,0]
    gradQ_z = [0,0,0,0,0]
    totalGrad_x = [0,0,0,0,0]
    totalGrad_y = [0,0,0,0,0]
    totalGrad_z = [0,0,0,0,0]
    testElementFaces = fmr.findFacesOfCell(cellIndex, parsedOwners, parsedNbours, parsedFaces)
    for i in range(len(testElementFaces)):
        area = geoc.calculateFaceArea(faceIndex, facesWCoords)
        faceNormalVec = geoc.calculateFaceNormalVector(faceIndex, cellIndex, parsedOwners, facesWCoords)
        n_x = faceNormalVec[0]
        n_y = faceNormalVec[1]
        n_z = faceNormalVec[2]
        for i in range(len(faceCenterValue)):
            gradQ_x[i] = (1/cellVolume)*faceCenterValue[i]*n_x*area
            gradQ_y[i] = (1/cellVolume)*faceCenterValue[i]*n_y*area
            gradQ_z[i] = (1/cellVolume)*faceCenterValue[i]*n_z*area
        totalGrad_x[i] = totalGrad_x[i] + gradQ_x[i]
        totalGrad_y[i] = totalGrad_y[i] + gradQ_y[i]
        totalGrad_z[i] = totalGrad_z[i] + gradQ_z[i]
    totalGrad = [totalGrad_x, totalGrad_y, totalGrad_z]
    totalGrad = np.array(totalGrad).T.tolist()
    return totalGrad

    
def interpolateValue(leftCellValue, rightCellValue):
    interpolatedValue = [0,0,0,0,0]
    for i in range(len(leftCellValue)):
        interpolatedValue[i] = (leftCellValue[i] + rightCellValue[i])/2
    return interpolatedValue

def calculateFaceLeftRightStateVector(leftCellCenterValue, rightCellCenterValue, leftCellIndex, 
                                  rightCellIndex, faceIndex, parsedOwners, parsedNbours, 
                                  parsedFaces, facesWCoords):
    ''' It reconstructs state vector on the left and right of the face '''
    Q_face_l = [0,0,0,0,0]
    Q_face_r = [0,0,0,0,0]
    Q_cell_l = leftCellCenterValue
    Q_cell_r = rightCellCenterValue
    Q_face = interpolateValue(Q_cell_l, Q_cell_r)
    gradQ_left = gradValue(Q_face, faceIndex, leftCellIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
    gradQ_right = gradValue(Q_face, faceIndex, rightCellIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
    r_left = geoc.calculateInterpolationVector(leftCellIndex, faceIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
    r_right = geoc.calculateInterpolationVector(rightCellIndex, faceIndex, parsedOwners, parsedNbours, parsedFaces, facesWCoords)
    prod_left = np.dot(gradQ_left, r_left)
    prod_right = np.dot(gradQ_right, r_right)
    for i in range(len(Q_face_l)):
        Q_face_l[i] = Q_cell_l[i] + prod_left[i]
        Q_face_r[i] = Q_cell_r[i] + prod_right[i]
    return Q_face_l, Q_face_r