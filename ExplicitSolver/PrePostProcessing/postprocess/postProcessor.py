import math
import numpy as np
import os

gamma = 1.4

def calculatePressure(Q):
    e_int = (Q[4] - 0.5*Q[0]*(Q[1]**2 + Q[2]**2 + Q[3]**2))/(Q[0])
    p = e_int * Q[0] * (gamma - 1)
    return p

def calculateEnergyFromPressure(Q,p):
    rho = Q[0]
    vel = calculateVelocity(Q)
    magVel = np.linalg.norm(vel)*np.linalg.norm(vel)
    E = p/(gamma - 1) + (0.5*rho*magVel)    
    return E

def calculateVelocity(Q):
    u = Q[1]/Q[0]
    v = Q[2]/Q[0]
    w = Q[3]/Q[0]
    vel = [u, v, w]
    return vel

def calculateMach(Q):
    p = calculatePressure(Q)
    rho = Q[0]
    vel = np.linalg.norm(calculateVelocity(Q))
    a = math.sqrt((p*gamma)/rho)
    Ma = vel/a
    return Ma

def calculateSpeedOfSound(Q):
    p = calculatePressure(Q)
    rho = Q[0]
    a = math.sqrt((p*gamma)/rho)
    return a



def createScalarFieldData(time, fieldName, dimensions, fieldDataList, totalCells, parsedBoundary):
    os.system(f"mkdir {time}")
    fieldFile = open(f"{time}/{fieldName}", "w")
    headerString = r""" /*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2406                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/"""
    foamHeaderString = r"""
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;"""
    locationString = f"\n    location    \"{time}\";\n"
    objectString = f"    object      {fieldName};\n"    
    dimensionsString = f"dimensions      {dimensions};\n\n"
    internalFieldHeader = "internalField      nonuniform List<scalar>\n"
    totalCellsString = f"{totalCells}\n"
    boundaryFieldHeader = "boundaryField\n" + "{\n"
    
    fieldFile.write(headerString)
    fieldFile.write(foamHeaderString)
    fieldFile.write(locationString)
    fieldFile.write(objectString)
    fieldFile.write("""}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
""")
    fieldFile.write(dimensionsString)
    fieldFile.write(internalFieldHeader)
    fieldFile.write(totalCellsString)
    fieldFile.write("(\n")
    for i in range(len(fieldDataList)):
        fieldFile.write(f"{fieldDataList[i]}" + "\n")
    fieldFile.write(")\n" + ";\n\n")
    fieldFile.write(boundaryFieldHeader)
    for i in range(len(parsedBoundary)):
        fieldFile.write(f"      {parsedBoundary[i][0]}\n")
        fieldFile.write("      {\n")
        fieldFile.write("            value      nonuniform List<scalar>\n")
        fieldFile.write(f"{parsedBoundary[i][2]}\n")
        fieldFile.write("(\n")
        for j in range(int(parsedBoundary[i][2])):
            fieldFile.write("1\n")
        fieldFile.write(")\n" + ";\n" + "      }\n")
    fieldFile.write("}\n")
    fieldFile.close()
    return

def createPressureFieldData(time, fieldDataList):
    createScalarFieldData(time, "p", "[1 -1 2 0 0 0 0]", fieldDataList)
    return

def createMachData(time, fieldDataList):
    createScalarFieldData(time, "Ma", "[0 0 0 0 0 0 0]", fieldDataList)
    return

def createUData(time, fieldDataList):
    createScalarFieldData(time, "u", "[0 1 -1 0 0 0 0]", fieldDataList)
    return

def createVData(time, fieldDataList):
    createScalarFieldData(time, "v", "[0 1 -1 0 0 0 0]", fieldDataList)
    return

def createAData(time, fieldDataList):
    createScalarFieldData(time, "a", "[0 1 -1 0 0 0 0]", fieldDataList)
    return