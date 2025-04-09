import math
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

def plot_shock_tube(Q_field_current, parsedPoints, totalCells, currentTime):
    """
    Plots U, rho, p vs. x (one-dimensional spatial coordinate) for Sod's shock tube problem.
    
    Parameters:
    - Q_field_current: List of conserved quantities at each cell center
    - parsedPoints: List of points representing the mesh coordinates
    - totalCells: Total number of cells in the mesh
    - gamma: Specific heat ratio (default is 1.4 for air)
    """
    # Initialize arrays for U, rho, and p
    rho = np.zeros(totalCells)
    U = np.zeros(totalCells)
    p = np.zeros(totalCells)
    x = np.zeros(totalCells)
    
    # Extract relevant quantities from the conserved quantities
    for cell in range(totalCells):
        # Unpack conserved quantities for the cell
        rho[cell] = Q_field_current[cell][0]  # rho
        u = Q_field_current[cell][1] / rho[cell]  # u = rho_u / rho
        U[cell] = u  # Velocity in the x-direction
        rhoE = Q_field_current[cell][4]  # rhoE
        # Calculate pressure using the ideal gas law (p = (gamma - 1) * (rhoE - 0.5*u^2))
        p[cell] = (gamma - 1) * (rhoE - 0.5 * u**2)
        # Get the x-coordinate of the cell center
        x[cell] = parsedPoints[cell][0]  # x-coordinate for 1D problem

    # Create the plots
    fig, axs = plt.subplots(3, 1, figsize=(10, 12))
    
    # Plot U vs x
    axs[0].plot(x, U, label="Velocity (U)", color='b')
    axs[0].set_title(f'Velocity (U) vs x, time = {currentTime} s')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('Velocity (U)')
    axs[0].grid(True)
    
    # Plot rho vs x
    axs[1].plot(x, rho, label="Density (rho)", color='g')
    axs[1].set_title('Density (rho) vs x')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('Density (rho)')
    axs[1].grid(True)
    
    # Plot p vs x
    axs[2].plot(x, p, label="Pressure (p)", color='r')
    axs[2].set_title('Pressure (p) vs x')
    axs[2].set_xlabel('x')
    axs[2].set_ylabel('Pressure (p)')
    axs[2].grid(True)
    
    # Adjust layout to avoid overlap
    plt.tight_layout()
    plt.show()

def createLogFile():
    os.system("cd ../..")
    f = open('logFile', 'w')
    f.write('\n')
    f.write('Solution has been started!\n')
    f.write('\n')
    f.close()
    os.system("cd PrePostProcessing/postprocess")

def logData(data):
    os.system("cd ../..")
    f = open('logFile', 'a+')
    f.write(data)
    f.write('\n')
    f.close()
    os.system("cd PrePostProcessing/postprocess")