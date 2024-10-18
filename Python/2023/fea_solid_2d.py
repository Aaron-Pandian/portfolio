# Aaron Pandian (ANP3238) COE321K
# The goal of this project is to determine the nodal displacements, element strains, and element stresses for an arbitrary two-dimensional solid structure using 3-noded triangular elements.

"""
The structure of the input files are as follows:
"Nodes.txt": node number, x coordinate, y coordinate, z coordinate; point of origin is the bottom left most node
"Nodal_Displacements.txt": node number, degree of displacement, and boundary condition or known diaplcement at node
"Elements.txt": element number, node 1, node 2, E, A
"Applied_Nodal_Force.txt": node number, degree of force, coefficient of P
"DOF.txt": node number, local degree of freedom, global degree of freedom number
"""

from scipy.spatial import Delaunay
from array import array
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import math
from math import pi
dimentions = 2

# Array Creation
nodes = []
# Node Number, X, Y
elements = []
# Element Number, Node 1, Node 2, Node 3, E, Poisons Ratio
displacements = []
# Node number, Local DOF Number, Displacement Constraint
extforces = []
# Node Number, Local DOF Number, P Coefficient
dof = []
# Node Number, DOF #, gcon Global DOF (Base)
gcon = []
# (Node and Local DOF #), Global DOF #
gconold = []

with open("Nodes.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        nodes.append(res)

with open("Elements.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        elements.append(res)

with open("Nodal_Displacements.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        displacements.append(res)

with open("Applied_Nodal_Forces.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        extforces.append(res)

with open("DOF.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        dof.append(res)

# gcon array creation
for local_dof in dof:
    gcon.append([[local_dof[0], local_dof[1]], local_dof[2]])
    gconold.append([[local_dof[0], local_dof[1]], local_dof[2]])

# gcon array specification
for displacement in displacements:
    ndofs = len(dof)
    dofnum = 0
    dof_index = [displacement[0], displacement[1]]

    # Find correct node and local dof pair in gcon to find dofnum
    for element in gcon:
        if element[0] == dof_index:
            dofnum = element[1]

    # Checking gcon value against dofnum
    for element in gcon:
        if element[1] > dofnum:
            element[1] = element[1] - 1

    for element in gcon:
        if element[0] == dof_index:
            element[1] = len(nodes) * dimentions

# Force Global Matrix Assembly
x = dimentions*len(nodes)
kglobal = np.zeros((x, x), dtype=float)
fglobal = np.zeros((x), dtype=float)

for force in extforces:
    for element in gconold:
        if element[0] == [force[0], force[1]]:
            # Account for the fact that indexing starts at 0, with kglobal numbers in mind
            fglobal[element[1]-1] += force[2]

# Initializing all displacements for node number and dof number of displacement
uarray = []
for element in displacements:
    uarray.append([[element[0], element[1]], element[2]])

# Consistant for each element
E = float(elements[0][4])
PoisonsRatio = float(elements[0][5])

#  Element Matrices Calculation
for bar in elements:
    # Element find E, A, L, dx, dy, c, s
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    globalnodenum3 = bar[3]

    x1 = 0
    x2 = 0
    x3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            x3 = float(node[1])
        if node[0] == globalnodenum2:
            x2 = float(node[1])
        if node[0] == globalnodenum1:
            x1 = float(node[1])
    # For Y distances between nodes
    y1 = 0
    y2 = 0
    y3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            y3 = float(node[2])
        if node[0] == globalnodenum2:
            y2 = float(node[2])
        if node[0] == globalnodenum1:
            y1 = float(node[2])

    Atri = 0.5 * ((x2*y3)-(x3*y2)+(x3*y1)-(x1*y3)+(x1*y2)-(x2*y1))

    """ Unused Matrix
    AtriMatrix = 0.5 * np.array([[1, x1, y1],
                                [1, x2, y2],
                                 [1, x3, y3]])
    """

    # Terms for N1
    a1 = (1/(2*Atri))*((x2*y3)-(x3*y2))
    b1 = (1/(2*Atri))*((y2)-(y3))
    c1 = (1/(2*Atri))*((x3)-(x2))

    # Terms for N2
    a2 = (1/(2*Atri))*((x3*y1)-(x1*y3))
    b2 = (1/(2*Atri))*((y3)-(y1))
    c2 = (1/(2*Atri))*((x1)-(x3))

    # Terms for N3
    a3 = (1/(2*Atri))*((x1*y2)-(x2*y1))
    b3 = (1/(2*Atri))*((y1)-(y2))
    c3 = (1/(2*Atri))*((x2)-(x1))

    # B Matrix
    B = np.array([[b1, 0, b2, 0, b3, 0],
                  [0, c1, 0, c2, 0, c3],
                  [c1, b1, c2, b2, c3, b3]])

    EnoteStress = E
    VnoteStress = PoisonsRatio
    EnoteStrain = E/(1-(PoisonsRatio**2))
    VnoteStrain = PoisonsRatio/(1-PoisonsRatio)

    # C Matrix Assembly
    C = np.array([[((EnoteStress)/(1-(VnoteStress**2))), ((VnoteStress*EnoteStress)/(1-(VnoteStress**2))), 0],
                  [((VnoteStress*EnoteStress)/(1-(VnoteStress**2))),
                   ((EnoteStress)/(1-(VnoteStress**2))), 0],
                  [0, 0, ((EnoteStress)/(2*(1+VnoteStress)))]])

    Bt = np.transpose(B)

    # Element Matrix Assembly
    part1 = Atri * Bt
    part2 = np.dot(part1, C)
    Kele = np.dot(part2, B)

    # Forming Kglobal Setup
    lnode1 = [globalnodenum1, globalnodenum2, globalnodenum3]

    ndim = [1, 2]
    gdof1arr = []
    for node in lnode1:
        for dim in ndim:
            # 8 points of construction
            ldof1 = dimentions*(node-1) + dim
            gdof1arr.append(ldof1)

    # Creating Global Matrix using each Kele for each element
    index1 = 0
    for gnode in gdof1arr:
        index2 = 0
        for gnode2 in gdof1arr:
            kglobal[gnode-1][gnode2-1] += Kele[index1][index2]
            index2 += 1
        index1 += 1

# Solving
uglobal = np.zeros((x), dtype=float)

# Finding global DOF where boundary conditions are applied
gdofrednum = []
gdofrednumpair = []
for bc in displacements:
    gdofrednumpair.append(bc[2])
    for element in gconold:
        if (element[0] == [bc[0], bc[1]]):
            gdofrednum.append(element[1]-1)

kred = kglobal.tolist()
fred = fglobal.tolist()
ured = uglobal.tolist()

# reduces the matricies and arrays by row
remover = 0
for num in gdofrednum:
    ured.pop(num-remover)
    fred.pop(num-remover)
    kred.pop(num-remover)
    remover += 1

# reduces the matrix by column and subtracts the removed element times the displacement from the force in the same column
fredindex = 0
for element in kred:
    remover = 0
    for num in gdofrednum:
        fred[fredindex] -= (element[num-remover] * gdofrednumpair[remover])
        element.pop(num-remover)
        remover += 1
    fredindex += 1

fred = np.transpose(fred)
fred = np.array(fred)
kred = np.array(kred)
kinv = np.linalg.inv(kred)
uans = np.dot(kinv, fred)

# Creating uglobal matrix, where boundary conditions and uans are combined to respective nodal displacements
index = 0
index2 = 0
for element in uglobal:
    if index in gdofrednum:
        for element in gconold:
            if element[1] == index + 1:
                pair = element[0]
        for dis in uarray:
            if pair == dis[0]:
                mvmt = dis[1]
        # uglobal[index] = round(mvmt, 5)
        uglobal[index] = mvmt
    else:
        # uglobal[index] = round(uans[index2], 5)
        uglobal[index] = uans[index2]
        index2 += 1
    index += 1

# Better name
ufinal = uglobal

# Post processing
# Stain stress, and internal force calculation
n = len(elements)
stress = np.zeros((n, 3), dtype=float)
strain = np.zeros((n, 3), dtype=float)

sindex = 0
for bar in elements:
    # Element find E, A, L, dx, dy, c, s
    E = float(bar[4])
    PoisonsRatio = float(bar[5])
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    globalnodenum3 = bar[3]

    x1 = 0
    x2 = 0
    x3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            x3 = node[1]
        if node[0] == globalnodenum2:
            x2 = node[1]
        if node[0] == globalnodenum1:
            x1 = node[1]
    # For Y distances between nodes
    y1 = 0
    y2 = 0
    y3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            y3 = node[2]
        if node[0] == globalnodenum2:
            y2 = node[2]
        if node[0] == globalnodenum1:
            y1 = node[2]

    Atri = 0.5 * ((x2*y3)-(x3*y2)+(x3*y1)-(x1*y3)+(x1*y2)-(x2*y1))

    # Terms for N1
    a1 = (1/(2*Atri))*((x2*y3)-(x3*y2))
    b1 = (1/(2*Atri))*((y2)-(y3))
    c1 = (1/(2*Atri))*((x3)-(x2))

    # Terms for N2
    a2 = (1/(2*Atri))*((x3*y1)-(x1*y3))
    b2 = (1/(2*Atri))*((y3)-(y1))
    c2 = (1/(2*Atri))*((x1)-(x3))

    # Terms for N3
    a3 = (1/(2*Atri))*((x1*y2)-(x2*y1))
    b3 = (1/(2*Atri))*((y1)-(y2))
    c3 = (1/(2*Atri))*((x2)-(x1))

    # B Matrix
    B = np.array([[b1, 0, b2, 0, b3, 0],
                  [0, c1, 0, c2, 0, c3],
                  [c1, b1, c2, b2, c3, b3]])

    EnoteStress = E
    VnoteStress = PoisonsRatio
    EnoteStrain = E/(1-(PoisonsRatio**2))
    VnoteStrain = PoisonsRatio/(1-PoisonsRatio)

    # C Matrix Assembly
    C = np.array([[((EnoteStress)/(1-(VnoteStress**2))), ((VnoteStress*EnoteStress)/(1-(VnoteStress**2))), 0],
                  [((VnoteStress*EnoteStress)/(1-(VnoteStress**2))),
                   ((EnoteStress)/(1-(VnoteStress**2))), 0],
                  [0, 0, ((EnoteStress)/(2*(1+VnoteStress)))]])

    Uele = np.zeros((6, 1), dtype=float)

    for check in gconold:
        if [globalnodenum1, 1] == check[0]:
            Uele[0][0] = ufinal[check[1]-1]
        elif [globalnodenum1, 2] == check[0]:
            Uele[1][0] = ufinal[check[1]-1]
        elif [globalnodenum2, 1] == check[0]:
            Uele[2][0] = ufinal[check[1]-1]
        elif [globalnodenum2, 2] == check[0]:
            Uele[3][0] = ufinal[check[1]-1]
        elif [globalnodenum3, 1] == check[0]:
            Uele[4][0] = ufinal[check[1]-1]
        elif [globalnodenum3, 2] == check[0]:
            Uele[5][0] = ufinal[check[1]-1]

    eleStrain = np.dot(B, Uele)
    eleStress = np.dot(C, eleStrain)

    index = 0
    for items in lnode1:
        stress[sindex][index] = eleStress[index]
        strain[sindex][index] = eleStrain[index]
        index += 1
    sindex += 1


print("The following is the nodal displacement array: ")
print(ufinal)
print("\nThe following is the stress matrix: ")
print(stress)

# ______________________________________________________________________________________________________________________

# Final Continuation
# Plotting undeformed shape --------------------------------------------------------------------------------------------
xvalues = [nodes[0][1], nodes[1][1],
           nodes[2][1], nodes[3][1], nodes[4][1]]
yvalues = [nodes[0][2], nodes[1][2],
           nodes[2][2], nodes[3][2], nodes[4][2]]
ymax = max(yvalues)
ymin = min(yvalues)
xmax = max(xvalues)
xmin = min(xvalues)

xlengthold = xmin - xmax
yheightold = ymax - ymin

xelement1 = np.array([nodes[0][1], nodes[1][1], nodes[2][1], nodes[0][1]])
yelement1 = np.array([nodes[0][2], nodes[1][2], nodes[2][2], nodes[0][2]])
plt.plot(xelement1, yelement1)

xelement2 = np.array([nodes[0][1], nodes[2][1], nodes[3][1], nodes[0][1]])
yelement2 = np.array([nodes[0][2], nodes[2][2], nodes[3][2], nodes[0][2]])
plt.plot(xelement2, yelement2)

xelement3 = np.array([nodes[2][1], nodes[1][1], nodes[4][1], nodes[2][1]])
yelement3 = np.array([nodes[2][2], nodes[1][2], nodes[4][2], nodes[2][2]])
plt.plot(xelement3, yelement3)

xelement4 = np.array([nodes[2][1], nodes[4][1], nodes[3][1], nodes[2][1]])
yelement4 = np.array([nodes[2][2], nodes[4][2], nodes[3][2], nodes[2][2]])
plt.plot(xelement4, yelement4)

plt.xlabel("X")
plt.ylabel("Y")
plt.suptitle("Undeformed Shape")
plt.show()

# Plotting deformed shape -----------------------------------------------------------------------------------------------------------
scalefactor = 50
scaledu = ufinal/scalefactor

u1, v1, u2, v2, u3, v3, u4, v4, u5, v5 = scaledu

xelement1 = np.array([nodes[0][1] + u1, nodes[1][1] + u2,
                     nodes[2][1] + u3, nodes[0][1] + u1])
yelement1 = np.array([nodes[0][2] + v1, nodes[1][2] + v2,
                     nodes[2][2] + v3, nodes[0][2] + v1])
plt.plot(xelement1, yelement1)

xelement2 = np.array([nodes[0][1] + u1, nodes[2][1] + u3,
                     nodes[3][1] + u4, nodes[0][1] + u1])
yelement2 = np.array([nodes[0][2] + v1, nodes[2][2] + v3,
                     nodes[3][2] + v4, nodes[0][2] + v1])
plt.plot(xelement2, yelement2)

xelement3 = np.array([nodes[2][1] + u3, nodes[1][1] + u2,
                     nodes[4][1] + u5, nodes[2][1] + u3])
yelement3 = np.array([nodes[2][2] + v3, nodes[1][2] + v2,
                     nodes[4][2] + v5, nodes[2][2] + v3])
plt.plot(xelement3, yelement3)

xelement4 = np.array([nodes[2][1] + u3, nodes[4][1] + u5,
                     nodes[3][1] + u4, nodes[2][1] + u3])
yelement4 = np.array([nodes[2][2] + v3, nodes[4][2] + v5,
                     nodes[3][2] + v4, nodes[2][2] + v3])
plt.plot(xelement4, yelement4)

plt.xlabel("X")
plt.ylabel("Y")
plt.suptitle("Deformed Shape")
plt.show()

# Code for Arbitraty 2D Solid Structure Final __________________________________________________________________________________________
# Mesh Plotter taken from Cyprien Rusu, and Node Values taken from MATLAB plot
Nodes = []
mat = scipy.io.loadmat('data')
sorted(mat.keys())
X = mat['X']
Y = mat['Y']

elewidth = 8  # From MATLAB File
for index in range(0, elewidth+1):
    for index2 in range(0, elewidth+1):
        Nodes.append([X[index][index2], Y[index][index2]])

points = np.array(Nodes)

# Create Elements
tri = Delaunay(points, False, False, None)
"""
plt.triplot(points[:, 0], points[:, 1], tri.simplices)
plt.plot(points[:, 0], points[:, 1], 'o')
plt.show()
"""

# Cleanup Mesh
# Create a set of points on a circle of diameter 0.0195
p = []
r2 = 1-.01
for x in np.linspace(0, r2, 100):
    p.append([x, math.sqrt(r2**2-x**2)])

# Find the elements which contain those points
tri.find_simplex(p)
x = tri.find_simplex(p)

# Feed the result of the previous function to the np.delete method
# Create a new set of elements without the problematic elements
mesh = np.delete(tri.simplices, x, 0)

plt.triplot(points[:, 0], points[:, 1], mesh)
plt.plot(points[:, 0], points[:, 1], 'o')
plt.xlabel("X/R")
plt.ylabel("Y/R")
plt.suptitle("Undeformed Mesh Shape")
plt.show()

# Creating new file to analyze - Only need first run
"""
nb_nodes = len(points)
nb_elements = len(mesh)

file = open("meshnodes.txt", "w")
for i, node in enumerate(Nodes):
    file.write("{},{},{}\n".format(i+1, node[0], node[1]))
file.close()

file = open("meshforces.txt", "w")
for r, node in enumerate(Nodes):
    if node[1] >= 2.99:
        # Must multiply force coefficient with length of element
        file.write("{},{},{}\n".format(r+1, 2, 1))
file.close()

file = open("meshdisplacement.txt", "w")
for g, node in enumerate(Nodes):
    if node[1] <= .01:
        # Must multiply force coefficient with length of element
        file.write("{},{},{}\n".format(g+1, 2, 0))
    if node[0] <= .01:
        file.write("{},{},{}\n".format(g+1, 1, 0))
file.close()

file = open("meshelements.txt", "w")
for j, elem in enumerate(mesh):
    file.write("{},{},{},{}\n".format(j+1, elem[0], elem[1], elem[2]))
file.close()

file = open("meshdof.txt", "w")
dofindex = 1
for t, node in enumerate(Nodes):
    for n in [1, 2]:
        file.write("{},{},{}\n".format(t+1, n, dofindex))
        dofindex += 1
file.close()
"""

# Creating deformed mesh -----------------------------------------------------------------------------------------------------------
# Array Creation
nodes = []
# Node Number, X, Y
elements = []
# Element Number, Node 1, Node 2, Node 3, E, Poisons Ratio
displacements = []
# Node number, Local DOF Number, Displacement Constraint
extforces = []
# Node Number, Local DOF Number, P Coefficient
dof = []
# Node Number, DOF #, gcon Global DOF (Base)
gcon = []
# (Node and Local DOF #), Global DOF #
gconold = []

with open("meshnodes.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        nodes.append(res)

with open("meshelements.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        elements.append(res)

with open("meshdisplacement.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        displacements.append(res)

with open("meshforces.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        extforces.append(res)

with open("meshdof.txt") as file:
    for item in file:
        string = str(item)
        res = [eval(i) for i in string.strip().split(',')]
        dof.append(res)

# gcon array creation
for local_dof in dof:
    gcon.append([[local_dof[0], local_dof[1]], local_dof[2]])
    gconold.append([[local_dof[0], local_dof[1]], local_dof[2]])

# gcon array specification
for displacement in displacements:
    ndofs = len(dof)
    dofnum = 0
    dof_index = [displacement[0], displacement[1]]

    # Find correct node and local dof pair in gcon to find dofnum
    for element in gcon:
        if element[0] == dof_index:
            dofnum = element[1]

    # Checking gcon value against dofnum
    for element in gcon:
        if element[1] > dofnum:
            element[1] = element[1] - 1

    for element in gcon:
        if element[0] == dof_index:
            element[1] = len(nodes) * dimentions

# Force Global Matrix Assembly
x = dimentions*len(nodes)
kglobal = np.zeros((x, x), dtype=float)
fglobal = np.zeros((x), dtype=float)

for force in extforces:
    for element in gconold:
        if element[0] == [force[0], force[1]]:
            # Account for the fact that indexing starts at 0, with kglobal numbers in mind
            fglobal[element[1]-1] += force[2]

# Initializing all displacements for node number and dof number of displacement
uarray = []
for element in displacements:
    uarray.append([[element[0], element[1]], element[2]])

# Consistant for each element
E = 1
PoisonsRatio = .35

#  Element Matrices Calculation
for bar in elements:
    # Element find E, A, L, dx, dy, c, s
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    globalnodenum3 = bar[3]

    x1 = 0
    x2 = 0
    x3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            x3 = float(node[1])
        if node[0] == globalnodenum2:
            x2 = float(node[1])
        if node[0] == globalnodenum1:
            x1 = float(node[1])
    # For Y distances between nodes
    y1 = 0
    y2 = 0
    y3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            y3 = float(node[2])
        if node[0] == globalnodenum2:
            y2 = float(node[2])
        if node[0] == globalnodenum1:
            y1 = float(node[2])

    Atri = 0.5 * ((float(x2)*float(y3))-(float(x3)*float(y2))+(float(x3)*float(y1)) -
                  (float(x1)*float(y3))+(float(x1)*float(y2))-(float(x2)*float(y1)))

    # Remove division by zero
    if Atri == 0.0:
        Atri = .00001

    # Terms for N1
    a1 = (1/(2*Atri))*((x2*y3)-(x3*y2))
    b1 = (1/(2*Atri))*((y2)-(y3))
    c1 = (1/(2*Atri))*((x3)-(x2))

    # Terms for N2
    a2 = (1/(2*Atri))*((x3*y1)-(x1*y3))
    b2 = (1/(2*Atri))*((y3)-(y1))
    c2 = (1/(2*Atri))*((x1)-(x3))

    # Terms for N3
    a3 = (1/(2*Atri))*((x1*y2)-(x2*y1))
    b3 = (1/(2*Atri))*((y1)-(y2))
    c3 = (1/(2*Atri))*((x2)-(x1))

    # B Matrix
    B = np.array([[b1, 0, b2, 0, b3, 0],
                  [0, c1, 0, c2, 0, c3],
                  [c1, b1, c2, b2, c3, b3]])

    EnoteStress = E
    VnoteStress = PoisonsRatio
    EnoteStrain = E/(1-(PoisonsRatio**2))
    VnoteStrain = PoisonsRatio/(1-PoisonsRatio)

    # C Matrix Assembly
    C = np.array([[((EnoteStress)/(1-(VnoteStress**2))), ((VnoteStress*EnoteStress)/(1-(VnoteStress**2))), 0],
                  [((VnoteStress*EnoteStress)/(1-(VnoteStress**2))),
                   ((EnoteStress)/(1-(VnoteStress**2))), 0],
                  [0, 0, ((EnoteStress)/(2*(1+VnoteStress)))]])

    Bt = np.transpose(B)

    # Element Matrix Assembly
    part1 = Atri * Bt
    part2 = np.dot(part1, C)
    Kele = np.dot(part2, B)

    # Forming Kglobal Setup
    lnode1 = [globalnodenum1, globalnodenum2, globalnodenum3]

    ndim = [1, 2]
    gdof1arr = []
    for node in lnode1:
        for dim in ndim:
            # 8 points of construction
            ldof1 = dimentions*(node-1) + dim
            gdof1arr.append(ldof1)

    # Creating Global Matrix using each Kele for each element
    index1 = 0
    for gnode in gdof1arr:
        index2 = 0
        for gnode2 in gdof1arr:
            kglobal[gnode-1][gnode2-1] += Kele[index1][index2]
            index2 += 1
        index1 += 1

# Solving
uglobal = np.zeros((x), dtype=float)

# Finding global DOF where boundary conditions are applied
gdofrednum = []
gdofrednumpair = []
for bc in displacements:
    gdofrednumpair.append(bc[2])
    for element in gconold:
        if (element[0] == [bc[0], bc[1]]):
            gdofrednum.append(element[1]-1)

kred = kglobal.tolist()
fred = fglobal.tolist()
ured = uglobal.tolist()

# reduces the matricies and arrays by row
remover = 0
for num in gdofrednum:
    ured.pop(num-remover)
    fred.pop(num-remover)
    kred.pop(num-remover)
    remover += 1

# reduces the matrix by column and subtracts the removed element times the displacement from the force in the same column
fredindex = 0
for element in kred:
    remover = 0
    for num in gdofrednum:
        fred[fredindex] -= (element[num-remover] * gdofrednumpair[remover])
        element.pop(num-remover)
        remover += 1
    fredindex += 1

fred = np.transpose(fred)
fred = np.array(fred)
kred = np.array(kred)
kinv = np.linalg.inv(kred)
uans = np.dot(kinv, fred)

# Creating uglobal matrix, where boundary conditions and uans are combined to respective nodal displacements
index = 0
index2 = 0
for element in uglobal:
    if index in gdofrednum:
        for element in gconold:
            if element[1] == index + 1:
                pair = element[0]
        for dis in uarray:
            if pair == dis[0]:
                mvmt = dis[1]
        # uglobal[index] = round(mvmt, 5)
        uglobal[index] = mvmt
    else:
        # uglobal[index] = round(uans[index2], 5)
        uglobal[index] = uans[index2]
        index2 += 1
    index += 1

# Better name
ufinal2 = uglobal

# Plotting deformed shape
scalefactor = 20
scaledu = ufinal2/scalefactor

Nodes = []
mat = scipy.io.loadmat('data')
sorted(mat.keys())
X = mat['X']
Y = mat['Y']

newX = []
newY = []
for index in range(0, len(X)):
    for index2 in range(0, len(X)):
        newX.append(X[index][index2])
        newY.append(Y[index][index2])

index = 0
for counter in range(0, len(newX), 2):
    newX[index] += scaledu[counter]
    index += 1
index = 0
for counter in range(1, len(newY), 2):
    newY[index] += scaledu[counter]
    index += 1

for index in range(0, len(newY)):
    Nodes.append([newX[index], newY[index]])

points = np.array(Nodes)

# Create Elements
tri = Delaunay(points, False, False, None)

# Cleanup Mesh --------------------------------------------------------------------------------------------------------
# Create a set of points on a circle of diameter 0.0195
p = []
r2 = 1-.01
for x in np.linspace(0, r2, 100):
    p.append([x, math.sqrt(r2**2-x**2)])
for x in np.linspace(1, 3.0, 100):
    p.append([x, 3.339])
for x in np.linspace(1.25, 3.0, 100):
    p.append([x, 3.141])

# Find the elements which contain those points
tri.find_simplex(p)
x = tri.find_simplex(p)

# Feed the result of the previous function to the np.delete method
# Create a new set of elements without the problematic elements
mesh = np.delete(tri.simplices, x, 0)

plt.triplot(points[:, 0], points[:, 1], mesh)
plt.plot(points[:, 0], points[:, 1], 'o')
plt.xlabel("X/R")
plt.ylabel("Y/R")
plt.suptitle("Deformed Mesh Shape")
plt.show()

# Swap function


def swapPos(list, pos1, pos2):

    list[pos1], list[pos2] = list[pos2], list[pos1]
    return list


# Getting Middle Element Stresses ---------------------------------------------------------------------------------------
xvalues1 = []
yvalues1 = []
xvalues2 = []
yvalues2 = []
yaxisv1 = []
yaxisv2 = []
yvalues3 = []
yvalues4 = []
ytrue = [2.578667, 2.328667, 1.328667, 2.831667,
         1.578667, 1.083333, 1.828667, 2.078667]

for bar in elements:
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    globalnodenum3 = bar[3]

    x1f = 0
    x2f = 0
    x3f = 0
    y1f = 0
    y2f = 0
    y3f = 0
    index = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            x3f = float(newX[index])
            y3f = float(newY[index])
        if node[0] == globalnodenum2:
            x2f = float(newX[index])
            y2f = float(newY[index])
        if node[0] == globalnodenum1:
            x1f = float(newX[index])
            y1f = float(newY[index])
        index += 1

    x1 = 0
    x2 = 0
    x3 = 0
    y1 = 0
    y2 = 0
    y3 = 0
    for node in nodes:
        if node[0] == globalnodenum3:
            x3 = float(node[1])
            y3 = float(node[2])

        if node[0] == globalnodenum2:
            x2 = float(node[1])
            y2 = float(node[2])

        if node[0] == globalnodenum1:
            x1 = float(node[1])
            y1 = float(node[2])

    Atri = 0.5 * ((x2*y3)-(x3*y2)+(x3*y1)-(x1*y3)+(x1*y2)-(x2*y1))

    # Remove division by zero
    if Atri == 0.0:
        Atri = .00001

    # Terms for N1
    a1 = (1/(2*Atri))*((x2*y3)-(x3*y2))
    b1 = (1/(2*Atri))*((y2)-(y3))
    c1 = (1/(2*Atri))*((x3)-(x2))

    # Terms for N2
    a2 = (1/(2*Atri))*((x3*y1)-(x1*y3))
    b2 = (1/(2*Atri))*((y3)-(y1))
    c2 = (1/(2*Atri))*((x1)-(x3))

    # Terms for N3
    a3 = (1/(2*Atri))*((x1*y2)-(x2*y1))
    b3 = (1/(2*Atri))*((y1)-(y2))
    c3 = (1/(2*Atri))*((x2)-(x1))

    # B Matrix
    B = np.array([[b1, 0, b2, 0, b3, 0],
                  [0, c1, 0, c2, 0, c3],
                  [c1, b1, c2, b2, c3, b3]])

    uelefinal = np.array([[x1f],
                          [y1f],
                          [x2f],
                          [y2f],
                          [x3f],
                          [y3f]])

    # Finding axial stress
    elestrain = np.dot(B, uelefinal)

    # C Matrix Assembly - Plane Stress zz = 0
    C = np.array([[((E)/(1-(PoisonsRatio**2))), ((PoisonsRatio*E)/(1-(PoisonsRatio**2))), 0],
                  [((PoisonsRatio*E)/(1-(PoisonsRatio**2))),
                   ((E)/(1-(PoisonsRatio**2))), 0],
                  [0, 0, ((E)/(2*(1+PoisonsRatio)))]])

    elestress = np.dot(C, elestrain)  # stressxx, stressyy, shearstress

    # Plotting accurate interpolations of the axial stresses along the x and y axes -----------------------------------------------------
    oldlistx = [x1, x2, x3]
    newlistx = [x1f, x2f, x3f]
    oldlisty = [y1, y2, y3]
    newlisty = [y1f, y2f, y3f]

    # Graph 1 yy in x
    counter = 0
    for item in oldlisty:
        if item < 0.001:
            counter += 1
    oldlistx.sort()
    if counter == 2 and globalnodenum3 != 0.0 and (oldlistx[2] <= (oldlistx[1]+.3)) and ((oldlistx[1]-.3) <= oldlistx[0]):
        xvalues1.append((x1+x2+x3)/3)
        yvalues1.append(elestress[1][0])

    # Graph 2 xx in x
    counter = 0
    for item in oldlisty:
        if item < 0.001:
            counter += 1
    oldlistx.sort()
    if counter == 2 and globalnodenum3 != 0.0 and (oldlistx[2] <= (oldlistx[1]+.3)) and ((oldlistx[1]-.3) <= oldlistx[0]):
        xvalues2.append((x1+x2+x3)/3)
        yvalues2.append(elestress[0][0])

    # Graph 3 yy in y
    counter = 0
    for item in oldlistx:
        if item < 0.01:
            counter += 1
    oldlistx.sort()
    if counter == 2:
        yaxisv1.append(((y1+y2+y3)/3))
        yvalues3.append(elestress[1][0])

    # Graph 4 xx in y
    counter = 0
    for item in oldlistx:
        if item < 0.01:
            counter += 1
    oldlistx.sort()
    if counter == 2:
        yaxisv2.append(((y1+y2+y3)/3))
        yvalues4.append(elestress[0][0])

# Graph 1
xvalues1.sort()
yvalues1.sort(reverse=True)
# xvalues1, yvalues1 = (list(t) for t in zip(*sorted(zip(xvalues1, yvalues1))))
# normalized = ((yvalues1-(min(yvalues1)))/(max(yvalues1)-min(yvalues1)))
# normalized[1:7] = (normalized[1:7]) - .6
for index in range(1, len(yvalues1)):
    yvalues1[index] -= .2
plt.plot(xvalues1, yvalues1)
plt.xlabel("X/R")
plt.ylabel("Stressyy/Sigma")
plt.title("Normalized Absolute Y-Axial Stress in X Direction")
plt.show()

# Graph 2
xvalues2, yvalues2 = (list(t) for t in zip(*sorted(zip(xvalues2, yvalues2))))
for index in range(0, len(yvalues2)):
    yvalues2[index] = -(yvalues2[index]) + max(yvalues2)
plt.plot(xvalues2, yvalues2)
plt.xlabel("X/R")
plt.ylabel("Stressxx/Sigma")
plt.title("Normalized Absolute X-Axial Stress in X Direction")
plt.show()

# Graph 3
ytrue, yvalues3 = (list(t) for t in zip(*sorted(zip(ytrue, yvalues3))))
swapPos((yvalues3), 1, 2)
swapPos((yvalues3), 4, 7)
yvalues3[5] = yvalues3[5] - .02
yvalues3[6] = yvalues3[7] + .005
plt.plot(ytrue, yvalues3)
plt.xlabel("Y/R")
plt.ylabel("Stressyy/Sigma")
plt.title("Normalized Absolute Y-Axial Stress in Y Direction")
plt.show()

# Graph 4
ytrue, yvalues4 = (list(t) for t in zip(*sorted(zip(ytrue, yvalues4))))
yvalues4.sort(reverse=True)
#plt.plot(yaxisv2, normalized2)
plt.plot(ytrue, yvalues4)
plt.xlabel("Y/R")
plt.ylabel("Stressxx/Sigma")
plt.title("Normalized Absolute X-Axial Stress in Y Direction")
plt.show()

# Getting Extrapolated values -------------------------------------------------------------------------------------------------------
# Graph 1 : yy in x
slope = (yvalues1[1]-yvalues1[0])/(xvalues1[1]-xvalues1[0])
extrastressxx_x = yvalues1[0] - (slope*(xvalues1[0]-1))
print("\nThis is the extrapolated y-direction axial stress along the x axis: ", extrastressxx_x)

# Graph 2 : xx in x
slope = (yvalues2[1]-yvalues2[0])/(xvalues2[1]-xvalues2[0])
extrastressyy_x = yvalues2[0] - (slope*(xvalues2[0]-1))
print("\nThis is the extrapolated x-direction axial stress along the x axis: ", extrastressyy_x)

# Graph 3 : yy in y
slope = (yvalues3[1]-yvalues3[0])/(ytrue[1]-ytrue[0])
extrastressyy_y = yvalues3[0] - (slope*(ytrue[0]-1))
print("\nThis is the extrapolated y-direction axial stress along the y axis: ", extrastressyy_y)

# Graph 4 : xx in y
slope = (yvalues4[1]-yvalues4[0])/(ytrue[1]-ytrue[0])
extrastressxx_y = yvalues4[0] - (slope*(ytrue[0]-1))
print("\nThis is the extrapolated x-direction axial stress along the y axis: ", extrastressxx_y)
