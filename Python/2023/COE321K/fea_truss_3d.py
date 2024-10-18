# Aaron Pandian (ANP3238) COE321K
# The goal of this project is to determine nodal displacements of arbitrary 3-d truss

"""
The structure of the input files are as follows:
"Nodes.txt": node number, x coordinate, y coordinate, z coordinate; point of origin is the bottom left most node
"Nodal_Displacements.txt": node number, degree of displacement, and boundary condition or known diaplcement at node
"Elements.txt": element number, node 1, node 2, E, A
"Applied_Nodal_Force.txt": node number, degree of force, coefficient of P
"DOF.txt": node number, local degree of freedom, global degree of freedom number
"""

from array import array
from math import pi
import numpy as np
dimentions = 3

# Array Creation
nodes = []
# Node Number, X, Y (real value)
elements = []
# Element Number, Node 1, Node 2, E, A (actual values)
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

#  Element Matrices Calculation
for bar in elements:
    # Element find E, A, L, dx, dy, c, s
    E = float(bar[3])
    A = float(bar[4])
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    # For X distances between nodes
    nx1 = 0
    nx2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            nx2 = node[1]
        if node[0] == globalnodenum1:
            nx1 = node[1]
        dx = nx2 - nx1
    # For Y distances between nodes
    ny1 = 0
    ny2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            ny2 = node[2]
        if node[0] == globalnodenum1:
            ny1 = node[2]
        dy = ny2 - ny1

    # For Z distances
    nz1 = 0
    nz2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            nz2 = node[3]
        if node[0] == globalnodenum1:
            nz1 = node[3]
        dz = nz2 - nz1

    # Element matrix components
    L = ((dx)**2+(dy)**2+(dz)**2)**.5
    if dx == 0:
        c1 = 0.0
    else:
        c1 = dx/L
    if dy == 0:
        c2 = 0.0
    else:
        c2 = dy/L
    if dz == 0:
        c3 = 0.0
    else:
        c3 = dz/L

    # Computing element matrix
    Kele = np.array([[c1**2, c1*c2, c1*c3, -(c1**2), -(c1*c2), -(c1*c3)],
                     [c1*c2, c2**2, c2*c3, -(c1*c2), -(c2**2), -(c2*c3)],
                     [c1*c3, c2*c3, c3**2, -(c1*c3), -(c2*c3), -(c3**2)],
                     [-(c1**2), -(c1*c2), -(c1*c3), c1**2, c1*c2, c1*c3],
                     [-(c1*c2), -(c2**2), -(c2*c3), c1*c2, c2**2, c2*c3],
                     [-(c1*c3), -(c2*c3), -(c3**2), c1*c3, c2*c3, c3**2]])

    Kele = ((E*A)/L)*Kele

    # Creating Global Matrix using each Kele for each element
    lnode1 = [globalnodenum1, globalnodenum2]

    ndim = [1, 2, 3]
    gdof1arr = []
    for node in lnode1:
        for dim in ndim:
            # 8 points of construction
            ldof1 = dimentions*(node-1) + dim
            gdof1arr.append(ldof1)

    index1 = 0
    for gnode in gdof1arr:
        index2 = 0
        for gnode2 in gdof1arr:
            # print(round(Kele[index1][index2], 3))
            kglobal[gnode-1][gnode2-1] += round(Kele[index1][index2], 3)
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
strain = []
stress = []

for bar in elements:
    # Element find E, A, L, dx, dy, c, s
    E = float(bar[3])
    A = float(bar[4])
    globalnodenum1 = bar[1]
    globalnodenum2 = bar[2]
    # For X distances between nodes
    nx1 = 0
    nx2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            nx2 = node[1]
        if node[0] == globalnodenum1:
            nx1 = node[1]
        dx = nx2 - nx1
    # For Y distances between nodes
    ny1 = 0
    ny2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            ny2 = node[2]
        if node[0] == globalnodenum1:
            ny1 = node[2]
        dy = ny2 - ny1

    # For Z distances
    nz1 = 0
    nz2 = 0
    for node in nodes:
        if node[0] == globalnodenum2:
            nz2 = node[3]
        if node[0] == globalnodenum1:
            nz1 = node[3]
        dz = nz2 - nz1

    # Element matrix components
    L = ((dx)**2+(dy)**2+(dz)**2)**.5
    if dx == 0:
        c1 = 0.0
    else:
        c1 = dx/L
    if dy == 0:
        c2 = 0.0
    else:
        c2 = dy/L
    if dz == 0:
        c3 = 0.0
    else:
        c3 = dz/L

    # Strain: positive strech, and negative contraction
    for check in gconold:
        if [globalnodenum1, 1] == check[0]:
            ux1 = ufinal[check[1]-1]
        elif [globalnodenum1, 2] == check[0]:
            uy1 = ufinal[check[1]-1]
        elif [globalnodenum1, 3] == check[0]:
            uz1 = ufinal[check[1]-1]
        elif [globalnodenum2, 1] == check[0]:
            ux2 = ufinal[check[1]-1]
        elif [globalnodenum2, 2] == check[0]:
            uy2 = ufinal[check[1]-1]
        elif [globalnodenum2, 3] == check[0]:
            uz2 = ufinal[check[1]-1]

    elestrain = ((ux2-ux1)*(c1/L))+((uy2-uy1)*(c2/L))+((uz2-uz1)*(c3/L))
    # strain.append(round(elestrain, 4))
    strain.append(round(elestrain, 4))
    elestress = elestrain * E
    # stress.append(round(elestress, 4))
    stress.append(round(elestress, 4))

ones = np.ones(len(elements), dtype=float)
onesmatrix = (np.diag(ones))*(E*A)
fbar = np.dot(onesmatrix, strain)

print(strain)
