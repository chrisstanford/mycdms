#!/usr/bin/python3.5
import numpy as np
from vtk import vtkXMLUnstructuredGridReader
from vtk.util import numpy_support as VN

reader = vtkXMLUnstructuredGridReader()
reader.SetFileName('potential000000.vtu')
reader.Update()

data = reader.GetOutput()
print(data)

n_cells = data.GetNumberOfCells()
n_points = data.GetNumberOfPoints()
print(n_cells)
print(n_points)
points_per_cell = data.GetCell(0).GetNumberOfPoints()
cellpointids = np.empty(shape=(n_cells,points_per_cell),dtype='int')
xyzfield = np.empty(shape=(n_points,4),dtype='float')
 
for i in range(n_cells):
    cell = data.GetCell(i)
    pointids = [cell.GetPointId(p) for p in range(points_per_cell)]
    cellpointids[i] = np.array(pointids)

for i in range(n_points):
    point = data.GetPoint(i)
#    if i==0: print(data.GetPointData())
    field = data.GetPointData().GetArray("f_1004").GetTuple1(i)
    xyzfield[i] = [(point[0]-5)/1000, (point[1]-5)/1000, (point[2])/1000, field/100]

np.savetxt('Epot_tetra.txt',cellpointids,fmt='%d')
np.savetxt('Epot_points.txt',xyzfield,fmt='%.5f')





