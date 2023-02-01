import mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
import numpy, scipy.io

path = "/media/rasoul/6aa24340-1780-4aae-8b9f-ffb7e2f3402c/Partial Nephrectomy/Codes/Geometry_Processing/Adrien_Idea/My_Optimisation/"
# path = "/media/rasoul/6aa24340-1780-4aae-8b9f-ffb7e2f3402c/Partial Nephrectomy/Codes/Geometry_Processing/Adrien_Idea/least-squares-course-master/src/ch6/5-lscm/"
m = mesh.Mesh(path + "test10.obj") # test10.obj load mesh

np.savetxt('F.out', m.T, delimiter=',')   # X is an arra
A = scipy.sparse.lil_matrix((2*m.ntriangles+4, 2*m.nverts)) # the variables are packed as u0,v0,u1,v1, ...
lock1, lock2 = 1,2#10324%m.nverts, 35492%m.nverts  # select two arbitrary vertices to pin
for (t,[i,j,k]) in enumerate(m.T): # for each triangle ijk
    [eij, ejk, eki] = mesh.project_triangle(m.V[i], m.V[j], m.V[k]) # project the triangle to a local 2d basis
    A[t*2+0, i*2  ] =  ejk[0] # (grad u)[0] = (grad v)[1]
    A[t*2+0, j*2  ] =  eki[0]
    A[t*2+0, k*2  ] =  eij[0]
    A[t*2+0, i*2+1] =  ejk[1]
    A[t*2+0, j*2+1] =  eki[1]
    A[t*2+0, k*2+1] =  eij[1]
    A[t*2+1, i*2  ] =  ejk[1] # (grad u)[1] = -(grad v)[0]
    A[t*2+1, j*2  ] =  eki[1]
    A[t*2+1, k*2  ] =  eij[1]
    A[t*2+1, i*2+1] = -ejk[0]
    A[t*2+1, j*2+1] = -eki[0]
    A[t*2+1, k*2+1] = -eij[0]
A[-1,lock2*2+1] = A[-2,lock2*2+0] = A[-3,lock1*2+1] = A[-4,lock1*2+0] = 1 # quadratic penalty
A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications

b = [0]*(2*m.ntriangles) + [0,0,1,1] # one pinned to (0,0), another to (1,1)
x, istop, itn, normr = lsmr(A, b)[0:4] # call the least squares solver

for v in range(m.nverts): # apply the computed flattening
    m.V[v] = np.array([x[v*2], x[v*2+1], 0])
print(m) # output the deformed mesh

filepath = path + "Out10.txt"

with open(filepath, 'w') as f:
    # f.write("# OBJ file\n")
    for v in range(m.nverts):
        f.write("v %.4f %.4f %.4f\n" % (m.V[v,0], m.V[v,1] , m.V[v,2]))
    for i in range(m.ntriangles):
        # f.write("f")
        # for i in range(m.ntriangles):
        f.write("f %.4f %.4f %.4f\n" % (m.T[i,0]+1, m.T[i,1] +1, m.T[i,2]+1))
    f.write("\n")
B = A.todense()
Mat_path = '/media/rasoul/6aa24340-1780-4aae-8b9f-ffb7e2f3402c/Partial Nephrectomy/Codes/Geometry_Processing/Adrien_Idea/My_Optimisation/A.mat';
scipy.io.savemat( Mat_path , mdict={'arr': B})

# C = b.todense()
Mat_path = '/media/rasoul/6aa24340-1780-4aae-8b9f-ffb7e2f3402c/Partial Nephrectomy/Codes/Geometry_Processing/Adrien_Idea/My_Optimisation/b.mat';
scipy.io.savemat( Mat_path , mdict={'arr': b})