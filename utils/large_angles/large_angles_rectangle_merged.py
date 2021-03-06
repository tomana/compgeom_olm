from dolfin import *
import numpy


#set_log_level(2);

Theta = pi/2
a, b = 0, 1.0
nr = 10  # divisions in r direction
nt = 50  # divisions in theta direction
mesh = RectangleMesh(a, 0, b, 1, nr, nt, 'crossed')


mesh = Mesh()
editor = DynamicMeshEditor()
editor.open(mesh, 2, 2, 2)
num_cells_x = 8;
num_cells_y = 40;

#editor.init_vertices(3*num_cells_x*num_cells_y) # test both versions of interface
#editor.init_cells(num_cells_x*num_cells_y)    # test both versions of interface

switchx = False;
switchy = False;

vert_list = []

for j in range(0,num_cells_y):
    switchy = not switchy;
    for i in range(0,num_cells_x):
# Add vertices
        switchx = not switchx;
        if(switchy): 
            switchx = not switchx
        spacer_x = 1.0/float(num_cells_x)
        spacer_y = 1.0/float(num_cells_y)
        if(switchx):
            if(i == 0):
                vert_list.append(((i+1)*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+2)*spacer_x, j*spacer_y))
            elif(i == num_cells_x - 1):
                vert_list.append((i*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, j*spacer_y))
            else:
                vert_list.append((i*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+2)*spacer_x, j*spacer_y))        
        else:
            if(i == 0):
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+2)*spacer_x, (j+1)*spacer_y))               
            elif(i == num_cells_x - 1):
                vert_list.append(((i)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
            else:
                vert_list.append(((i)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+2)*spacer_x, (j+1)*spacer_y))
        
        if(switchy):
            switchx = not switchx


print vert_list

set(vert_list)

new_list = list(set(vert_list))

print new_list

cell = 0;
for points in new_list:
    editor.add_vertex(cell, points[0], points[1])
    cell += 1

switchx = False;
switchy = False;

for j in range(0,num_cells_y):
    switchy = not switchy;
    for i in range(0,num_cells_x):
        vert_list = []
# Add vertices
        switchx = not switchx;
        if(switchy): 
            switchx = not switchx
        spacer_x = 1.0/float(num_cells_x)
        spacer_y = 1.0/float(num_cells_y)
        if(switchx):
            if(i == 0):
                vert_list.append(((i+1)*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+2)*spacer_x, j*spacer_y))
            elif(i == num_cells_x - 1):
                vert_list.append((i*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, j*spacer_y))
            else:
                vert_list.append((i*spacer_x, j*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+2)*spacer_x, j*spacer_y))                
        else:
            if(i == 0):
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+2)*spacer_x, (j+1)*spacer_y))                   
            elif(i == num_cells_x - 1):
                vert_list.append(((i)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j+1)*spacer_y))
            else:
                vert_list.append(((i)*spacer_x, (j+1)*spacer_y))
                vert_list.append(((i+1)*spacer_x, (j)*spacer_y))
                vert_list.append(((i+2)*spacer_x, (j+1)*spacer_y))
        
        if(switchy):
            switchx = not switchx
# Add cell
        
        editor.add_cell(j*num_cells_x + i, new_list.index(vert_list[0]), new_list.index(vert_list[1]), new_list.index(vert_list[2]))

# Close editor
editor.close()

mesh.order()
mesh.init()

print mesh.ordered()

file = File("mesh.xml")
file << mesh

plot(mesh)

V = FunctionSpace(mesh, "Lagrange", 1)

spacer_x = 1.0/float(num_cells_x)
spacer_y = 1.0/float(num_cells_y)
        
# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    #return x[0] < DOLFIN_EPS + spacer_x or x[0] > 1.0 - DOLFIN_EPS - spacer_x or x[1] < DOLFIN_EPS - spacer_y or x[1] > 1.0 - DOLFIN_EPS - spacer_y
    return x[0] < DOLFIN_EPS + spacer_x or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS  or x[1] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
g = Expression("sin(5*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution
plot(u)

interactive()
