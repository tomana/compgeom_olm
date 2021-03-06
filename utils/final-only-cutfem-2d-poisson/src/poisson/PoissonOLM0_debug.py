#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
cell = triangle

V = FiniteElement("CG", cell, 1)

u = TrialFunction(V)
v = TestFunction(V)
f = Coefficient(V)

n = FacetNormal(cell)
h = 2.0*Circumradius(cell)

dxq = dc(1, metadata={"num_cells": 1})

gamma = Coefficient(V)
g = Function(V)

# Standard assembly over uncut cells
a = inner(grad(u), grad(v))*dx(0)
# Nitsche term to enforce weak boundary conditions
# on the (matching) boundary.
a += -inner(dot(grad(u),n), v)*ds
a += -inner(dot(grad(v),n), u)*ds
a += gamma/h*inner(u,v)*ds

# Integrate same form over cut cells (marked with 1)
a += inner(grad(u), grad(v))*dxq

L = f*v*dx(0) + f*v*dxq 
