## Front advancing algorithm for detecting the intersection between two meshes

Based on Martin Ganders original algorithm, described in the paper:

[An algorithm for non-matching grid projections with linear complexity](http://archive-ouverte.unige.ch/download/unige:6553/ATTACHMENT01)

Adapted from implementations in MeshKit 1.0 and MOAB.

[MeshKit 1.0](https://trac.mcs.anl.gov/projects/fathom/wiki/MeshKit) implementation: [MeshKit 1.0 ProjectShell.cpp](http://www.mcs.anl.gov/~fathom/meshkit-docs/html/ProjectShell_8cpp_source.html) 

[Mesh Oriented dAtaBase (MOAB)](http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB) implementation: [MOAB Intx2Mesh.cpp](http://www.mcs.anl.gov/~fathom/moab-docs/html/Intx2Mesh_8cpp_source.html)

### Core modifications

Adapted to [Dolfin](https://bitbucket.org/fenics-project/dolfin) datastructures and modified to attain linear complexity.
