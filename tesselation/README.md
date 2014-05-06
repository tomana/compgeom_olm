## Greedy tesselation algorithm for triangulating triangles with respect to multiple intersecting facets.

Novel implementation of an algorithm _find references_, extended to handle 
cases where more than one facet is intersecting a triangle. Using [Dolfin](https://bitbucket.org/fenics-project/dolfin) datastructures.

### Short description:

#### Classes:
_TesselateMeshIntersect.cpp/h_ is the main class tesselating the partially overlapped elements of one mesh with respect to the boundary of another.

_TesselateTriangle.cpp/h_ is the class triangulating each individual triangle.

#### Functions:

init(mesh, collidingMesh) takes the two input meshes, _mesh_ beign the _background mesh_ and _collidingMesh_ being the overlapping mesh.

```cpp triangulate_first(...) ``` calls a collision detecting algorithm to fill the map
```cpp
boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirstMesh
```
associating indexes in the first input _mesh_ to a vector of indexes with the colliding entities in _collidingMesh_.


