## Bounding volume hierarchy with different bounding volumes for meshes.

Implementation of a Bounding Volume Hierarchy using [Dolfin](https://bitbucket.org/fenics-project/dolfin) datastructures and classes.

Many of the algorithms are based on excellent examples in the book [Real-time collision detecting](http://realtimecollisiondetection.net/) by Christer Ericson.

### Short description:
The tree construction algorithm available is top-down. The hierarchy is constructed by _merging_ bounding volumes of leaf nodes. A selection of bounding volumes are available, namely _Bounding Spheres_, _Axis Aligned Bounding Boxes_ (AABB) and _k Discrete Oriented Polytopes_ (kDOPs), all more or less capable of being merged in a satisfying manner.

#### Classes:
_BoundingVolumeTree.cpp/h_ is an interface class with virtual functions.

_BoundingVolumeTreeImplementation.cpp/h_ is the implementation containing the actual BVH.

#### Usage:
'''cpp
dolfin::Mesh * mesh;
dolfin::Mesh * collidingMesh;
boost::unordered\_map<size_t, std::vector< size_t > > intersectionset;

dolfin::BoundingVolumeTree * tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
tree->all_intersected_entities(*collidingMesh, intersectionset);
'''

The unordered\_map _intersectionset_ associates indexes in the first input _mesh_ to a vector of indexes with the colliding entities in _collidingMesh_.

Some options are available in the call to the constructor. In the above example "SimpleCartesian" denotes a primtive test (i.e. an triangle-triangle intersection test) from the [CGAL](www.cgal.org) computational geometry library, other available are "ExactPredicates", using the exact geometric predicates kernel in CGAL, and "Dolfin" using the primtive tests available in [Dolfin](https://bitbucket.org/fenics-project/dolfin). The argument "DOP26" specifies that a Discrete Oriented Polytopes with 26 planes is used, the other bounding volumes bounding volumes are "BoundingSphere", "AABB" and "DOPx" for x = 6, 8, 12, 14, 18, 26.
