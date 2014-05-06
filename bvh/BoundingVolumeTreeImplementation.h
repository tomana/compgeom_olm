// Copyright (C) 2009 Andre Massing
//
// This file is part of DOLFIN.

// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Johannes Ring, 2009.
//
// First added:  2009-09-11
// Last changed: 2013-09-16

#ifndef BOUNDINGVOLUMETREEIMPLEMENTATION_H
#define BOUNDINGVOLUMETREEIMPLEMENTATION_H

#include <vector>
#include <utility>
#include <queue>
#include <stack>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>

#include <dolfin/geometry/Point.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/SubsetIterator.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <ctime>
#include <set>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/NoDeleter.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#define HAS_CGAL 1

#ifdef HAS_CGAL

#include "cutfem_geometry.h"
#include "cgal_includes.h"

typedef CGAL::Simple_cartesian<double> SCK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;

namespace dolfin
{

/// @brief Interface class for the actual implementation of the intersection operations.
///
/// @internal This is neccessary since the search tree has a dimension dependent type, hence encapsulates in the
/// inheriting BoundingVolumeTreeImplementation_d! It provides the glue level between the dimension independent implementation
/// of the mesh class and the dimension dependent search structures in CGAL.
class BoundingVolumeTreeImplementation
{
public:

    // Only default constructor, since the search tree has a dimension dependent type, hence encapsulates in the
    // inheriting BoundingVolumeTreeImplementation_d!

    // Virtual destructor
    virtual ~BoundingVolumeTreeImplementation() {}

    virtual void all_intersected_entities(const Point & point,  boost::unordered_map<size_t, std::vector< size_t> > & ids_result) const = 0;
    virtual void all_intersected_entities(const std::vector<Point> & points,  boost::unordered_map<size_t, std::vector< size_t> > & ids_result) const = 0;

    virtual void all_intersected_entities(const MeshEntity & entity,  boost::unordered_map<size_t, std::vector< size_t> > & ids_result) const = 0;
    virtual void all_intersected_entities(const std::vector<MeshEntity> & entities,  boost::unordered_map<size_t, std::vector< size_t> > & ids_result) const = 0;

    virtual void all_intersected_entities(const Mesh & another_mesh,  boost::unordered_map<size_t, std::vector< size_t> > & ids_result) const = 0;
    virtual int any_intersected_entity(const Point & point) const = 0;
    virtual Point closest_point(const Point & point) const = 0;
    virtual std::size_t closest_cell(const Point & point) const = 0;
    virtual std::pair<Point,std::size_t> closest_point_and_cell(const Point & point) const = 0;
    virtual double distance(const Point & point) const = 0;
};


/// MY IMPLEMENTATION
class BoundingVolumeTreeImplementation_d : public BoundingVolumeTreeImplementation
{
public:

    /// Constructor

    BoundingVolumeTreeImplementation_d(std::shared_ptr<const Mesh> mesh, std::string bvtype, std::string kerneltype, bool usedop);

    BoundingVolumeTreeImplementation_d(const MeshFunction<std::size_t> labels, std::size_t label);

    /// Destructor which ensures to reset AABB tree
    virtual ~BoundingVolumeTreeImplementation_d();

    virtual void all_intersected_entities(const Point& point,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const;
    virtual void all_intersected_entities(const std::vector<Point>& points,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const;

    virtual void all_intersected_entities(const MeshEntity& entity,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const;
    virtual void all_intersected_entities(const std::vector<MeshEntity>& entities,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const;

    virtual void all_intersected_entities(const Mesh& another_mesh,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const;

    virtual  int any_intersected_entity(const Point& point) const;

    virtual Point closest_point(const Point& point) const;
    virtual std::size_t closest_cell(const Point& point) const;
    virtual std::pair<Point, std::size_t> closest_point_and_cell(const Point& point) const;
    virtual double distance(const Point & point) const;

    ///Topological dimension of the mesh.
    //static const std::size_t dim = PT::dim;
    mutable bool usedopbool;

private:

    mutable bool point_search_tree_constructed;
    mutable std::string bvtypestring = "BoundingSphere";
    mutable std::string kerneltypestring = "SCK";


    /* Tree
    *  RTCD OPTIMALISERING kap
    *
    */
    template<size_t N>
    struct KDOP
    {
        float dist_[N];
        float r;
    };


    template <class T>
    struct Node
    {
        T BV;
        KDOP<26> LEAFBV;
        int type;
        int numObjects;
        Point midpoint;
        size_t object;
        Node *left;
        Node *right;
    };

    /* Bounding Volume Sphere
    *
    *
    */
    struct BVsphere
    {
        Point c;
        float r;
    };

    int testBVBV(BVsphere &a, BVsphere &b) const;
    void mostSeparatedPointsOnAABB(int &minimum, int &maximum, Point points[], int numPts) const;
    void sphereFromDistantPoints(BVsphere &s, Point points[], int numPts) const;
    void sphereOfSphereAndPoint(BVsphere &s, Point &p) const;
    void ritterSphere(BVsphere &s, Point points[], int numPts) const;

    void computeBV(BVsphere &sphere, std::vector< Point > v) const;

    bool DescendA(Node<BVsphere> **a, Node<BVsphere> **b) const
    {
    return ((**b).type == 0) || (!((**a).type == 0) && (**a).BV.r >= (**b).BV.r);
    }

    /* AABB
    *
    *
    */
    struct AABB
    {
        Point xmin;
        Point xmax;
        float r;
    };

    int testBVBV(AABB &a, AABB &b) const;

    void computeBV(AABB &s, std::vector< Point > v) const;

    bool DescendA(Node<AABB> **a, Node<AABB> **b) const
    {
    return ((**b).type == 0) || (!((**a).type == 0) && (**a).BV.r >= (**b).BV.r);
    }
    /* k Discrete Oriented Polytope
    *
    *
    */

    template<std::size_t N>
    int testBVBV(KDOP<N> &a, KDOP<N> &b) const;

    void computeBV(KDOP<6> &dop6,  std::vector< Point > v) const;
    void computeBV(KDOP<8> &dop8,  std::vector< Point > v) const;
    void computeBV(KDOP<12> &dop12,  std::vector< Point > v) const;
    void computeBV(KDOP<14> &dop14,  std::vector< Point > v) const;
    void computeBV(KDOP<18> &dop18,  std::vector< Point > v) const;
    void computeBV(KDOP<26> &dop26,  std::vector< Point > v) const;
    void evalsize(float value, float &smallest, float &biggest) const;

    // MERGETEST

    void mergeNodes(BVsphere &BV, Node<BVsphere> nodes[], int numObjects) const;
    void mergeNodes(AABB &BV, Node<AABB> nodes[], int numObjects) const;
    template<std::size_t N>
    void mergeNodes(KDOP<N> &BV, Node<KDOP<N> > nodes[], int numObjects) const;

    template<size_t N>
    bool DescendA(Node<KDOP<N> > **a, Node<KDOP<N> > **b) const
    {
    return ((**b).type == 0) || (!((**a).type == 0) && (**a).BV.r >= (**b).BV.r);
    }

    /*
    bool DescendA(Node<KDOP<18> > **a, Node<KDOP<18> > **b) const
    {
    return ((**b).type == 0) || (!((**a).type == 0) && (**a).BV.r >= (**b).BV.r);
    }
    bool DescendA(Node<KDOP<12> > **a, Node<KDOP<12> > **b) const
    {
    return ((**b).type == 0) || (!((**a).type == 0) && (**a).BV.r >= (**b).BV.r);
    }
    */

    Node<BVsphere> **bvspheretree;
    Node<AABB> **aabbtree;
    Node<KDOP<6> > **kdop6tree;
    Node<KDOP<8> > **kdop8tree;
    Node<KDOP<12> > **kdop12tree;
    Node<KDOP<14> > **kdop14tree;
    Node<KDOP<18> > **kdop18tree;
    Node<KDOP<26> > **kdop26tree;

    Mesh _mesh;

    mutable int counterqueries;
    mutable int countertraversals;
    mutable int countertrue;

    template <class T>
    void topDownBVTreeFromMesh(Node<T> **tree, const Mesh &mesh) const;

    template <class T>
    void topDownBVTreeFromMeshEntities(Node<T> **tree, const std::vector<MeshEntity>& entities) const;

    template <class T>
    void topDownBVTreeFromPoints(Node<T> **tree, const std::vector<Point>& points) const;

    template <class T, class Container>
    void topDownBVTree(Node<T> **tree, const Container &container, size_t object[], int tdim, int numObjects) const;

    template <class T, class Container>
    void topDownBVTree(Node<T> **tree, const Container &container, std::vector <T> boundingvolumes, int tdim, int numObjects) const;

    template <class T, class Container>
    void topDownBVTree(Node<T> **tree, const Container &container, Node<T> initialnodes[], int numObjects) const;

    std::vector<Point> pointsFromContainer(const Mesh & mesh, size_t object[], int numObjects) const;

    std::vector<Point> pointsFromContainer(const std::vector<MeshEntity> & entities, size_t object[], int numObjects) const;

    std::vector<Point> pointsFromContainer(const std::vector<Point> & points, size_t object[], int numObjects) const;

    std::vector<Point> pointsFromContainer(const MeshEntity& entity) const;

    int partitionObjects(size_t * object, const Mesh & mesh, int numObjects, int axis) const;
    int partitionObjects(size_t * object, const std::vector<MeshEntity> & entities, int numObjects, int axis) const;
    int partitionObjects(size_t * object, const std::vector<Point> & points, int numObjects, int axis) const;

    template <class T>
    int partitionNodes(Node<T> initialnodes[], int numObjects, int axis) const;

    void mostSeparatedAxisOnAABB(int &axis, Point points[], int numPts) const;

    template <class T>
    void mostSeparatedMidpointsOfNodes(int &axis, Node<T> initialnodes[], int numObjects) const;

    template <class T>
    void midPointOfNodes(Point &midpoint, Node<T> initialnodes[], int numObjects) const;

    template<class T>
    void treeCollidingWithTree(std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > r, const Mesh &mesh, Node<T> **a, Node< T> **b, std::string kernel) const;

    template<class T, class Container >
    void treeCollidingWithTreeNonRecursive(std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > r, const Container &container, Node<T> **a, Node<T> **b, std::string kernel, bool usedop26) const;

    template<class T, class Container >
    void treeCollidingWithTreeNonRecursive( boost::unordered_map<size_t, std::vector< size_t> > &r, const Container &container, Node<T> **a, Node<T> **b, std::string kernel, bool usedop26) const;


    template<class T, class Container >
    void treeCollidingWithEntity(std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > r, const Container &container, Node<T> **a, MeshEntity b, std::string kernel) const;

    template<class T>
    bool BVOverlap(Node<T> **a, Node<T> **b) const;


    struct sortmeshx
    {
      const Mesh & mesh;

      sortmeshx(const Mesh & mesh): mesh(mesh) {}

      inline bool operator()(size_t i, size_t j)
      {
        dolfin::Cell pickcell(mesh, i);
        dolfin::VertexIterator v(pickcell);

        dolfin::Cell pickcell2(mesh, j);
        dolfin::VertexIterator v2(pickcell2);

        return (v[0].point().x() < v2[0].point().x());
      }
    };

    struct sortmeshy
    {
      const Mesh & mesh;

      sortmeshy(const Mesh & mesh): mesh(mesh) {}

      inline bool operator()(size_t i, size_t j)
      {
        dolfin::Cell pickcell(mesh, i);
        dolfin::VertexIterator v(pickcell);

        dolfin::Cell pickcell2(mesh, j);
        dolfin::VertexIterator v2(pickcell2);

        return (v[0].point().y() < v2[0].point().y());
      }
    };

    struct sortmeshz
    {
      const Mesh & mesh;

      sortmeshz(const Mesh & mesh): mesh(mesh) {}

      inline bool operator()(size_t i, size_t j)
      {
        dolfin::Cell pickcell(mesh, i);
        dolfin::VertexIterator v(pickcell);

        dolfin::Cell pickcell2(mesh, j);
        dolfin::VertexIterator v2(pickcell2);

        return (v[0].point().z() < v2[0].point().z());

      }
    };

    template <class T>
    struct sortkdopx
    {
      Node<T> * initialnodes;
      sortkdopx<T>(Node<T> *initialnodes): initialnodes(initialnodes) {}
      inline bool operator()(Node<T> i, Node<T> j)
      {

        return (i.midpoint.x() < j.midpoint.x());
      }
    };

    template <class T>
    struct sortkdopy
    {
      Node<T> * initialnodes;
      sortkdopy<T>(Node<T> *initialnodes): initialnodes(initialnodes) {}
      inline bool operator()(Node<T> i, Node<T> j)
      {

        return (i.midpoint.y() < j.midpoint.y());
      }
    };

    template <class T>
    struct sortkdopz
    {
      Node<T> * initialnodes;
      sortkdopz<T>(Node<T> *initialnodes): initialnodes(initialnodes) {}
      inline bool operator()(Node<T> i, Node<T> j)
      {

        return (i.midpoint.z() < j.midpoint.z());
      }
    };

    template<class T>
    bool CollidePrimitivesEPICK(const Mesh &mesh, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesEPICK(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesEPICK(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const;

    template<class T>
    bool CollidePrimitivesSCK(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesSCK(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesSCK(const Mesh &mesh, Node<T> **a, Node<T> **b) const;

    template<class T>
    bool CollidePrimitivesDOLFIN(const Mesh &mesh, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesDOLFIN(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const;
    template<class T>
    bool CollidePrimitivesDOLFIN(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const;
};
}

#endif

#endif
