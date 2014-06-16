#pragma once

#include <string>
#include <utility>
#include <boost/scoped_ptr.hpp>
#include <dolfin/common/types.h>

#include "BoundingVolumeTree.h"
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/function/dolfin_function.h>
#include <cutfem/mesh/CutMesh.h>
#include "TesselateTriangle.h"
#include <dolfin/refinement/refine.h>

#define HAS_CGAL 1

#include "cgal_includes.h"
#include "cutfem_geometry.h"

class TesselateMeshIntersect
{
public:

    TesselateMeshIntersect();
    void init(dolfin::Mesh *mesh_input, dolfin::Mesh *collidingMesh_input);

    const dolfin::Mesh * mesh;
    dolfin::Mesh edit_mesh;
    const dolfin::Mesh * collidingMesh;

    dolfin::BoundingVolumeTree * tree;
    dolfin::BoundingVolumeTree * collidingTree;

    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirstMesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetSecondMesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirst_0_Mesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirst_1_Mesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirst_2_Mesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetSecond_0_Mesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetSecond_1_Mesh;
    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetSecond_2_Mesh;

    boost::unordered_map<size_t, std::vector< std::pair<size_t, size_t> > > border_triangle_facet_set_first;
    boost::unordered_map<size_t, std::vector< size_t > > border_triangle_point_outside_set_first;
    boost::unordered_map<size_t, std::vector< dolfin::Mesh > > border_triangles;

     // MIDPOINT DRAW LIST
    std::vector< std::pair<dolfin::Point, std::vector< size_t > > > midpoint_triangle_indexes;

    std::shared_ptr<cutfem::CutMesh> overlapped_triangles_triangulation;
    std::shared_ptr<cutfem::CutMesh> facets_triangulation;
    std::vector<double> interface_normals;

    void clearsets();

    // TRIANGULATE
    void markborder();
    void markborderfacets();
    void triangulate_first();

    // LINE LINE INTERSECTION
    dolfin::Point intersection_two_lines(dolfin::Point p1, dolfin::Point p2, dolfin::Point p3, dolfin::Point p4);

    // CUTFEM CONFORMER
    void cutfem_conformer();

    std::shared_ptr< dolfin::CellFunction<size_t> > _domain_marker;

    void mark_ghost_penalty_facets();
    std::vector< boost::shared_ptr<dolfin::FacetFunction<size_t> > > _interior_facet_marker;


    /*
    // Refine
    void refinemesh();
    void refinemeshborder();
    void refinecollidingmeshborder();
    */


};
