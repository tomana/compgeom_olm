// Copyright (C) 2009 Andre Massing
//
// This file is part of DOLFIN.
//
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
// First added:  2009-09-01
// Last changed: 2011-11-15

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/NoDeleter.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include "BoundingVolumeTree.h"
#include "BoundingVolumeTreeImplementation.h"


#define HAS_CGAL 1

/*
#ifdef HAS_CGAL
#include "MeshPrimitive.h"
#include "PrimitiveTraits.h"
#endif
*/

using namespace dolfin;

BoundingVolumeTree::BoundingVolumeTree(const std::string& kernel_type, const std::string& bv_type,
                                       const bool& usedop)
    : _mesh(new Mesh()),
      _label(0),
      _use_labels(false),
      _kernel_type(kernel_type),
      _usedop(usedop)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
BoundingVolumeTree::BoundingVolumeTree(const Mesh& mesh,
                                       const std::string& kernel_type,
                                       const std::string& bv_type,
                                       const bool& usedop)
    : _mesh(reference_to_no_delete_pointer(mesh)),
      _label(0),
      _use_labels(false),
      _kernel_type(kernel_type),
      _bv_type(bv_type),
      _usedop(usedop)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
BoundingVolumeTree::BoundingVolumeTree(std::shared_ptr<const Mesh> mesh,
                                       const std::string& kernel_type,
                                       const std::string& bv_type,
                                       const bool& usedop)
    : _mesh(mesh),
      _labels(new MeshFunction<std::size_t>()),
      _label(0),
      _use_labels(false),
      _kernel_type(kernel_type),
      _bv_type(bv_type),
      _usedop(usedop)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
BoundingVolumeTree::BoundingVolumeTree(const MeshFunction<std::size_t>& labels,
                                       std::size_t label,
                                       const std::string& kernel_type,
                                       const std::string& bv_type,
                                       const bool& usedop)
    : _mesh(new Mesh()),
      _labels(reference_to_no_delete_pointer(labels)),
      _label(label),
      _use_labels(true),
      _kernel_type(kernel_type),
      _bv_type(bv_type),
      _usedop(usedop)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
BoundingVolumeTree::BoundingVolumeTree(std::shared_ptr<const MeshFunction<std::size_t> > labels,
                                       std::size_t label,
                                       const std::string& kernel_type,
                                       const std::string& bv_type,
                                       const bool& usedop)
    : _mesh(new Mesh()),
      _labels(labels),
      _label(label),
      _use_labels(true),
      _kernel_type(kernel_type),
      _bv_type(bv_type),
      _usedop(usedop)
{
    // Do nothing
}
//-----------------------------------------------------------------------------
BoundingVolumeTree::~BoundingVolumeTree()
{
    // Do nothing
}
//-----------------------------------------------------------------------------
void BoundingVolumeTree::all_intersected_entities(const Point& point,
         boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    rImpl().all_intersected_entities(point, ids_result);
}

//-----------------------------------------------------------------------------
void BoundingVolumeTree::all_intersected_entities(const std::vector<Point>& points,
         boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    rImpl().all_intersected_entities(points, ids_result);
}

//-----------------------------------------------------------------------------
void BoundingVolumeTree::all_intersected_entities(const MeshEntity & entity,
         boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    rImpl().all_intersected_entities(entity,ids_result);
}
//-----------------------------------------------------------------------------
void BoundingVolumeTree::all_intersected_entities(const std::vector<MeshEntity>& entities,
         boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    rImpl().all_intersected_entities(entities, ids_result);
}

//-----------------------------------------------------------------------------
void BoundingVolumeTree::all_intersected_entities(const Mesh& another_mesh,
         boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    rImpl().all_intersected_entities(another_mesh, ids_result);
}

const BoundingVolumeTreeImplementation& BoundingVolumeTree::rImpl() const
{
    if (!_pImpl)
        _pImpl.reset(const_cast<BoundingVolumeTree *>(this)->create_intersection_operator(_kernel_type));
    return *_pImpl;
}

BoundingVolumeTreeImplementation* BoundingVolumeTree::create_intersection_operator(
    const std::string& kernel_type = "SimpleCartesian")
{
    return new BoundingVolumeTreeImplementation_d(_mesh, _bv_type, _kernel_type, _usedop);
}

void BoundingVolumeTree::clear()
{
    _pImpl.reset();
}
//-----------------------------------------------------------------------------
const Mesh& BoundingVolumeTree::mesh() const
{
    dolfin_assert(_mesh);
    return *_mesh;
}
//-----------------------------------------------------------------------------
/*
int BoundingVolumeTree::any_intersected_entity(const Point& point) const
{
    return rImpl().any_intersected_entity(point);
}
//-----------------------------------------------------------------------------
Point BoundingVolumeTree::closest_point(const Point& point) const
{
    return rImpl().closest_point(point);
}
//-----------------------------------------------------------------------------
std::size_t BoundingVolumeTree::closest_cell(const Point& point) const
{
    return rImpl().closest_cell(point);
}
//-----------------------------------------------------------------------------
std::pair<Point, std::size_t>
BoundingVolumeTree::closest_point_and_cell(const Point& point) const
{
    return rImpl().closest_point_and_cell(point);
}
//-----------------------------------------------------------------------------
double BoundingVolumeTree::distance(const Point & point) const
{
    return rImpl().distance(point);
}
//-----------------------------------------------------------------------------
// IMPLEMENT LATER
/*
void BoundingVolumeTree::reset_kernel(const std::string& kernel_type)
{
    _pImpl.reset(create_intersection_operator(kernel_type));
}
//-----------------------------------------------------------------------------
void BoundingVolumeTree::clear()
{
    _pImpl.reset();
}
//-----------------------------------------------------------------------------
const Mesh& BoundingVolumeTree::mesh() const
{
    dolfin_assert(_mesh);
    return *_mesh;
}
//-----------------------------------------------------------------------------
const BoundingVolumeTreeImplementation& BoundingVolumeTree::rImpl() const
{
    if (!_pImpl)
        _pImpl.reset(const_cast<BoundingVolumeTree *>(this)->create_intersection_operator(_kernel_type));
    return *_pImpl;
}
//-----------------------------------------------------------------------------

#ifdef HAS_CGAL
BoundingVolumeTreeImplementation*
BoundingVolumeTree::create_intersection_operator(
    const std::string& kernel_type = "SimpleCartesian")
{
    if (!_use_labels)
    {
        if (kernel_type == "ExactPredicates")
        {
            switch(_mesh->type().cell_type())
            {
            case CellType::point      :
                return new BoundingVolumeTreeImplementation_d<PointCell, EPICK>(_mesh);
            case CellType::interval   :
                return new BoundingVolumeTreeImplementation_d<IntervalCell, EPICK>(_mesh);
            case CellType::triangle   :
                return new BoundingVolumeTreeImplementation_d<TriangleCell, EPICK>(_mesh);
            case CellType::tetrahedron:
                return new BoundingVolumeTreeImplementation_d<TetrahedronCell, EPICK>(_mesh);
            default:
                dolfin_error("BoundingVolumeTree.cpp",
                             "create intersection operator",
                             "Cell type of mesh is not known. Allowed cell types are point, interval, triangle and tetrahedron");
                return 0;
            }
        }
        else  // Default is SimpleCartesion
        {
            switch( _mesh->type().cell_type())
            {
            case CellType::point      :
                return new BoundingVolumeTreeImplementation_d< PointCell, SCK  >(_mesh);
            case CellType::interval   :
                return new BoundingVolumeTreeImplementation_d< IntervalCell, SCK  >(_mesh);
            case CellType::triangle   :
                return new BoundingVolumeTreeImplementation_d< TriangleCell, SCK >(_mesh);
            case CellType::tetrahedron:
                return new BoundingVolumeTreeImplementation_d< TetrahedronCell, SCK  >(_mesh);
            default:
            {
                dolfin_error("BoundingVolumeTree.cpp",
                             "create intersection operator",
                             "Cell type of mesh is not known. Allowed cell types are point, interval, triangle and tetrahedron");
                return 0;
            }
            }
        }
    }
    else
    {
        if (kernel_type == "ExactPredicates")
        {
            switch(_labels->dim())
            {
            case CellType::point      :
                return new BoundingVolumeTreeImplementation_d<PointCell, EPICK>(*_labels, _label);
            case CellType::interval   :
                return new BoundingVolumeTreeImplementation_d<IntervalCell, EPICK>(*_labels, _label);
            case CellType::triangle   :
                return new BoundingVolumeTreeImplementation_d<TriangleCell, EPICK>(*_labels, _label);
            case CellType::tetrahedron:
                return new BoundingVolumeTreeImplementation_d<TetrahedronCell, EPICK>(*_labels, _label);
            default:
                dolfin_error("BoundingVolumeTree.cpp",
                             "create intersection operator",
                             "Cell type of mesh is not known. Allowed cell types are point, interval, triangle and tetrahedron");
                return 0;
            }
        }
        else  // Default is SimpleCartesion
        {
            switch(_labels->dim())
            {
            case CellType::point      :
                return new BoundingVolumeTreeImplementation_d< PointCell, SCK  >(*_labels, _label);
            case CellType::interval   :
                return new BoundingVolumeTreeImplementation_d< IntervalCell, SCK  >(*_labels, _label);
            case CellType::triangle   :
                return new BoundingVolumeTreeImplementation_d< TriangleCell, SCK >(*_labels, _label);
            case CellType::tetrahedron:
                return new BoundingVolumeTreeImplementation_d< TetrahedronCell, SCK  >(*_labels, _label);
            default:
                dolfin_error("BoundingVolumeTree.cpp",
                             "create intersection operator",
                             "Cell type of mesh is not known. Allowed cell types are point, interval, triangle and tetrahedron");
                return 0;
            }
        }
    }
}


#else
//If CGAL support is not available, throw an exception.
BoundingVolumeTreeImplementation*
BoundingVolumeTree::create_intersection_operator(
    const std::string & kernel_type = "SimpleCartesian")
{
    dolfin_error("BoundingVolumeTree.cpp",
                 "create intersection operator",
                 "BoundingVolumeTreeImplementation is not available, DOLFIN has been compiled without CGAL");

    // Use to member variable to avoid fuss compiler warnings
    dolfin_assert(_label);
    dolfin_assert(_use_labels);

    return 0;
}

#endif
*/
//-----------------------------------------------------------------------------
