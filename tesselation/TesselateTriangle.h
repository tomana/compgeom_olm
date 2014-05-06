#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

class triangulate
{
public:
    struct TriFacet
    {
        dolfin::Point fv0;
        dolfin::Point fv1;
        dolfin::Point normal;
    };
    //std::vector<dolfin::Point> compute_triangulation_all_facets(std::vector< dolfin::Point> triangle, std::vector<TriFacet> facets);
    // helper for illustration
    static std::pair<std::vector<dolfin::Point>,  std::vector<std::vector<size_t > > > compute_triangulation_all_facets(std::vector< dolfin::Point> triangle, std::vector<TriFacet> facets);
    static std::pair<std::vector<dolfin::Point>,  std::vector<std::vector<size_t > > >  compute_triangulation_single_facet(std::vector<dolfin::Point> triangle, TriFacet facet, std::vector<std::vector<size_t > >);
    static std::pair<std::vector<dolfin::Point>,  std::vector<size_t > > compute_facet_triangulation(dolfin::Point v0, dolfin::Point v1, dolfin::Point v2, TriFacet facet);
    static double level_set_facet(dolfin::Point x, dolfin::Point p, dolfin::Point n);
    static dolfin::Point linear_interp(dolfin::Point p_0, dolfin::Point p_1, TriFacet facet);
};
