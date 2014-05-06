#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <boost/unordered/unordered_map.hpp>
#include <queue>
#include <set>

class AdvancingFront
{
public:

    void get_ordered_neighbours(dolfin::Mesh * mesh_n, int index, int neighbors[3]);
    void getXYcoords(dolfin::Mesh * mesh_n, size_t index, double triangle[6]);
    int borderPointsOfXinY2(double * X, double * Y, double * P);
    int EdgeIntersections2(double * red, double * blue, int mark[3], double * points, int & nPoints);
    void computeIntersectionBetweenRedAndBlue(dolfin::Mesh * red_mesh, dolfin::Mesh * blue_mesh, int red,
            int blue, double * P, int & nP, int mark[3]);
    int swap2(double * p, double * q);
    int SortAndRemoveDoubles2(double * P, int & nP);
    double dist2(double * d1, double * d2);

    void compute_advancing(dolfin::Mesh * blue_mesh, dolfin::Mesh * red_mesh,size_t seed_red, size_t seed_blue);

    boost::unordered_map<size_t, std::vector< size_t > > intersectionsetFirstMesh;
    std::vector<std::vector<dolfin::Point> > intersection_points;
    boost::unordered_map<size_t, std::vector<std::vector<dolfin::Point> > > indexed_points;
};
