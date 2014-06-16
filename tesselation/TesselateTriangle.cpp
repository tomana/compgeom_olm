#include "TesselateTriangle.h"

/*
std::vector<dolfin::Point> triangulate::compute_triangulation_all_facets(std::vector<dolfin::Point> triangle, std::vector<TriFacet> facets)
{

    std::vector< std::vector<size_t> > triangle_index;
    std::vector<size_t> test;
    triangle_index.push_back(test);
    std::pair<std::vector<dolfin::Point>, std::vector< std::vector<size_t > > > triangle_pair = std::make_pair(triangle, triangle_index);
    for (int i = 0; i < facets.size(); i++)
    {
        triangle_pair = compute_triangulation_single_facet(triangle_pair.first, facets.at(i), triangle_pair.second);
    }

    std::vector<dolfin::Point> triangles_marked_1;
    for (int i = 0; i < triangle_pair.second.size(); i++)
    {

        bool addtriangle = false;
        //std::cout << "Indexes" << std::endl;
        for(int j = 0; j < triangle_pair.second.at(i).size(); j ++){
        //std::cout << triangle_pair.second.at(i).at(j) << std::endl;
            if(triangle_pair.second.at(i).at(j) == 1) {
                addtriangle = true;
            }
        }
        //if(triangle_pair.second.at(i) == 1) {
        if(addtriangle) {
        triangles_marked_1.push_back(triangle_pair.first.at(i*3));
        triangles_marked_1.push_back(triangle_pair.first.at(i*3 + 1));
        triangles_marked_1.push_back(triangle_pair.first.at(i*3 + 2));
        }
        //}
    }

    return triangles_marked_1;
}
*/


std::pair <std::vector<dolfin::Point>, std::vector< std::vector<size_t > > > TesselateTriangle::compute_triangulation_all_facets(std::vector<dolfin::Point> triangle, std::vector<TriFacet> facets)
{

    std::vector< std::vector<size_t> > triangle_index;
    std::vector<size_t> test;
    triangle_index.push_back(test);
    std::pair<std::vector<dolfin::Point>, std::vector< std::vector<size_t > > > triangle_pair = std::make_pair(triangle, triangle_index);
    for (int i = 0; i < facets.size(); i++)
    {
        triangle_pair = compute_triangulation_single_facet(triangle_pair.first, facets.at(i), triangle_pair.second);
    }

    std::vector<dolfin::Point> triangles_marked_1;
    for (int i = 0; i < triangle_pair.second.size(); i++)
    {

        bool addtriangle = false;
        //std::cout << "Indexes" << std::endl;
        for(int j = 0; j < triangle_pair.second.at(i).size(); j ++){
        //std::cout << triangle_pair.second.at(i).at(j) << std::endl;
            if(triangle_pair.second.at(i).at(j) == 1) {
                addtriangle = true;
            }
        }
        //if(triangle_pair.second.at(i) == 1) {
        if(addtriangle) {
        triangles_marked_1.push_back(triangle_pair.first.at(i*3));
        triangles_marked_1.push_back(triangle_pair.first.at(i*3 + 1));
        triangles_marked_1.push_back(triangle_pair.first.at(i*3 + 2));
        }
        //}
    }
    return triangle_pair;
}


std::pair<std::vector<dolfin::Point>,  std::vector< std::vector<size_t > > > TesselateTriangle::compute_triangulation_single_facet(std::vector<dolfin::Point> triangle, TriFacet facet, std::vector< std::vector <size_t > > triangle_index)
{
    std::vector<dolfin::Point> triangle_new;
    std::vector< std::vector<size_t > > triangle_new_index;

    for (int i = 0; i < triangle_index.size(); i++)
    {
        std::pair<std::vector<dolfin::Point>,  std::vector<size_t > > triangle_temp = compute_facet_triangulation(triangle.at(i*3), triangle.at(i*3+1), triangle.at(i*3+2), facet);
        for(int j = 0; j < triangle_temp.second.size(); j++)
        {
            triangle_new.push_back(triangle_temp.first.at(j*3));
            triangle_new.push_back(triangle_temp.first.at(j*3 + 1));
            triangle_new.push_back(triangle_temp.first.at(j*3 + 2));

            std::vector<size_t > new_inner_index = triangle_index.at(i);
            new_inner_index.push_back(triangle_temp.second.at(j));
            triangle_new_index.push_back(new_inner_index);

             /*
            if(triangle_temp.second.at(j) == 1)
            {

            } else {
                triangle_new.push_back(triangle_temp.first.at(j*3));
                triangle_new.push_back(triangle_temp.first.at(j*3 + 1));
                triangle_new.push_back(triangle_temp.first.at(j*3 + 2));
                triangle_new_index.push_back(0);
            }
            */
        }
        //triangle_new.insert(triangle_new.end(), triangle_temp.begin(), triangle_temp.end());
    }
    return std::make_pair(triangle_new, triangle_new_index);
}

std::pair<std::vector<dolfin::Point>,  std::vector<size_t > > TesselateTriangle::compute_facet_triangulation(dolfin::Point v0, dolfin::Point v1, dolfin::Point v2, TriFacet facet)
{
    std::vector<dolfin::Point > triangles;
    std::vector<size_t > triangles_index;

    size_t triindex = 0;
    if (level_set_facet(v0, facet.fv0, facet.normal) < 0) triindex |= 1;
    if (level_set_facet(v1, facet.fv0, facet.normal) < 0) triindex |= 2;
    if (level_set_facet(v2, facet.fv0, facet.normal) < 0) triindex |= 4;
    switch(triindex)
    {
    case 0:
    {
        triangles.push_back(v0);
        triangles.push_back(v1);
        triangles.push_back(v2);
        triangles_index.push_back(1);
    }
    break;
    case 7:
    {
        triangles.push_back(v0);
        triangles.push_back(v1);
        triangles.push_back(v2);
        triangles_index.push_back(0);
    }
    break;
    case 1:
    {
        dolfin::Point v3 = linear_interp(v0, v1, facet);
        dolfin::Point v4 = linear_interp(v0, v2, facet);
        // +
        triangles.push_back(v0);
        triangles.push_back(v3);
        triangles.push_back(v4);
        triangles_index.push_back(0);

        // -
        triangles.push_back(v3);
        triangles.push_back(v2);
        triangles.push_back(v4);
        triangles_index.push_back(1);

        // -
        triangles.push_back(v3);
        triangles.push_back(v1);
        triangles.push_back(v2);
        triangles_index.push_back(1);
    }
    break;
    case 6:
    {
        dolfin::Point v3 = linear_interp(v0, v1, facet);
        dolfin::Point v4 = linear_interp(v0, v2, facet);
        // -
        triangles.push_back(v0);
        triangles.push_back(v3);
        triangles.push_back(v4);
        triangles_index.push_back(1);

        // +
        triangles.push_back(v3);
        triangles.push_back(v2);
        triangles.push_back(v4);
        triangles_index.push_back(0);

        // +
        triangles.push_back(v3); triangles.push_back(v1); triangles.push_back(v2);
        triangles_index.push_back(0);
    }
    break;

    case 2:
    {
        dolfin::Point v3 = linear_interp(v0, v1, facet);
        dolfin::Point v4 = linear_interp(v1, v2, facet);
        // +
        triangles.push_back(v0);
        triangles.push_back(v3);
        triangles.push_back(v2);
        triangles_index.push_back(1);

        // -
        triangles.push_back(v3);
        triangles.push_back(v4);
        triangles.push_back(v2);
        triangles_index.push_back(1);

        // -
        triangles.push_back(v3);
        triangles.push_back(v1);
        triangles.push_back(v4);
        triangles_index.push_back(0);
    }
    break;

    case 5:
    {
        dolfin::Point v3 = linear_interp(v0, v1, facet);
        dolfin::Point v4 = linear_interp(v1, v2, facet);
        // -
        triangles.push_back(v0);
        triangles.push_back(v3);
        triangles.push_back(v2);
        triangles_index.push_back(0);

        // +
        triangles.push_back(v3);
        triangles.push_back(v4);
        triangles.push_back(v2);
        triangles_index.push_back(0);

        // +
        triangles.push_back(v3);
        triangles.push_back(v1);
        triangles.push_back(v4);
        triangles_index.push_back(1);
    } break;

    case 3:
    case 4:
    {
        dolfin::Point v3 = linear_interp(v1, v2, facet);
        dolfin::Point v4 = linear_interp(v0, v2, facet);
        triangles.push_back(v0);
        triangles.push_back(v1);
        triangles.push_back(v4);

        triangles.push_back(v1);
        triangles.push_back(v3);
        triangles.push_back(v4);

        triangles.push_back(v4);
        triangles.push_back(v3);
        triangles.push_back(v2);

        if(triindex == 3) {
        triangles_index.push_back(0); triangles_index.push_back(0); triangles_index.push_back(1);
        } else {
        triangles_index.push_back(1);
        triangles_index.push_back(1);
        triangles_index.push_back(0);
        }
    } break;
    default:
    {
        std::cout << "warning: cutting case not considered" << std::endl;
    }
    break;
    }
    /*
    triangles.push_back(v0);
        triangles.push_back(v1);
        triangles.push_back(v2);
        */
    return std::make_pair(triangles, triangles_index);

}

double TesselateTriangle::level_set_facet(dolfin::Point x, dolfin::Point p, dolfin::Point n)
{

    return (x - p).dot(n);

    //return test.x()*n.x() + test.y()*n.y();
    /*
    std::cout << "dot product" << std::endl;
    std::cout << x << std::endl;
     std::cout << p << std::endl;
      std::cout << n << std::endl;
    std::cout << test << std::endl;

    return test;
    */
}

dolfin::Point TesselateTriangle::linear_interp(dolfin::Point p_0, dolfin::Point p_1, TriFacet facet)
{
    double phi_0 = level_set_facet(p_0, facet.fv0, facet.normal);
    double phi_1 = level_set_facet(p_1, facet.fv0, facet.normal);

    return p_0 + (p_0 - p_1)*(phi_0)/(phi_1-phi_0);
}
