#include "TesselateMeshIntersect.h"

TesselateMeshIntersect::TesselateMeshIntersect() :
    overlapped_triangles_triangulation(new cutfem::CutMesh(1))
    , facets_triangulation(new cutfem::CutMesh(2)),
    _interior_facet_marker(2)
{
// Do nothing
}

void TesselateMeshIntersect::init(dolfin::Mesh *mesh_input, dolfin::Mesh *collidingMesh_input)
{
    mesh = mesh_input;
    collidingMesh = collidingMesh_input;
    dolfin::CellFunction<size_t> test(*mesh, 1);
    _domain_marker = boost::shared_ptr< dolfin::CellFunction<size_t> >(new dolfin::CellFunction<size_t>(*mesh,1));

    _interior_facet_marker[0] = boost::shared_ptr<dolfin::FacetFunction<size_t> >(new dolfin::FacetFunction<size_t>(*mesh,2));
    _interior_facet_marker[1] = boost::shared_ptr<dolfin::FacetFunction<size_t> >(new dolfin::FacetFunction<size_t>(*mesh,2));
}

void TesselateMeshIntersect::triangulate_first()
{
    collidingMesh->init();
    mesh->init();
    clearsets();
    markborder();
    markborderfacets();
    border_triangles.clear();
    midpoint_triangle_indexes.clear();

    for (boost::unordered_map<size_t, std::vector< size_t > >::iterator i = border_triangle_point_outside_set_first.begin(); i != border_triangle_point_outside_set_first.end(); i++)
    {
        size_t celln = (*i).first;

        dolfin::Cell pickcell(*mesh, celln);
        dolfin::VertexIterator v(pickcell);


        std::vector<triangulate::TriFacet > facets;
        for (int j = 0; j < border_triangle_facet_set_first[(*i).first].size(); j++)
        {
            std::pair<size_t, size_t> cpair = border_triangle_facet_set_first[(*i).first][j];

            dolfin::Cell cell(*collidingMesh, cpair.first);
            dolfin::Facet f(cell.mesh(), cell.entities(1)[cpair.second]);

            // Get global index of vertices on the facet
            const std::size_t v0 = f.entities(0)[0];
            const std::size_t v1 = f.entities(0)[1];

            // Get mesh geometry
            const dolfin::MeshGeometry& geometry = cell.mesh().geometry();

            // Get the coordinates of the two vertices
            const double* p0 = geometry.x(v0);
            const double* p1 = geometry.x(v1);

            // Facet vertices as points.

            dolfin::Point fv0 = dolfin::Point(p0[0], p0[1]);
            dolfin::Point fv1 = dolfin::Point(p1[0], p1[1]);

            triangulate::TriFacet newfacet;

            newfacet.fv0 = fv0;
            newfacet.fv1 = fv1;

            size_t first_index = border_triangle_point_outside_set_first[(*i).first][0];

            newfacet.normal = f.normal();
            bool addbool3 = true;


            for(int k = 0; k < facets.size(); k++)
            {

                if(((facets[k].fv0.x() == fv0.x() && facets[k].fv1.x() == fv1.x()) && (facets[k].fv0.y() == fv0.y() && facets[k].fv1.y() == fv1.y()))
                        || ((facets[k].fv1.x() == fv0.x() && facets[k].fv0.x() == fv1.x()) && (facets[k].fv1.y() == fv0.y() && facets[k].fv0.y() == fv1.y()))
                  )
                {
                    addbool3 = false;
                }

            }


            if(addbool3)
            {
                facets.push_back(newfacet);
            }
        }

        size_t first_index = border_triangle_point_outside_set_first[(*i).first][0];
        size_t second_index = (first_index + 1) % 3;
        size_t third_index = (first_index + 2) % 3;



        bool addbool = true;
        dolfin::Point first_point = v[first_index].point();
        dolfin::Point second_point = v[second_index].point();
        dolfin::Point third_point = v[third_index].point();


        // HELPER
        std::pair <std::vector<dolfin::Point>, std::vector< std::vector<size_t > > > triangle;

        triangle.first.push_back(first_point);
        triangle.first.push_back(second_point);
        triangle.first.push_back(third_point);

        triangle = triangulate::compute_triangulation_all_facets(triangle.first, facets);

        for (int j = 0; j < triangle.second.size(); j ++)
        {
            dolfin::Mesh localmesh;
            dolfin::MeshEditor subtriangulation_editor;
            subtriangulation_editor.open(localmesh, dolfin::CellType::triangle, 2, 2);
            subtriangulation_editor.init_vertices(3);
            subtriangulation_editor.init_cells(1);
            subtriangulation_editor.add_vertex(0, triangle.first[j*3]);
            subtriangulation_editor.add_vertex(1, triangle.first[j*3+1]);
            subtriangulation_editor.add_vertex(2, triangle.first[j*3+2]);
            subtriangulation_editor.add_cell(0, 0, 1, 2);
            subtriangulation_editor.close();
            dolfin::Point midpoint_tri;
            bool addbool2 = true;

            bool addtriangle = false;
            for(int k = 0; k < triangle.second.at(j).size(); k ++)
            {
                if(triangle.second.at(j).at(k) == 1)
                {
                    addtriangle = true;
                }
            }
            std::pair<dolfin::Point, std::vector< size_t > > midpoint_pair = std::make_pair(midpoint_tri, triangle.second.at(j));

            if(addbool2)
            {
                midpoint_triangle_indexes.push_back(midpoint_pair);
                if(addtriangle) border_triangles[(*i).first].push_back(localmesh);
            }
        }
    }
    cutfem_conformer();
    mark_ghost_penalty_facets();
}

void TesselateMeshIntersect::cutfem_conformer()
{

    dolfin::MeshEditor subtriangulation_editor;
    subtriangulation_editor.open(*overlapped_triangles_triangulation, dolfin::CellType::triangle, 2, 2);


    size_t totalsize = 0;
    for (boost::unordered_map<size_t, std::vector< dolfin::Mesh > >::iterator i = border_triangles.begin(); i != border_triangles.end(); i++)
    {
        for (int j = 0; j < (*i).second.size(); j++)
        {
            for (dolfin::CellIterator cellit((*i).second[j]); !cellit.end(); ++cellit)
            {
                totalsize += 1;
            }
        }
    }
    subtriangulation_editor.init_vertices(totalsize*3);
    subtriangulation_editor.init_cells(totalsize);

    std::vector<size_t> bg_map;

    bg_map.resize(totalsize);

    overlapped_triangles_triangulation->parent_entity_maps() = std::vector< std::vector<size_t> >(1);
    overlapped_triangles_triangulation->parent_entity_maps()[0] = bg_map;

    size_t k = 0;

    for (boost::unordered_map<size_t, std::vector< dolfin::Mesh > >::iterator i = border_triangles.begin(); i != border_triangles.end(); i++)
    {
        //overlapped_triangles_triangulation->parent_entity_maps()[k] = std::vector< std::size_t >((*i).second.size() +);
        for (int j = 0; j < (*i).second.size(); j++)
        {
            overlapped_triangles_triangulation->parent_entity_maps()[0][k] = (*i).first;

            for (dolfin::CellIterator cellit((*i).second[j]); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell((*i).second[j], cellit->global_index());
                dolfin::VertexIterator v(pickcell);

                subtriangulation_editor.add_vertex(k*3, v[0].point());
                subtriangulation_editor.add_vertex(k*3+1, v[1].point());
                subtriangulation_editor.add_vertex(k*3+2, v[2].point());
                subtriangulation_editor.add_cell(k, k*3, k*3 + 1, k*3 + 2);

                k += 1;
            }
        }

    }
    subtriangulation_editor.close();


    // SUMMING BECAUSE DYNAMIC MESH EDITOR DOES NOT WORK FOR CELL TYPE INTERVAL
    size_t l = 0;
    for (boost::unordered_map<size_t, std::vector< std::pair<size_t, size_t> > >::iterator i = border_triangle_facet_set_first.begin(); i != border_triangle_facet_set_first.end(); i++)
    {
        //size_t cell = (*i).first;

        //dolfin::Cell cell2(*mesh, (*i).first);

        size_t cell = (*i).first;
        //cout << cell << endl;

        dolfin::Cell pickcell(*mesh, cell);
        dolfin::VertexIterator v(pickcell);
        for (int j = 0; j < (*i).second.size(); j++)
        {
            std::pair<size_t, size_t> cpair = (*i).second[j];
            dolfin::Cell cell(*collidingMesh, cpair.first);
            //dolfin::FacetIterator v(cell);
            //drawLines2D(v[0].point(),v[1].point());

            // Create facet from the mesh and local facet number
            dolfin::Facet f(cell.mesh(), cell.entities(1)[cpair.second]);

            // Get global index of vertices on the facet
            const std::size_t v0 = f.entities(0)[0];
            const std::size_t v1 = f.entities(0)[1];

            // Get mesh geometry
            const dolfin::MeshGeometry& geometry = cell.mesh().geometry();

            // Get the coordinates of the two vertices
            const double* p0 = geometry.x(v0);
            const double* p1 = geometry.x(v1);


            // BOTH POINTS INSIDE CELL, ADD
            if(cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                l += 1;
            }
            else if(!cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                l += 1;
            }
            else if(cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && !cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                l += 1;
            }
            else if(!cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && !cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                l += 1;
            }

        }
    }

    dolfin::MeshEditor interface_editor;

    std::vector<size_t> first_map;
    std::vector<size_t> second_map;

    first_map.resize(l);
    second_map.resize(l);

    facets_triangulation->parent_entity_maps() = std::vector< std::vector<size_t> >(2);

    facets_triangulation->parent_entity_maps()[0] = first_map;
    facets_triangulation->parent_entity_maps()[1] = second_map;

    interface_editor.open(*facets_triangulation, dolfin::CellType::interval, 1, 2);

    size_t m = 0;
    interface_editor.init_vertices(l*2);
    interface_editor.init_cells(l);

    for (boost::unordered_map<size_t, std::vector< std::pair<size_t, size_t> > >::iterator i = border_triangle_facet_set_first.begin(); i != border_triangle_facet_set_first.end(); i++)
    {
        //size_t cell = (*i).first;

        //dolfin::Cell cell2(*mesh, (*i).first);

        size_t cell = (*i).first;

        //cout << cell << endl;

        dolfin::Cell pickcell(*mesh, cell);
        dolfin::VertexIterator v(pickcell);
        for (int j = 0; j < (*i).second.size(); j++)
        {
            std::pair<size_t, size_t> cpair = (*i).second[j];
            dolfin::Cell cell(*collidingMesh, cpair.first);
            //dolfin::FacetIterator v(cell);
            //drawLines2D(v[0].point(),v[1].point());

            // Create facet from the mesh and local facet number
            dolfin::Facet f(cell.mesh(), cell.entities(1)[cpair.second]);

            // Get global index of vertices on the facet
            const std::size_t v0 = f.entities(0)[0];
            const std::size_t v1 = f.entities(0)[1];


            // Get mesh geometry
            const dolfin::MeshGeometry& geometry = cell.mesh().geometry();

            // Get the coordinates of the two vertices
            const double* p0 = geometry.x(v0);
            const double* p1 = geometry.x(v1);


            if(cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                interface_editor.add_vertex(m*2, dolfin::Point(p0[0], p0[1]));
                interface_editor.add_vertex(m*2+1, dolfin::Point(p1[0], p1[1]));
                interface_editor.add_cell(m, m*2, m*2 + 1);
                facets_triangulation->parent_entity_maps()[0][m] = (*i).first;
                facets_triangulation->parent_entity_maps()[1][m] = cpair.first;
                m += 1;

            }
            else if(!cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                // COMPUTE INTERSECTION p1 to p0, iterate over cell facets..

                dolfin::Point i0;
                for (int i = 0; i < 3; i++)
                {
                    dolfin::Facet trifacet(pickcell.mesh(), pickcell.entities(1)[i]);
                    if(cutfem::GeometryKernel::do_intersect(trifacet,f))
                    {

                        // THIS IS NASTY, BUT THIS IS HOW ITS DONE IN OTHER CODES IN DOLFIN
                        // Get global index of vertices on the facet
                        const std::size_t v0 = trifacet.entities(0)[0];
                        const std::size_t v1 = trifacet.entities(0)[1];
                        const dolfin::MeshGeometry& geometry = pickcell.mesh().geometry();
                        // Get the coordinates of the two vertices
                        const double* f0 = geometry.x(v0);
                        const double* f1 = geometry.x(v1);

                        dolfin::Point fv0 = dolfin::Point(f0[0], f0[1]);
                        dolfin::Point fv1 = dolfin::Point(f1[0], f1[1]);

                        dolfin::Point pv0 = dolfin::Point(p0[0], p0[1]);
                        dolfin::Point pv1 = dolfin::Point(p1[0], p1[1]);

                        i0 = intersection_two_lines(fv0, fv1, pv0, pv1);
                    }
                }

                // here, i0 is coordinate of intersection

                interface_editor.add_vertex(m*2, i0);
                interface_editor.add_vertex(m*2+1, dolfin::Point(p1[0], p1[1]));
                interface_editor.add_cell(m, m*2, m*2 + 1);
                facets_triangulation->parent_entity_maps()[0][m] = (*i).first;
                facets_triangulation->parent_entity_maps()[1][m] = cpair.first;
                m += 1;

            }
            else if(cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && !cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                // COMPUTE INTERSECTION p0 to p1, iterate over cell facets..
                dolfin::Point i1;
                for (int i = 0; i < 3; i++)
                {
                    dolfin::Facet trifacet(pickcell.mesh(), pickcell.entities(1)[i]);
                    if(cutfem::GeometryKernel::do_intersect(trifacet,f))
                    {

                        // THIS IS NASTY, BUT THIS IS HOW ITS DONE IN OTHER CODES IN DOLFIN
                        // Get global index of vertices on the facet
                        const std::size_t v0 = trifacet.entities(0)[0];
                        const std::size_t v1 = trifacet.entities(0)[1];
                        const dolfin::MeshGeometry& geometry = pickcell.mesh().geometry();
                        // Get the coordinates of the two vertices
                        const double* f0 = geometry.x(v0);
                        const double* f1 = geometry.x(v1);

                        dolfin::Point fv0 = dolfin::Point(f0[0], f0[1]);
                        dolfin::Point fv1 = dolfin::Point(f1[0], f1[1]);

                        dolfin::Point pv0 = dolfin::Point(p0[0], p0[1]);
                        dolfin::Point pv1 = dolfin::Point(p1[0], p1[1]);

                        i1 = intersection_two_lines(fv0, fv1, pv0, pv1);
                    }
                }

                // here, i1 is coordinate of intersection
                interface_editor.add_vertex(m*2, dolfin::Point(p0[0], p0[1]));
                interface_editor.add_vertex(m*2+1, dolfin::Point(i1[0], i1[1]));
                interface_editor.add_cell(m, m*2, m*2 + 1);
                facets_triangulation->parent_entity_maps()[0][m] = (*i).first;
                facets_triangulation->parent_entity_maps()[1][m] = cpair.first;
                m += 1;

            }
            else if(!cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p0[0], p0[1]))
                    && !cutfem::GeometryKernel::do_intersect(pickcell, dolfin::Point(p1[0], p1[1])))
            {
                /// HERE SOMETHING NEEDS TO CHANGE IF WE ARE GONNA DO THE WHOLE THING WITHOUT
                /// PRIMTIVE-PRIMITIVE INTERSECTION, THEN WE HAVE TO CHECK IF THE LINE INTERSECTS
                /// THE PRIMTIVIE AT ALL. IF IT DOESN'T, CAN WE SAFELY ADD IT TO T_1 or T_2?
                /// KEEPING TRACK OF WHICH FACETS ARE AT ALL INTERSECTED SHOULD NOT BE A PROBLEM.
                // COMPUTE INTERSECTION p0 to p1, iterate over cell facets..
                size_t intersection_index= 0;
                dolfin::Point i0;
                dolfin::Point i1;
                for (int i = 0; i < 3; i++)
                {
                    dolfin::Facet trifacet(pickcell.mesh(), pickcell.entities(1)[i]);
                    if(cutfem::GeometryKernel::do_intersect(trifacet,f))
                    {

                        // THIS IS NASTY, BUT THIS IS HOW ITS DONE IN OTHER CODES IN DOLFIN
                        // Get global index of vertices on the facet
                        const std::size_t v0 = trifacet.entities(0)[0];
                        const std::size_t v1 = trifacet.entities(0)[1];
                        const dolfin::MeshGeometry& geometry = pickcell.mesh().geometry();
                        // Get the coordinates of the two vertices
                        const double* f0 = geometry.x(v0);
                        const double* f1 = geometry.x(v1);

                        dolfin::Point fv0 = dolfin::Point(f0[0], f0[1]);
                        dolfin::Point fv1 = dolfin::Point(f1[0], f1[1]);

                        dolfin::Point pv0 = dolfin::Point(p0[0], p0[1]);
                        dolfin::Point pv1 = dolfin::Point(p1[0], p1[1]);
                        if(intersection_index == 0) {
                        i0 = intersection_two_lines(fv0, fv1, pv0, pv1);
                        intersection_index += 1;
                        } else if(intersection_index == 1) {
                        i1 = intersection_two_lines(fv0, fv1, pv0, pv1);
                        intersection_index += 1;
                        } else if(intersection_index == 2) {
                        std::cout << "Something wierd is going on, a line intersects a triangle on 3 sides" << std::endl;
                        }
                    }
                }

                // here, i1 is coordinate of intersection
                interface_editor.add_vertex(m*2, dolfin::Point(i0[0], i0[1]));
                interface_editor.add_vertex(m*2+1, dolfin::Point(i1[0], i1[1]));
                interface_editor.add_cell(m, m*2, m*2 + 1);
                facets_triangulation->parent_entity_maps()[0][m] = (*i).first;
                facets_triangulation->parent_entity_maps()[1][m] = cpair.first;
                m += 1;

            }
        }
    }
    //std::cout << m << std::endl;
    interface_editor.close();
}

void TesselateMeshIntersect::markborder()
{
    mesh->init();
    collidingMesh->init();
    //std::cout << "Marking border" << std::endl;
    tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
    intersectionsetFirstMesh.clear();
    tree->all_intersected_entities(*collidingMesh, intersectionsetFirstMesh);

    dolfin::CellFunction<bool> cell_markers(*mesh);
    cell_markers.set_all(true);

    dolfin::CellFunction<size_t> & domain_marker = *_domain_marker;
    domain_marker.set_all(0);

    intersectionsetSecond_1_Mesh.clear();
    intersectionsetFirst_0_Mesh.clear();
    intersectionsetFirst_1_Mesh.clear();
    intersectionsetFirst_2_Mesh.clear();
    border_triangle_point_outside_set_first.clear();

    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirstMesh.begin(); i != intersectionsetFirstMesh.end(); i++)
    {
        if (intersectionsetFirstMesh[(*i).first].size() > 0)
        {
            cell_markers[(*i).first] = false;

            // SKETCHY, BUT WORKS FOR FULLY SUBMERGED COLLIDING MESH

            size_t cell = (*i).first;
            dolfin::MeshEntity pickcell(*mesh, mesh->topology().dim(), cell);
            bool completelyoverlapped = true;
            for (int j = 0; j < pickcell.num_entities(mesh->topology().dim()); j++)
            {
                if (intersectionsetFirstMesh.find(pickcell.entities(mesh->topology().dim())[j]) == intersectionsetFirstMesh.end())
                {
                    cell_markers[cell] = false;
                    completelyoverlapped = false;
                    for(int k = 0; k < intersectionsetFirstMesh[(*i).first].size(); k++)
                    {
                        intersectionsetFirst_1_Mesh[(*i).first].push_back(intersectionsetFirstMesh[(*i).first][k]);
                        std::vector<size_t> test;
                        intersectionsetSecond_1_Mesh[intersectionsetFirstMesh[(*i).first][k]] = test;
                        domain_marker[(*i).first] = 1;
                    }
                }
            }
            if(completelyoverlapped)
            {
                std::vector<size_t> test;
                intersectionsetFirst_2_Mesh[(*i).first] = test;
                cell_markers[cell] = false;
                domain_marker[(*i).first] = 2;
            }
        }
    }

    for (dolfin::CellIterator cellit(*mesh); !cellit.end(); ++cellit)
    {
        if (cell_markers[cellit->index()] == true)
        {

            std::vector<size_t> test;
            intersectionsetFirst_0_Mesh[cellit->index()] = test;

        }
    }

    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirst_1_Mesh.begin(); i != intersectionsetFirst_1_Mesh.end(); i++)
    {
        size_t cell = (*i).first;
        //cout << cell << endl;
        dolfin::Cell pickcell(*mesh, cell);
        dolfin::VertexIterator v(pickcell);

        // HERE SOME SELECTION GOES WRONG... OR, MAYBE NOT.
        for (dolfin::VertexIterator v(pickcell); !v.end(); ++v)
        {
            dolfin::Vertex testv = (*v);
            size_t localindex = v.pos();

            for (int j = 0; j < testv.num_entities(mesh->topology().dim()); ++j)
            {
                dolfin::Cell pickcell2(*mesh, testv.entities(mesh->topology().dim())[j]);

                if (intersectionsetFirst_0_Mesh.find(pickcell2.index()) != intersectionsetFirst_0_Mesh.end())
                {
                    border_triangle_point_outside_set_first[(*i).first].push_back(localindex);
                }
            }
        }
    }

}

void TesselateMeshIntersect::markborderfacets()
{
    border_triangle_facet_set_first.clear();
    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirst_1_Mesh.begin(); i != intersectionsetFirst_1_Mesh.end(); i++)
    {
        for(int k = 0; k < intersectionsetFirst_1_Mesh[(*i).first].size(); k++)
        {
            size_t cell = intersectionsetFirst_1_Mesh[(*i).first][k];
            //cout << cell << endl;
            dolfin::MeshEntity pickcell(*collidingMesh, collidingMesh->topology().dim(), cell);
            for (dolfin::FacetIterator facet(pickcell); !facet.end(); ++facet)
            {
                if(facet->exterior())
                {
                    dolfin::MeshEntity pickcell2(*mesh, mesh->topology().dim(), (*i).first);
                    if (cutfem::GeometryKernel::do_intersect(pickcell2, *facet))
                    {
                        border_triangle_facet_set_first[(*i).first].push_back(std::make_pair(cell, pickcell.index(*facet)));
                    }
                }
            }
        }
    }

    for (boost::unordered_map<size_t, std::vector< std::pair<size_t, size_t> > >::iterator i = border_triangle_facet_set_first.begin(); i != border_triangle_facet_set_first.end(); i++)
    {

    sort( (*i).second.begin(), (*i).second.end() );
            (*i).second.erase( unique( (*i).second.begin(), (*i).second.end() ), (*i).second.end() );
    }
}

void TesselateMeshIntersect::clearsets()
{
    intersectionsetFirstMesh.clear();
    intersectionsetSecondMesh.clear();
    intersectionsetFirst_0_Mesh.clear();
    intersectionsetFirst_1_Mesh.clear();
    intersectionsetFirst_2_Mesh.clear();
    intersectionsetSecond_0_Mesh.clear();
    intersectionsetSecond_1_Mesh.clear();
    intersectionsetSecond_2_Mesh.clear();
    border_triangle_facet_set_first.clear();
    border_triangle_point_outside_set_first.clear();
    border_triangles.clear();
}

void TesselateMeshIntersect::refinemesh()
{
    mesh->init();
    collidingMesh->init();
    std::cout << "Building tree for refinement" << std::endl;
    tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
    intersectionsetFirst_1_Mesh.clear();
    tree->all_intersected_entities(*collidingMesh, intersectionsetFirst_1_Mesh);

    dolfin::CellFunction<bool> cell_markers(*mesh);
    cell_markers.set_all(false);

    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirst_1_Mesh.begin(); i != intersectionsetFirst_1_Mesh.end(); i++)
    {
        size_t cell = (*i).first;
        //cout << cell << endl;
        cell_markers[cell] = true;
    }
//    dolfin::compute_connectivity(*mesh, mesh->topology().dim(), mesh->topology().dim());
    dolfin::Mesh test = dolfin::refine(*mesh, cell_markers);
    test.init();
    mesh = new dolfin::Mesh(test);
    tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
    intersectionsetFirst_1_Mesh.clear();
    tree->all_intersected_entities(*collidingMesh, intersectionsetFirst_1_Mesh);
}

void TesselateMeshIntersect::refinemeshborder()
{
    mesh->init();
    collidingMesh->init();
    std::cout << "Building trees for refinement" << std::endl;
    tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
    intersectionsetFirstMesh.clear();
    tree->all_intersected_entities(*collidingMesh, intersectionsetFirstMesh);

    dolfin::CellFunction<bool> cell_markers(*mesh);
    cell_markers.set_all(false);

    dolfin::CellFunction<bool> cell_markers2(*collidingMesh);
    cell_markers2.set_all(false);

    intersectionsetFirst_0_Mesh.clear();
    intersectionsetFirst_1_Mesh.clear();
    intersectionsetFirst_2_Mesh.clear();
    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirstMesh.begin(); i != intersectionsetFirstMesh.end(); i++)
    {
        if (intersectionsetFirstMesh[(*i).first].size() > 0)
        {
            /*
            if(intersectionsetFirstMesh[(*i).first])
            for(int j = 0; j < intersectionsetFirstMesh[(*i).first].size(); j++) {
                if(intersectionsetSecondMesh(intersectionsetFirstMesh[(*i).first].at(j))){


                }

            }
            */
            //intersectionsetFirst_1_Mesh[(*i).first].push_back(0);
            size_t cell = (*i).first;
            dolfin::MeshEntity pickcell(*mesh, mesh->topology().dim(), cell);
            bool test = false;
            for (int j = 0; j < pickcell.num_entities(mesh->topology().dim()); j++)
            {
                if ( intersectionsetFirstMesh.find(pickcell.entities(mesh->topology().dim())[j]) == intersectionsetFirstMesh.end() )
                {
                    cell_markers[cell] = true;
                    for(int k = 0; k < intersectionsetFirstMesh[(*i).first].size(); k++)
                    {
                        cell_markers2[intersectionsetFirstMesh[(*i).first][k]] = true;
                    }
                }
                else
                {

                }
            }


        }


    }

    dolfin::Mesh test = dolfin::refine(*mesh, cell_markers);
    test.init();
    mesh = new dolfin::Mesh(test);
    markborder();
    markborderfacets();

}


void TesselateMeshIntersect::refinecollidingmeshborder()
{
    mesh->init();
    collidingMesh->init();
    std::cout << "Building trees for refinement" << std::endl;
    tree = new dolfin::BoundingVolumeTree(*mesh, "SimpleCartesian", "DOP26", true);
    intersectionsetFirstMesh.clear();
    tree->all_intersected_entities(*collidingMesh, intersectionsetFirstMesh);

    dolfin::CellFunction<bool> cell_markers(*mesh);
    cell_markers.set_all(false);

    dolfin::CellFunction<bool> cell_markers2(*collidingMesh);
    cell_markers2.set_all(false);

    intersectionsetFirst_0_Mesh.clear();
    intersectionsetFirst_1_Mesh.clear();
    intersectionsetFirst_2_Mesh.clear();
    for (boost::unordered_map<size_t, std::vector< size_t> >::iterator i = intersectionsetFirstMesh.begin(); i != intersectionsetFirstMesh.end(); i++)
    {
        if (intersectionsetFirstMesh[(*i).first].size() > 0)
        {
            size_t cell = (*i).first;
            dolfin::MeshEntity pickcell(*mesh, mesh->topology().dim(), cell);
            bool test = false;
            for (int j = 0; j < pickcell.num_entities(mesh->topology().dim()); j++)
            {
                if ( intersectionsetFirstMesh.find(pickcell.entities(mesh->topology().dim())[j]) == intersectionsetFirstMesh.end() )
                {
                    cell_markers[cell] = true;
                    for(int k = 0; k < intersectionsetFirstMesh[(*i).first].size(); k++)
                    {
                        cell_markers2[intersectionsetFirstMesh[(*i).first][k]] = true;
                    }
                }
                else
                {

                }
            }
        }
    }

    dolfin::Mesh test2 = dolfin::refine(*collidingMesh, cell_markers2);
    test2.init();
    collidingMesh = new dolfin::Mesh(test2);
    markborder();
    markborderfacets();
}

/*
dolfin::Point TesselateMeshIntersect::line_plane_intersection(dolfin::Point l_0, dolfin::Point l_1, dolfin::Point p_0, dolfin::Point normal)
{
    dolfin::Point l = l_0 -l_1 ;
    double d = ((p_0 - l_0).dot(normal))/(l.dot(normal));
    return (d*l) + l_0;
}

bool TesselateMeshIntersect::point_on_line(dolfin::Point linePointA, dolfin::Point linePointB, dolfin::Point point)
{
    const double EPSILON = 0.0000000001f;

    double a = (linePointB.y() - linePointA.y()) / (linePointB.x() - linePointA.x());

    double b = linePointA.y() - a * linePointA.x();

    double maxx = std::max(linePointA.x(), linePointB.x());
    double minx = std::min(linePointA.x(), linePointB.x());
    double maxy = std::max(linePointA.y(), linePointB.y());
    double miny = std::min(linePointA.y(), linePointB.y());

    if ( abs(point.y() - (a*point.x()+b)) < EPSILON)
    {
        if(point.y() <= maxy & point.y() >= miny & point.x() <= maxx & point.x() >= minx)
        {
            return true;
        }
    }
    return false;
}

void add_unique_point(std::vector< dolfin::Point> * vertices, dolfin::Point point)
{
    if (vertices->size() == 0)
    {
        vertices->push_back(point);
    }
    else
    {
        for (int i = 0; i < vertices->size(); i++)
        {
            if(vertices->at(i).squared_distance(point) < 0.00001)
            {
                return;
            }
        }
        vertices->push_back(point);
    }
}
*/

dolfin::Point TesselateMeshIntersect::intersection_two_lines(dolfin::Point p1, dolfin::Point p2, dolfin::Point p3, dolfin::Point p4)
{
// Store the values for fast access and easy
// equations-to-code conversion
    double x1 = p1.x(), x2 = p2.x(), x3 = p3.x(), x4 = p4.x();
    double y1 = p1.y(), y2 = p2.y(), y3 = p3.y(), y4 = p4.y();

    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
// If d is zero, there is no intersection
    if (d == 0) return NULL;

// Get the x and y
    double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
    double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;

// Check if the x and y coordinates are within both lines
    if ( x < std::min(x1, x2) || x > std::max(x1, x2) ||
            x < std::min(x3, x4) || x > std::max(x3, x4) ) return NULL;
    if ( y < std::min(y1, y2) || y > std::max(y1, y2) ||
            y < std::min(y3, y4) || y > std::max(y3, y4) ) return NULL;

// Return the point of intersection
    dolfin::Point ret = dolfin::Point(x, y);
    return ret;
}

void TesselateMeshIntersect::mark_ghost_penalty_facets()
{
  dolfin::CellFunction<size_t> domain_marker = *_domain_marker;
  size_t _gdim = 2;
  mesh->init(_gdim - 1, _gdim);

  for (dolfin::FacetIterator facet(*mesh); !facet.end(); ++facet)
  {
    if (facet->num_entities(_gdim) != 1)
    {
      const size_t marker0 =  domain_marker[ facet->entities(_gdim)[0] ];
      const size_t marker1 =  domain_marker[ facet->entities(_gdim)[1] ];
      // 1. Mark facet incident with a cut cell
      if ((marker0 == 1 && marker1 == 1))
      {
        (*_interior_facet_marker[0])[*facet] = 1;
        (*_interior_facet_marker[1])[*facet] = 1;
      }
      if((marker0 == 1 && marker1 == 0) || (marker0 == 0 && marker1 == 1))
      {
        (*_interior_facet_marker[0])[*facet] = 1;
        (*_interior_facet_marker[1])[*facet] = 0;
      }
      if((marker0 == 1 && marker1 == 2) || (marker0 == 2 && marker1 == 1))
        (*_interior_facet_marker[1])[*facet] = 1;
      // 2. Mark interior facets
      if ((marker0 == 0 && marker1 == 0))
      {
        (*_interior_facet_marker[0])[*facet] = 0;
        (*_interior_facet_marker[1])[*facet] = 0;
      }
      if ((marker0 == 2 && marker1 == 2))
      {
        (*_interior_facet_marker[0])[*facet] = 2;
        (*_interior_facet_marker[1])[*facet] = 2;
      }
    }
  }
}
