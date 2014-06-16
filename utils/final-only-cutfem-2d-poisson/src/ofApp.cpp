#include "ofApp.h"

//#include <typeinfo>
class InnerDomain : public dolfin::SubDomain
{
public:
    virtual bool inside(const dolfin::Array<double>& x, bool on_boundary) const
    {
        return std::abs(x[0]) < 1.0 + DOLFIN_EPS
               and std::abs(x[1]) < 1.0 + DOLFIN_EPS;
    }
};
//--------------------------------------------------------------
void ofApp::setup()
{
    size_t N = 1;
    // Build background mesh
    double a = 3.0;
    double b = -3.0;

    std::shared_ptr<dolfin::Mesh> mesh0(
        new dolfin::RectangleMesh(a, a, b, b, 3 * N, 3 * N));
    mesh0->init(mesh0->topology().dim() - 1, mesh0->topology().dim());

    // Create sub mesh
    std::shared_ptr<dolfin::CellFunction<std::size_t> > domain_marker_0
    (new dolfin::CellFunction<std::size_t>(*mesh0, 0));

    InnerDomain().mark(*domain_marker_0, 2);
    std::shared_ptr<dolfin::SubMesh> mesh1(
        new dolfin::SubMesh(*mesh0, *domain_marker_0, 2));

    double scalefactor0 = 1/6.0;
    dolfin::MeshGeometry& geometry0 = mesh0->geometry();
    for (std::size_t i = 0; i < geometry0.size(); i++)
    {
        // Get coordinate
        double* x = geometry0.x(i);
        // Scale
        const double x0 = x[0]*scalefactor0 + 0.5;
        const double x1 = x[1]*scalefactor0 + 0.5;
        // Store coordinate
        x[0] = x0;
        x[1] = x1;
    }

    double scalefactor = 1.0001* 1/6.0;
    dolfin::MeshGeometry& geometry = mesh1->geometry();
    for (std::size_t i = 0; i < geometry.size(); i++)
    {
        // Get coordinate
        double* x = geometry.x(i);
        // Scale
        const double x0 = x[0]*scalefactor + 0.5;
        const double x1 = x[1]*scalefactor + 0.5;
        // Store coordinate
        x[0] = x0;
        x[1] = x1;
    }
    mesh = new dolfin::Mesh(*mesh0);
    collidingMesh = new dolfin::Mesh(*mesh1);

    //collidingMesh->rotate(40,2);

    //mesh = new dolfin::Mesh("data/mesh_f.xml");
    //collidingMesh = new dolfin::Mesh("data/mesh2_f.xml");
    mesh->init();
    collidingMesh->init();

    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes_new(new cutfem::TriangulatedOverlappingMeshes);
    overlapping_meshes = overlapping_meshes_new;

    std::shared_ptr<dolfin::Mesh> mesh_mesh(new dolfin::Mesh(*mesh));
    std::shared_ptr<dolfin::Mesh> mesh_collidingMesh(new dolfin::Mesh(*collidingMesh));
    overlapping_meshes->add(mesh_mesh);
    overlapping_meshes->add(mesh_collidingMesh);

    // Computing intersections
    overlapping_meshes->compute_collisions();
    overlapping_meshes->compute_intersections();


    drawsomething = true;
    outputeps = false;
    outputtex = false;
    cam.setNearClip(0.000000000001);
    savedPose = ofMatrix4x4(0.999998, 1.06809e-09, -0.00241552,        0,
                            0.00222934, 0.384979, 0.922922,        0,
                            0.000929923, -0.922926, 0.384978,        0,
                            -612.275, -1651.12,  185.784,        1
                           );



    cam.setTransformMatrix(savedPose);

    cam.enableOrtho();
    checksets = true;

    /// SETTING MESHES TO EXAMPLE SET REMOVE IF USING MESHES ABOVE
    //mesh2dbuttonPressed();

    gui.setup("Setup mesh"); // most of the time you don't need a name
    gui.setPosition(10, 10);
    gui2.setup("Draw solution");
    gui2.setPosition(220, 10);

    mesh2dbutton.addListener(this, &ofApp::mesh2dbuttonPressed);
    mesh2dgenbutton.addListener(this, &ofApp::mesh2dgenbuttonPressed);

    gui.add(draw_indexes.setup("draw_indexes", true));
    gui.add(draw_first_mesh.setup("draw_first_mesh", true));
    gui.add(draw_first_mesh_u.setup("draw_first_mesh_u", true));
    gui.add(draw_first_mesh_t_0_u.setup("draw_first_mesh_t_0_u", true));
    gui.add(draw_first_mesh_t_0_t_1_u.setup("draw_first_mesh_t_0_t_1_u", true));
    gui.add(draw_second_mesh.setup("draw_second_mesh", true));
    gui.add(draw_second_mesh_u.setup("draw_second_mesh_u", true));
    gui.add(draw_t_0_1.setup("draw_t_0_1", true));
    gui.add(draw_t_0_gamma.setup("draw_t_0_gamma", true));
    gui.add(draw_t_0_2.setup("draw_t_0_2", true));
    gui.add(draw_border_first.setup("draw_border_first", true));
    gui.add(draw_triangulation_mesh.setup("draw_triangulation_mesh", true));
    gui.add(draw_triangulation_mesh_u.setup("draw_triangulation_mesh_u", true));
    gui.add(draw_index_triangulation.setup("draw_index_triangulation", true));

    gui.add(mesh2dbutton.setup("2d mesh"));
    gui.add(mesh2dgenbutton.setup("2d generated mesh"));

    gui.add(meshu_r.setup( "mesh Uniform", 8, 1, 100 ));
    gui.add(meshx_r.setup( "mesh res X", 8, 1, 100 ));
    gui.add(meshy_r.setup( "mesh res Y", 8, 1, 100 ));
    gui.add(meshz_r.setup( "mesh res Z", 8, 1, 100 ));

    gui.add(cmeshu_r.setup( "cmesh Uniform", 8, 1, 100 ));
    gui.add(cmeshx_r.setup( "cmesh res X", 8, 1, 100 ));
    gui.add(cmeshy_r.setup( "cmesh res Y", 8, 1, 100 ));
    gui.add(cmeshz_r.setup( "cmesh res Z", 8, 1, 100 ));
    gui.add(gamma_u.setup("gamma_u", 10, 0, 100 ));

    gui2.add(zaxisscaler.setup( "zaxisscaler", 0.5, 0, 10 ));
    gui2.add(transparency.setup( "transparency", 150, 0, 255 ));
    gui2.add(svg_n_slider.setup( "svg_n_slider", 0.5, 0, 10 ));
    gui2.add(firstcolor.set("firstcolor",ofColor(255,128,128),ofColor(0,0),ofColor(255,255)));
    gui2.add(secondcolor.set("secondcolor",ofColor(128,128,255),ofColor(0,0),ofColor(255,255)));
    gui2.add(allscaler.setup("allscaler", 555, 0, 5000));
    gui2.add(rotatez.setup("rotatez", 225, 0, 360 ));
    gui2.add(trans_u_z.setup("trans_u_z", 0, -2, 2 ));
    gui2.add(trans_first_mesh.setup("trans_first_mesh", 0, -2, 2 ));
    gui2.add(trans_second_mesh.setup("trans_second_mesh", 0, -2, 2 ));
    gui2.add(trans_t_0_1.setup("trans_t_0_1", 0, -2, 2 ));
    gui2.add(trans_t_0_gamma.setup("trans_t_0_gamma", 0, -2, 2 ));
    gui2.add(trans_t_0_2.setup("trans_t_0_2", 0, -2, 2 ));
    gui2.add(trans_triangulation_mesh.setup("trans_triangulation_mesh", 0, -2, 2 ));
    gui2.add(trans_border_first.setup("trans_border_first", 0, -2, 2 ));



    /*
    size_t ints_mesh[] = {3,2,4,203};
    size_t ints_collidingMesh[] = {33,52,53,39};
    size_t ints_subtriang_mesh[] = {3};
    size_t ints_facettriang_collidingMesh[] = {4,203};
    */

    size_t ints_mesh[] = {};
    size_t ints_collidingMesh[] = {};
    size_t ints_subtriang_mesh[] = {};
    size_t ints_facettriang_collidingMesh[] = {};
    std::vector<size_t> ints_mesh_vec (ints_mesh, ints_mesh + sizeof(ints_mesh) / sizeof(size_t) );
    std::vector<size_t> ints_collidingMesh_vec (ints_collidingMesh, ints_collidingMesh + sizeof(ints_collidingMesh) / sizeof(size_t) );
    std::vector<size_t> ints_subtriang_mesh_vec (ints_subtriang_mesh, ints_subtriang_mesh + sizeof(ints_subtriang_mesh) / sizeof(size_t) );
    std::vector<size_t> ints_facettriang_collidingMesh_vec (ints_facettriang_collidingMesh, ints_facettriang_collidingMesh + sizeof(ints_facettriang_collidingMesh) / sizeof(size_t) );
    draw_indexes_mesh_vec = ints_mesh_vec;
    draw_indexes_collidingMesh_vec = ints_collidingMesh_vec;
    draw_indexes_subtriang_mesh_vec = ints_subtriang_mesh_vec;
    draw_indexes_facettriang_collidingMesh_vec = ints_facettriang_collidingMesh_vec;
    cout << "DONE SETUP" << endl;
}

void ofApp::mesh2dbuttonPressed()
{
    checksets = true;
    mesh = new dolfin::Mesh("data/mesh.xdmf");
    collidingMesh = new dolfin::Mesh("data/mesh2.xdmf");

    mesh->init();
    collidingMesh->init();

    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes_new(new cutfem::TriangulatedOverlappingMeshes);
    overlapping_meshes = overlapping_meshes_new;

    std::shared_ptr<dolfin::Mesh> mesh_mesh(new dolfin::Mesh(*mesh));
    std::shared_ptr<dolfin::Mesh> mesh_collidingMesh(new dolfin::Mesh(*collidingMesh));
    overlapping_meshes->add(mesh_mesh);
    overlapping_meshes->add(mesh_collidingMesh);

    // Computing intersections
    overlapping_meshes->compute_collisions();
    overlapping_meshes->compute_intersections();


    savedPose = ofMatrix4x4(                           0.999998, 1.06809e-09, -0.00241552,        0,
                0.00222934, 0.384979, 0.922922,        0,
                0.000929923, -0.922926, 0.384978,        0,
                -870.933, -2049.85,  446.463,        1
                           );


    cam.setTransformMatrix(savedPose);
}

void ofApp::mesh2dgenbuttonPressed()
{
    checksets = true;

    mesh = new dolfin::RectangleMesh(0, 0, 1, 1, meshx_r, meshy_r);
    collidingMesh = new dolfin::RectangleMesh(0.3001, 0.3001, 0.6001, 0.6001,cmeshx_r, cmeshy_r);

    mesh->init();
    collidingMesh->init();

    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes_new(new cutfem::TriangulatedOverlappingMeshes);
    overlapping_meshes = overlapping_meshes_new;

    std::shared_ptr<dolfin::Mesh> mesh_mesh(new dolfin::Mesh(*mesh));
    std::shared_ptr<dolfin::Mesh> mesh_collidingMesh(new dolfin::Mesh(*collidingMesh));
    overlapping_meshes->add(mesh_mesh);
    overlapping_meshes->add(mesh_collidingMesh);

    // Computing intersections
    overlapping_meshes->compute_collisions();
    overlapping_meshes->compute_intersections();

    savedPose = ofMatrix4x4(0.999998, 1.06809e-09, -0.00241552,        0,
                            0.00222934, 0.384979, 0.922922,        0,
                            0.000929923, -0.922926, 0.384978,        0,
                            -612.275, -1651.12,  185.784,        1
                           );
    cam.setTransformMatrix(savedPose);

}

void ofApp::meshuniform()
{
    int x = meshu_r;
    meshx_r = x;
    int y = meshu_r;
    meshy_r = y;
    int z = meshu_r;
    meshz_r = z;
}

void ofApp::cmeshuniform()
{
    int x = cmeshu_r;
    cmeshx_r = x;
    int y = cmeshu_r;
    cmeshy_r = y;
    int z = cmeshu_r;
    cmeshz_r = z;
}

//--------------------------------------------------------------
void ofApp::update()
{

    if (meshu_r != lastmeshu_r)
    {
        meshuniform();
    }
    if (cmeshu_r != lastcmeshu_r)
    {
        cmeshuniform();
    }
    lastmeshu_r = meshu_r;
    lastcmeshu_r = cmeshu_r;

    /*
    if(checksets == true)
    {

        if (triangulation.intersectionsetFirst_1_Mesh.size() > 0)
        {
            draw_first_mesh = false;
            draw_second_mesh = true;
            draw_first_0 = true;
            draw_first_1 = true;
            draw_first_2 = true;
            draw_second_1 = true;
            draw_second_0 = false;
            draw_second_2 = false;
        }
        else
        {
            draw_first_mesh = true;
            draw_second_mesh = true;
            draw_first_1 = true;
            draw_second_1 = true;
            draw_first_0 = false;
            draw_first_2 = false;
            draw_second_0 = false;
            draw_second_2 = false;
        }
        checksets = false;
    }
    */
}

//--------------------------------------------------------------
void ofApp::draw()
{

    ofBackground(255,255,255,255);

    if(drawsomething)
    {
        ofSetColor(255);
        float scalefactor = allscaler;
        glEnable(GL_BLEND);
        cam.begin();
        ofRotate(rotatez, 0,0,1);
        ofScale(scalefactor,scalefactor,scalefactor);
        glEnable(GL_DEPTH_TEST);

        // DRAW BACKGROUND MESH
        if(draw_first_mesh)
        {

            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            ofPushMatrix();
            ofTranslate(0,0,trans_first_mesh);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {
                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    float x = pickcell.midpoint().x();
                    float y = pickcell.midpoint().y();
                    if(draw_indexes) draw_or_save_string(x, y, indexstring);

                    drawTetLines2D(v[0].point(),v[1].point(),v[2].point());

                }
            }
            ofPopMatrix();
        }

        if(true)
        {

            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {
                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    float x = pickcell.midpoint().x();
                    float y = pickcell.midpoint().y();
                    if(draw_t_0_1)
                    {
                        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 0)
                        {
                            ofPushMatrix();
                            ofTranslate(0,0,trans_t_0_1);
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            drawTetLines2D(v[0].point(),v[1].point(),v[2].point());
                            ofPopMatrix();
                        }
                    }
                    if(draw_t_0_gamma)
                    {
                        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 1)
                        {
                            ofPushMatrix();
                            ofTranslate(0,0,trans_t_0_gamma);
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            drawTetLines2D(v[0].point(),v[1].point(),v[2].point());
                            ofPopMatrix();
                        }
                    }
                    if(draw_t_0_2)
                    {
                        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 2)
                        {
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            ofPushMatrix();
                            ofTranslate(0,0,trans_t_0_2);
                            drawTetLines2D(v[0].point(),v[1].point(),v[2].point());
                            ofPopMatrix();
                        }
                    }


                }
            }
        }


        // DRAW SOLUTION ON BACKGROUND MESH
        if(draw_first_mesh_u)
        {
            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            gl2psEnable(GL2PS_BLEND);
            ofPushMatrix();
            ofTranslate(0,0,trans_u_z);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {


                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {

                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    dolfin::Point v0_u;
                    dolfin::Point v1_u;
                    dolfin::Point v2_u;
                    if(vert_values_1.size() > 0)
                    {
                        //v[2].entities(0)[0]][0]
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), zaxisscaler*vert_values_1[v[0].index()]);
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), zaxisscaler*vert_values_1[v[1].index()]);
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), zaxisscaler*vert_values_1[v[2].index()]);
                        ofColor r = firstcolor;
                        ofColor b = secondcolor;
                        ofColor c0;
                        ofColor c1;
                        ofColor c2;
                        if(u_min_1 == 0 && u_max_2 == 0)
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], 0, 0.1, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], 0, 0.1, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], 0, 0.1, 0, 1));
                        }
                        else
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], u_min_1, u_max_2, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], u_min_1, u_max_2, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], u_min_1, u_max_2, 0, 1));
                        }

                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        //if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 0) //| overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 1)
                        //{
                        if(draw_indexes) draw_or_save_string(x, y, indexstring);
                        int alpha = transparency;
                        c0.a = alpha;
                        c1.a = alpha;
                        c2.a = alpha;
                        drawTriangleColor2D(v0_u,v1_u,v2_u,c0, c1, c2);

                        alpha = 255;
                        c0.a = alpha;
                        c1.a = alpha;
                        c2.a = alpha;
                        ofColor bl = ofColor::black;
                        drawTriangleLinesColor2D(v0_u,v1_u,v2_u, bl, bl, bl);
                        //}

                    }
                    else
                    {
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), v[0].point().z());
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), v[1].point().z());
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), v[2].point().z());

                        drawTetLines2D(v0_u,v1_u,v2_u);
                    }



                }
            }
            ofPopMatrix();
        }

        // DRAW SOLUTION ON BACKGROUND MESH
        if(draw_first_mesh_t_0_t_1_u)
        {
            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            gl2psEnable(GL2PS_BLEND);
            ofPushMatrix();
            ofTranslate(0,0,trans_u_z);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {


                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {

                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    dolfin::Point v0_u;
                    dolfin::Point v1_u;
                    dolfin::Point v2_u;
                    if(vert_values_1.size() > 0)
                    {
                        //v[2].entities(0)[0]][0]
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), zaxisscaler*vert_values_1[v[0].index()]);
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), zaxisscaler*vert_values_1[v[1].index()]);
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), zaxisscaler*vert_values_1[v[2].index()]);
                        ofColor r = firstcolor;
                        ofColor b = secondcolor;
                        ofColor c0;
                        ofColor c1;
                        ofColor c2;
                        if(u_min_1 == 0 && u_max_2 == 0)
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], 0, 0.1, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], 0, 0.1, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], 0, 0.1, 0, 1));
                        }
                        else
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], u_min_1, u_max_2, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], u_min_1, u_max_2, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], u_min_1, u_max_2, 0, 1));
                        }

                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 0 || overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 1)
                        {
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            int alpha = transparency;
                            c0.a = alpha;
                            c1.a = alpha;
                            c2.a = alpha;
                            drawTriangleColor2D(v0_u,v1_u,v2_u,c0, c1, c2);

                            alpha = 255;
                            c0.a = alpha;
                            c1.a = alpha;
                            c2.a = alpha;
                            ofColor bl = ofColor::black;
                            drawTriangleLinesColor2D(v0_u,v1_u,v2_u, bl, bl, bl);
                        }

                    }
                    else
                    {
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), v[0].point().z());
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), v[1].point().z());
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), v[2].point().z());

                        drawTetLines2D(v0_u,v1_u,v2_u);
                    }



                }
            }
            ofPopMatrix();
        }

        // DRAW SOLUTION ON BACKGROUND MESH
        if(draw_first_mesh_t_0_u)
        {
            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            gl2psEnable(GL2PS_BLEND);
            ofPushMatrix();
            ofTranslate(0,0,trans_u_z);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {


                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {

                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    dolfin::Point v0_u;
                    dolfin::Point v1_u;
                    dolfin::Point v2_u;
                    if(vert_values_1.size() > 0)
                    {
                        //v[2].entities(0)[0]][0]
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), zaxisscaler*vert_values_1[v[0].index()]);
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), zaxisscaler*vert_values_1[v[1].index()]);
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), zaxisscaler*vert_values_1[v[2].index()]);
                        ofColor r = firstcolor;
                        ofColor b = secondcolor;
                        ofColor c0;
                        ofColor c1;
                        ofColor c2;
                        if(u_min_1 == 0 && u_max_2 == 0)
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], 0, 0.1, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], 0, 0.1, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], 0, 0.1, 0, 1));
                        }
                        else
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_1[v[0].index()], u_min_1, u_max_2, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_1[v[1].index()], u_min_1, u_max_2, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_1[v[2].index()], u_min_1, u_max_2, 0, 1));
                        }


                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 0) //| overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 1)
                        {
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            int alpha = transparency;
                            c0.a = alpha;
                            c1.a = alpha;
                            c2.a = alpha;
                            drawTriangleColor2D(v0_u,v1_u,v2_u,c0, c1, c2);

                            alpha = 255;
                            c0.a = alpha;
                            c1.a = alpha;
                            c2.a = alpha;
                            ofColor bl = ofColor::black;
                            drawTriangleLinesColor2D(v0_u,v1_u,v2_u, bl, bl, bl);
                        }

                    }
                    else
                    {
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), v[0].point().z());
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), v[1].point().z());
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), v[2].point().z());

                        drawTetLines2D(v0_u,v1_u,v2_u);
                    }



                }
            }
            ofPopMatrix();
        }

        // DRAW SOLUTION ON OVERLAPPING MESH

        if(draw_second_mesh_u)
        {
            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            gl2psEnable(GL2PS_BLEND);
            ofPushMatrix();
            ofTranslate(0,0,trans_u_z);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(1)); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell(*overlapping_meshes->mesh(1), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {

                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    string indexstring;
                    if(!draw_indexes_mesh_vec.empty())
                    {
                        if(std::find(draw_indexes_mesh_vec.begin(), draw_indexes_mesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_mesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    dolfin::Point v0_u;
                    dolfin::Point v1_u;
                    dolfin::Point v2_u;
                    if(vert_values_2.size() > 0)
                    {
                        //v[2].entities(0)[0]][0]
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), zaxisscaler*vert_values_2[v[0].index()]);
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), zaxisscaler*vert_values_2[v[1].index()]);
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), zaxisscaler*vert_values_2[v[2].index()]);
                        ofColor r = firstcolor;
                        ofColor b = secondcolor;
                        ofColor c0;
                        ofColor c1;
                        ofColor c2;

                        if(u_min_1 == 0 && u_max_2 == 0)
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_2[v[0].index()], 0, 0.1, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_2[v[1].index()], 0, 0.1, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_2[v[2].index()], 0, 0.1, 0, 1));
                        }
                        else
                        {
                            c0 = b.getLerped(r, ofMap(vert_values_2[v[0].index()], u_min_1, u_max_2, 0, 1));
                            c1 = b.getLerped(r, ofMap(vert_values_2[v[1].index()], u_min_1, u_max_2, 0, 1));
                            c2 = b.getLerped(r, ofMap(vert_values_2[v[2].index()], u_min_1, u_max_2, 0, 1));
                        }

                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        if(draw_indexes) draw_or_save_string(x, y, indexstring);


                        int alpha = transparency;
                        c0.a = alpha;
                        c1.a = alpha;
                        c2.a = alpha;
                        drawTriangleColor2D(v0_u,v1_u,v2_u,c0, c1, c2);

                        alpha = 255;
                        c0.a = alpha;
                        c1.a = alpha;
                        c2.a = alpha;
                        ofColor bl = ofColor::black;
                        drawTriangleLinesColor2D(v0_u,v1_u,v2_u, bl, bl, bl);

                    }
                    else
                    {
                        v0_u = dolfin::Point(v[0].point().x(), v[0].point().y(), v[0].point().z());
                        v1_u = dolfin::Point(v[1].point().x(), v[1].point().y(), v[1].point().z());
                        v2_u = dolfin::Point(v[2].point().x(), v[2].point().y(), v[2].point().z());
                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        if(draw_indexes) draw_or_save_string(x, y, indexstring);

                        drawTetLines2D(v0_u,v1_u,v2_u);

                        string markstring;
                        size_t test = overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()];
                        markstring = ofToString(test);
                        draw_or_save_string(x, y, markstring);
                    }



                }
            }
            ofPopMatrix();
        }





        // DRAW SECOND MESH
        if(draw_second_mesh)
        {
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            ofSetColor(255,0,0,255);
            ofPushMatrix();
            ofTranslate(0,0,trans_second_mesh);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(1)); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell(*overlapping_meshes->mesh(1), cellit->global_index());
                dolfin::VertexIterator v(pickcell);
                if(mesh->geometry().dim() == 3)
                {
                    drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                }
                else
                {
                    float x = pickcell.midpoint().x();
                    float y = pickcell.midpoint().y();
                    string indexstring;
                    if(!draw_indexes_collidingMesh_vec.empty())
                    {
                        if(std::find(draw_indexes_collidingMesh_vec.begin(), draw_indexes_collidingMesh_vec.end(),
                                     cellit->global_index()) != draw_indexes_collidingMesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                    }
                    if(draw_indexes) draw_or_save_string(x, y, indexstring);
                    drawTetLines2D(v[0].point(),v[1].point(),v[2].point());
                }

            }
            ofPopMatrix();
        }

        // DRAW BORDER
        if(draw_border_first)
        {
            ofSetColor(0,255,0,255);
            glLineWidth(3.0);
            gl2psLineWidth(3.0);
            ofPushMatrix();
            ofTranslate(0,0,trans_border_first);
            if(overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->size(1) > 1)
            {
                for (dolfin::CellIterator cellit(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first); !cellit.end(); ++cellit)
                {
                    dolfin::Cell pickcell(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first, cellit->global_index());
                    dolfin::VertexIterator v(pickcell);

                    string indexstring;
                    if(!draw_indexes_facettriang_collidingMesh_vec.empty())
                    {
                        if(std::find(draw_indexes_facettriang_collidingMesh_vec.begin(), draw_indexes_facettriang_collidingMesh_vec.end(),
                                     (size_t) overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->parent_entity_maps()[0][cellit->global_index()]) != draw_indexes_facettriang_collidingMesh_vec.end())
                        {
                            indexstring += ofToString(cellit->global_index());
                            indexstring += "_{" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->parent_entity_maps()[0][cellit->global_index()]);
                            indexstring += ", " + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->parent_entity_maps()[1][cellit->global_index()]) + "}";

                        }
                    }
                    else
                    {
                        indexstring += ofToString(cellit->global_index());
                        indexstring += "_{" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->parent_entity_maps()[0][cellit->global_index()]);
                        indexstring += ", " + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(1).first->parent_entity_maps()[1][cellit->global_index()]) + "}";
                    }

                    float x = pickcell.midpoint().x();
                    float y = pickcell.midpoint().y();
                    if(draw_indexes) draw_or_save_string(x, y, indexstring);
                    drawLines2D(v[0].point(),v[1].point());
                    glPointSize(6.0);
                    gl2psPointSize(6.0);
                    drawBoxPoint2D(v[0].point());
                    drawBoxPoint2D(v[1].point());

                }
            }
            ofPopMatrix();
        }

        // DRAW TRIANGULATION OF T_(0,gamma)
        if(draw_triangulation_mesh)
        {
            glLineWidth(2.0);
            gl2psLineWidth(2.0);
            ofSetColor(0,0,0,255);
            ofPushMatrix();
            ofTranslate(0,0,trans_triangulation_mesh);
            if(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->size(2) > 1)
            {
                for (dolfin::CellIterator cellit(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first); !cellit.end(); ++cellit)
                {
                    dolfin::Cell pickcell(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first, cellit->global_index());
                    dolfin::VertexIterator v(pickcell);
                    if(overlapping_meshes->mesh(0)->geometry().dim() == 3)
                    {
                        drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                    }
                    else
                    {

                        string indexstring;
                        if(!draw_indexes_subtriang_mesh_vec.empty())
                        {
                            if(std::find(draw_indexes_subtriang_mesh_vec.begin(), draw_indexes_subtriang_mesh_vec.end(),
                                         (size_t) overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]) != draw_indexes_subtriang_mesh_vec.end())
                            {
                                indexstring += ofToString(cellit->global_index());
                                indexstring += "_" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]);
                            }
                        }
                        else
                        {
                            indexstring += ofToString(cellit->global_index());
                            indexstring += "_" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]);
                        }
                        float x = pickcell.midpoint().x();
                        float y = pickcell.midpoint().y();
                        if(draw_indexes) draw_or_save_string(x, y, indexstring);
                        drawTetLines2D(v[0].point(),v[1].point(),v[2].point());
                    }
                }
            }
            ofPopMatrix();
        }


        if(draw_triangulation_mesh_u)
        {
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            ofSetColor(0,0,0,255);
            ofPushMatrix();
            ofTranslate(0,0,trans_u_z);
            if(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->size(2) > 1)
            {
                for (dolfin::CellIterator cellit(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first); !cellit.end(); ++cellit)
                {
                    dolfin::Cell pickcell(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first, cellit->global_index());
                    dolfin::VertexIterator v(pickcell);
                    if(overlapping_meshes->mesh(0)->geometry().dim() == 3)
                    {
                        drawTetLines3D(v[0].point(),v[1].point(),v[2].point(),v[3].point());
                    }
                    else
                    {

                        string indexstring;
                        if(!draw_indexes_subtriang_mesh_vec.empty())
                        {
                            if(std::find(draw_indexes_subtriang_mesh_vec.begin(), draw_indexes_subtriang_mesh_vec.end(),
                                         (size_t) overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]) != draw_indexes_subtriang_mesh_vec.end())
                            {
                                indexstring += ofToString(cellit->global_index());
                                indexstring += "_" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]);
                            }
                        }
                        else
                        {
                            indexstring += ofToString(cellit->global_index());
                            indexstring += "_" + ofToString(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]);
                        }
                        dolfin::Point v_u[3];

                        if(vert_values_1.size() > 0)
                        {
                            dolfin::Cell pickcell_first_mesh(*overlapping_meshes->mesh(0), overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->parent_entity_maps()[0][cellit->global_index()]);
                            dolfin::VertexIterator v_first_mesh(pickcell_first_mesh);

                            //v[2].entities(0)[0]][0]
                            //v0_u = dolfin::Point(v_first_mesh[0].point().x(), v_first_mesh[0].point().y(), zaxisscaler*vert_values_1[v_first_mesh[0].index()]);
                            //v1_u = dolfin::Point(v_first_mesh[1].point().x(), v_first_mesh[1].point().y(), zaxisscaler*vert_values_1[v_first_mesh[1].index()]);
                            //v2_u = dolfin::Point(v_first_mesh[2].point().x(), v_first_mesh[2].point().y(), zaxisscaler*vert_values_1[v_first_mesh[2].index()]);
                            ofColor r = firstcolor;
                            ofColor b = secondcolor;

                            double ip_u[3];
                            ip_u[0] =  ofApp::interpolate_triangle(v_first_mesh[0].point(),
                                                                   vert_values_1[v_first_mesh[0].index()],
                                                                   v_first_mesh[1].point(),
                                                                   vert_values_1[v_first_mesh[1].index()],
                                                                   v_first_mesh[2].point(),
                                                                   vert_values_1[v_first_mesh[2].index()],
                                                                   v[0].point());
                            ip_u[1] =  ofApp::interpolate_triangle(v_first_mesh[0].point(),
                                                                   vert_values_1[v_first_mesh[0].index()],
                                                                   v_first_mesh[1].point(),
                                                                   vert_values_1[v_first_mesh[1].index()],
                                                                   v_first_mesh[2].point(),
                                                                   vert_values_1[v_first_mesh[2].index()],
                                                                   v[1].point());
                            ip_u[2] =  ofApp::interpolate_triangle(v_first_mesh[0].point(),
                                                                   vert_values_1[v_first_mesh[0].index()],
                                                                   v_first_mesh[1].point(),
                                                                   vert_values_1[v_first_mesh[1].index()],
                                                                   v_first_mesh[2].point(),
                                                                   vert_values_1[v_first_mesh[2].index()],
                                                                   v[2].point());



                            v_u[0] = dolfin::Point(v[0].point().x(), v[0].point().y(), zaxisscaler*ip_u[0]);
                            v_u[1] = dolfin::Point(v[1].point().x(), v[1].point().y(), zaxisscaler*ip_u[1]);
                            v_u[2] = dolfin::Point(v[2].point().x(), v[2].point().y(), zaxisscaler*ip_u[2]);

                            ofColor c0;
                            ofColor c1;
                            ofColor c2;
                            if(u_min_1 == 0 && u_max_2 == 0)
                            {
                                c0 = b.getLerped(r, ofMap(ip_u[v[0].index()], 0, 0.1, 0, 1));
                                c1 = b.getLerped(r, ofMap(ip_u[v[1].index()], 0, 0.1, 0, 1));
                                c2 = b.getLerped(r, ofMap(ip_u[v[2].index()], 0, 0.1, 0, 1));
                            }
                            else
                            {
                                c0 = b.getLerped(r, ofMap(ip_u[0], u_min_1, u_max_2, 0, 1));
                                c1 = b.getLerped(r, ofMap(ip_u[1], u_min_1, u_max_2, 0, 1));
                                c2 = b.getLerped(r, ofMap(ip_u[2], u_min_1, u_max_2, 0, 1));
                            }

                            float x = pickcell.midpoint().x();
                            float y = pickcell.midpoint().y();
                            if(draw_indexes) draw_or_save_string(x, y, indexstring);
                            int alpha = transparency;
                            c0.a = alpha;
                            c1.a = alpha;
                            c2.a = alpha;
                            drawTriangleColor2D(v_u[0],v_u[1],v_u[2],c0, c1, c2);

                            for (int i = 0; i < 2; i++)
                            {
                                for (int j = 0; j < 2; j++)
                                {
//                            std::cout << "check" << std::endl;
//                            std::cout << i << std::endl;
//                          std::cout << j << std::endl;
//                          std::cout << i+1%3 << std::endl;
//                          std::cout << j+1%3 << std::endl;

                                    if(point_on_line(v_first_mesh[i].point(), v_first_mesh[i+1%3].point(), v[j].point())
                                            && point_on_line(v_first_mesh[i].point(), v_first_mesh[i+1%3].point(), v[j+1%3].point()))
                                    {
                                        //if(do_intersect_two_lines(v[j].point(), v[j+1%2].point(), v_first_mesh[i].point(), v_first_mesh[i+1%2].point())) {
                                        alpha = 255;
                                        c0.a = alpha;
                                        c1.a = alpha;
                                        c2.a = alpha;
                                        ofColor bl = ofColor::black;
                                        ofSetColor(bl);
                                        drawLines2D(v_u[j],v_u[j+1%3]);
                                    }
                                }
                            }
                        }
                        else
                        {
                            v_u[0] = dolfin::Point(v[0].point().x(), v[0].point().y(), v[0].point().z());
                            v_u[1] = dolfin::Point(v[1].point().x(), v[1].point().y(), v[1].point().z());
                            v_u[2] = dolfin::Point(v[2].point().x(), v[2].point().y(), v[2].point().z());

                            drawTetLines2D(v_u[0],v_u[1],v_u[2]);
                        }

                    }
                }
            }
            ofPopMatrix();
        }
        // DRAW CELL MARKERS 0,1,2
        if(draw_index_triangulation)
        {
            ofSetColor(0,0,255,255);
            glLineWidth(1.0);
            gl2psLineWidth(1.0);
            gl2psEnable(GL2PS_BLEND);
            ofPushMatrix();
            ofTranslate(0,0,trans_triangulation_mesh);
            for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
            {
                dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());

                float x = pickcell.midpoint().x();
                float y = pickcell.midpoint().y();

                string markstring;
                size_t test = overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()];
                markstring = ofToString(test);
                draw_or_save_string(x, y, markstring);
            }
            ofPopMatrix();
        }
        cam.end();
    }
    if(!outputeps)
    {
        gui.draw();
        gui2.draw();
    }
}

void ofApp::drawTet(dolfin::Point v0,dolfin::Point v1,dolfin::Point v2, dolfin::Point v3)
{
    glBegin(GL_TRIANGLE_STRIP);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v1.x(), v1.y(), v1.z());
    glVertex3f(v2.x(), v2.y(), v2.z());
    glVertex3f(v3.x(), v3.y(), v3.z());
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v1.x(), v1.y(), v1.z());
    glEnd();
}

void ofApp::drawTetLines3D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2, dolfin::Point v3)
{

    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v1.x(), v1.y(), v1.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v3.x(), v3.y(), v3.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glVertex3f(v3.x(), v3.y(), v3.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v2.x(), v2.y(), v2.z());
    glVertex3f(v3.x(), v3.y(), v3.z());
    glEnd();

}

void ofApp::drawTetLines2D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2)

{

    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v1.x(), v1.y(), v1.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();
}

void ofApp::drawTriangleColor2D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2, ofColor c0, ofColor c1, ofColor c2)

{
    glBegin(GL_TRIANGLE_STRIP);
    glColor4f(c0.r/255.0, c0.g/255.0, c0.b/255.0, c0.a/255.0);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glColor4f(c1.r/255.0, c1.g/255.0, c1.b/255.0, c1.a/255.0);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glColor4f(c2.r/255.0, c2.g/255.0, c2.b/255.0, c2.a/255.0);
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();
}

void ofApp::drawTriangleLinesColor2D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2, ofColor c0, ofColor c1, ofColor c2)

{
    glBegin(GL_LINES);
    glColor4f(c0.r/255.0, c0.g/255.0, c0.b/255.0, c0.a/255.0);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glColor4f(c1.r/255.0, c1.g/255.0, c1.b/255.0, c1.a/255.0);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glEnd();

    glBegin(GL_LINES);
    glColor4f(c0.r/255.0, c0.g/255.0, c0.b/255.0, c0.a/255.0);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glColor4f(c2.r/255.0, c2.g/255.0, c2.b/255.0, c2.a/255.0);
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();

    glBegin(GL_LINES);
    glColor4f(c1.r/255.0, c1.g/255.0, c1.b/255.0, c1.a/255.0);
    glVertex3f(v1.x(), v1.y(), v1.z());
    glColor4f(c2.r/255.0, c2.g/255.0, c2.b/255.0, c2.a/255.0);
    glVertex3f(v2.x(), v2.y(), v2.z());
    glEnd();
}

void ofApp::drawLines2D(dolfin::Point v0,dolfin::Point v1)

{
    glBegin(GL_LINES);
    glVertex3f(v0.x(), v0.y(), v0.z());
    glVertex3f(v1.x(), v1.y(), v1.z());
    glEnd();
}

void ofApp::drawBoxPoint2D(dolfin::Point v0)
{
    glBegin(GL_POINTS);
    glVertex3f(v0.x(),v0.y(),v0.z());
    glEnd();
}


//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{

    if(key == ' ')
    {
        //triangulation.refinemeshborder();
        checksets = true;
    }
    if(key == 'm')
    {
    timing_bvh.clear();
    timing_triang.clear();
    timing_wo_cg.clear();
    timing_total.clear();
    size_mesh1.clear();
    size_mesh2.clear();
    size_triangulation.clear();
    size_interface.clear();

    for (int i = 0; i < 4; i++) {
        newStructuredMeshes(80 + i*80);
        runPoisson();
        runPoissonOLM();
    }

    std::cout << "\\addplot[stack plots=y, area style,black!80!white,fill=black!80!white, mark=square*] coordinates {";
    for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(timing_wo_cg[i].first) + ", " + ofToString(timing_wo_cg[i].second) + ") ";
    }
    std::cout << "} \\closedcycle;" << std::endl;
    std::cout << "\\addplot[stack plots=y, area style,black!50!white,fill=black!50!white, mark=square*] coordinates {";
    for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(timing_triang[i].first) + ", " + ofToString(timing_triang[i].second) + ") ";
    }
    std::cout << "} \\closedcycle;" << std::endl;
    std::cout << "\\addplot[stack plots=y, area style,black!30!white,fill=black!30!white, mark=square*] coordinates {";
    for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(timing_bvh[i].first) + ", " + ofToString(timing_bvh[i].second) + ") ";
    }
    std::cout << "} \\closedcycle;" << std::endl;
    std::cout << "\\addplot[black!30!white, mark=square, mark options=solid, dashed] coordinates {";
    for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(timing_standard_poisson[i].first) + ", " + ofToString(timing_standard_poisson[i].second) + ") ";
    }
    std::cout << "};" << std::endl;
    std::cout << "\\addplot[black!100!white,fill=black!80!white] coordinates {";
        for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(size_mesh1[i].first) + ", " + ofToString(size_mesh1[i].second) + ") ";
    }
    std::cout << "};" << std::endl;
    std::cout << "\\addplot[black!100!white,fill=black!40!white] coordinates {";
        for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(size_mesh2[i].first) + ", " + ofToString(size_mesh2[i].second) + ") ";
    }
    std::cout << "};" << std::endl;
    std::cout << "\\addplot[black!100!white,fill=black!20!white] coordinates {";
        for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(size_triangulation[i].first) + ", " + ofToString(size_triangulation[i].second) + ") ";
    }
    std::cout << "};" << std::endl;
    std::cout << "\\addplot[black!100!white,fill=black!0!white] coordinates {";
    for (int i = 0; i < timing_bvh.size(); i++) {
        std::cout << "(" + ofToString(size_interface[i].first) + ", " + ofToString(size_interface[i].second) + ") ";
    }
    std::cout << "};" << std::endl;
    }
    if(key == 'd')
    {
        drawsomething = !drawsomething;
    }

    if(key == 'x')
    {
        //triangulation.refinemesh();
        //triangulation.markborder();
        checksets = true;
    }
    if(key == 'c')
    {
        //triangulation.markborder();
        checksets = true;
    }

    if(key == 'e')
    {

        if (true)
        {
            outputeps = true;
            FILE *fp = fopen("data/myfile.pdf", "wb");
            GLint buffsize = 0, state = GL2PS_OVERFLOW;
            GLint viewport[4];

            glGetIntegerv(GL_VIEWPORT, viewport);

            while( state == GL2PS_OVERFLOW )
            {
                buffsize += 1024*1024;
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, svg_n_slider, svg_n_slider, svg_n_slider, buffsize,
                                 fp, "data/myfile.pdf" );

                gl2psEnable(GL2PS_BLEND);
                draw();


                state = gl2psEndPage();
            }



            fclose(fp);
            //system("ps2eps -f data/myfile data/myfile");
            outputeps = false;
        }

        if (true)
        {
            outputeps = true;
            FILE *fp = fopen("data/myfile.svg", "wb");
            GLint buffsize = 0, state = GL2PS_OVERFLOW;
            GLint viewport[4];

            glGetIntegerv(GL_VIEWPORT, viewport);

            while( state == GL2PS_OVERFLOW )
            {
                buffsize += 1024*1024;
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_SVG, GL2PS_BSP_SORT, GL2PS_SILENT | GL2PS_BEST_ROOT | GL2PS_SIMPLE_LINE_OFFSET,
                                 GL_RGBA, 0, NULL, svg_n_slider, svg_n_slider, svg_n_slider, buffsize,
                                 fp, "data/myfile.svg" );
                /*
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT |
                                 GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
                                 GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                 fp, "data/myfile" );
                                 */
                gl2psEnable(GL2PS_BLEND);
                draw();


                state = gl2psEndPage();
            }



            fclose(fp);

            outputeps = false;
        }
        if(true)
        {
            outputtex = true;
            FILE *fp = fopen("data/myfile.tex", "wb");
            GLint buffsize = 0, state = GL2PS_OVERFLOW;
            GLint viewport[4];

            glGetIntegerv(GL_VIEWPORT, viewport);

            while( state == GL2PS_OVERFLOW )
            {
                buffsize += 1024*1024;
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_TEX, GL2PS_SIMPLE_SORT, GL2PS_SILENT |
                                 GL2PS_NO_BLENDING |
                                 GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                 fp, "data/myfile.tex" );
                /*
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT |
                                 GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
                                 GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                 fp, "data/myfile" );
                                 */
                draw();


                state = gl2psEndPage();
            }



            fclose(fp);
            system("./data/bashclip.sh");
            outputtex = false;
        }
        if(true)
        {
            outputtex = true;
            outputeps = true;
            FILE *fp = fopen("data/myfile.eps", "wb");
            GLint buffsize = 0, state = GL2PS_OVERFLOW;
            GLint viewport[4];

            glGetIntegerv(GL_VIEWPORT, viewport);

            while( state == GL2PS_OVERFLOW )
            {
                buffsize += 1024*1024;
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_EPS, GL2PS_SIMPLE_SORT, GL2PS_SILENT |
                                 GL2PS_NO_BLENDING |
                                 GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                 fp, "data/myfile.eps" );
                /*
                gl2psBeginPage ( "Mesh", "DOLFIN_OF_TOMNAERLAND", viewport,
                                 GL2PS_PDF, GL2PS_BSP_SORT, GL2PS_SILENT |
                                 GL2PS_SIMPLE_LINE_OFFSET | GL2PS_NO_BLENDING |
                                 GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT,
                                 GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
                                 fp, "data/myfile" );
                                 */
                draw();


                state = gl2psEndPage();
            }



            fclose(fp);
            //system("ps2eps -f data/myfile data/myfile");
            outputeps = true;
            outputtex = false;
        }
    }
    if(key == 'w')
    {
        testVolume();
    }

    if(key == 'l')
    {
        dolfin::File mesh_file("gendata/vesselaneurysmSaved.xml");
        mesh_file << *mesh;
    }

    if(key == 'b')
    {
        savedPose = cam.getGlobalTransformMatrix();
        cout << savedPose << endl;
    }



    if(key == 'k')
    {
        cout << cam.getNearClip() << endl;
        cam.setNearClip(cam.getNearClip()/10.0);
        cout << cam.getFarClip() << endl;
        cam.setFarClip(cam.getFarClip()*2.0);
    }
    if(key == 'b')
    {
        runPoissonOLM_patch_test_p1();
    }
    if(key == 'v')
    {
        runPoissonOLM_patch_test_p2();
    }
    if(key == 't')
    {
        runPoissonOLM_convergence_p1();
    }

    if(key == 'z')
    {
        n_structured++;
        newStructuredMeshes(pow(n_structured,2));
        //newStructuredMeshes(n_structured);
    }

    if(key == 'n')
    {
        runPoissonOLM_convergence_p2();
    }

    if(key == 'u')
    {
        refineMeshes();
    }
    if(key == 'p')
    {
        runPoissonOLM();
    }
    if(key == 'o')
    {
        runPoisson();
    }
    if(key == 'i')
    {
        runPoissonNitscheBC();
    }

    if(key == 'g')
    {
        t_0_strings.clear();
        t_2_strings.clear();
        combined_strings.clear();

        for(int i = 0; i < 8; i++)
        {
            std::cout << "CURRENT STEP" << std::endl;
            std::cout << i << std::endl;
            //newStructuredMeshes(pow(2,i));
            refineMeshes();
            refineMeshes();
            runPoissonOLM_convergence_p1();
        }
        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings[i].first) + ", " + ofToString(t_0_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings[i].first) + ", " + ofToString(t_2_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings[i].first) + ", " + ofToString(combined_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings.size(); i++)
        {

            string eoc = " ";
            if(i < combined_strings.size() -1)
            {
                eoc = ofToString((log(combined_strings[i+1].second) - log(combined_strings[i].second))/(log(combined_strings[i+1].first) - log(combined_strings[i].first)));
            }
            std::cout << ofToString(combined_strings[i].first) + "& " + ofToString(combined_strings[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;

    }

    if(key == 'f')
    {
        t_0_strings.clear();
        t_2_strings.clear();
        combined_strings.clear();

        for(int i = 0; i < 4; i++)
        {
            std::cout << "CURRENT STEP" << std::endl;
            std::cout << i << std::endl;
            //newStructuredMeshes(pow(2,i));
            refineMeshes();
            refineMeshes();
            runPoissonOLM_convergence_p2();
        }
        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings[i].first) + ", " + ofToString(t_0_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings[i].first) + ", " + ofToString(t_2_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings[i].first) + ", " + ofToString(combined_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings.size(); i++)
        {

            string eoc = " ";
            if(i < combined_strings.size() -1)
            {
                eoc = ofToString(std::abs((std::log(combined_strings[i+1].second) - std::log(combined_strings[i].second)))/std::abs(std::log(combined_strings[i+1].first) - std::log(combined_strings[i].first)));
            }
            std::cout << ofToString(combined_strings[i].first) + "& " + ofToString(combined_strings[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;

    }

    if(key == 'c')
    {
        t_0_strings.clear();
        t_2_strings.clear();
        combined_strings.clear();
        t_0_strings_p2.clear();
        t_2_strings_p2.clear();
        combined_strings_p2.clear();


        for(int i = 1; i < 5; i++)
        {
            std::cout << "CURRENT STEP" << std::endl;
            std::cout << i << std::endl;

            //newStructuredMeshes(pow(i,2));
            refineMeshes();
            refineMeshes();
            runPoissonOLM_convergence_p1();
            runPoissonOLM_convergence_p2();

            //newStructuredMeshes(i*2);
            //newStructuredMeshes(pow(2,i));

            /*
            refineMeshes();

            testVolume();
            runPoissonOLM_convergence_p2();
            */
                /*
            newStructuredMeshes(i*2);
            //pow(2,i)
            //refineMeshes();

            */
        }
        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings[i].first) + ", " + ofToString(t_0_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings[i].first) + ", " + ofToString(t_2_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings[i].first) + ", " + ofToString(combined_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings.size(); i++)
        {

            string eoc = " ";
            if(i > 0)
            {
                eoc = ofToString((log(combined_strings[i].second) - log(combined_strings[i-1].second))/(log(combined_strings[i].first) - log(combined_strings[i-1].first)));
            }
            std::cout << ofToString(combined_strings[i].first) + "& " + ofToString(combined_strings[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;


        ///


        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings_p2[i].first) + ", " + ofToString(t_0_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings_p2[i].first) + ", " + ofToString(t_2_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings_p2[i].first) + ", " + ofToString(combined_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings_p2.size(); i++)
        {

            string eoc = " ";
            if(i > 0)
            {
                eoc = ofToString((log(combined_strings_p2[i].second) - log(combined_strings_p2[i-1].second))/(log(combined_strings_p2[i].first) - log(combined_strings_p2[i-1].first)));
            }
            std::cout << ofToString(combined_strings_p2[i].first) + "& " + ofToString(combined_strings_p2[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;


    }

    if(key == 'x')
    {
        t_0_strings.clear();
        t_2_strings.clear();
        combined_strings.clear();
        t_0_strings_p2.clear();
        t_2_strings_p2.clear();
        combined_strings_p2.clear();


        for(int i = 0; i < 7; i++)
        {
            std::cout << "CURRENT STEP" << std::endl;
            std::cout << i << std::endl;

            newStructuredMeshes(pow(2,i));
            //refineMeshes();
            //refineMeshes();
            runPoissonOLM_convergence_p1();
            runPoissonOLM_convergence_p2();

            //newStructuredMeshes(i*2);
            //newStructuredMeshes(pow(2,i));

            /*
            refineMeshes();

            testVolume();
            runPoissonOLM_convergence_p2();
            */
                /*
            newStructuredMeshes(i*2);
            //pow(2,i)
            //refineMeshes();

            */
        }
        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings[i].first) + ", " + ofToString(t_0_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings[i].first) + ", " + ofToString(t_2_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings[i].first) + ", " + ofToString(combined_strings[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings.size(); i++)
        {

            string eoc = " ";
            if(i > 0)
            {
                eoc = ofToString((log(combined_strings[i].second) - log(combined_strings[i-1].second))/(log(combined_strings[i].first) - log(combined_strings[i-1].first)));
            }
            std::cout << ofToString(combined_strings[i].first) + "& " + ofToString(combined_strings[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;


        ///


        std::cout << "###########################GRAPH###########################" << std::endl;
        std::cout << "\\begin{figure}[!h]\\centering\\begin{tikzpicture}\\begin{axis}[ xmode=log, ymode=log, log basis x=10, xlabel=$h_{max}$, ylabel= $\\text{H}^1 \\text{norm}$, legend style={at={(0.34,0.92)}}, width=130mm, height=80mm, max space between ticks=90pt]" << std::endl;
        std::cout << "\\addplot[color=blue,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_0_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_0_strings_p2[i].first) + ", " + ofToString(t_0_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_0 + \\mathcal{T}_\\gamma$}" << std::endl;
        std::cout << "\\addplot[color=red,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < t_2_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(t_2_strings_p2[i].first) + ", " + ofToString(t_2_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{$\\mathcal{T}_2$}" << std::endl;
        std::cout << "\\addplot[color=black,mark=o]coordinates {" << std::endl;
        for(int i = 0; i < combined_strings_p2.size(); i++)
        {
            std::cout <<  "(" + ofToString(combined_strings_p2[i].first) + ", " + ofToString(combined_strings_p2[i].second) + ")" << std::endl;
        }
        std::cout << "};\\addlegendentry{combined}" << std::endl;
        std::cout << "\\end{axis}\\end{tikzpicture}\\caption{Convergence in $H^1$ norm} \\label{adv_plot}\\end{figure}" << std::endl;


        std::cout << "###########################TABLE###########################" << std::endl;
        std::cout << "\\begin{table}[h!]\\small\\noindent\\makebox[\\textwidth]{\\begin{tabular}{|r |r |r|}\\hline $h_{max}$ & $ ||\\cdot||_{\\text{H}_1} $ \\\\ \\hline \\hline" << std::endl;
        /*
        for(int i = 0; i < t_0_strings.size(); i++) {
        std::cout << ofToString(t_0_strings[i].first) + "& " + ofToString(t_0_strings[i].second) + "\\\\ " << std::endl;
        }
        for(int i = 0; i < t_2_strings.size(); i++) {
        std::cout << ofToString(t_2_strings[i].first) + "& " + ofToString(t_2_strings[i].second) + "\\\\ " << std::endl;
        }
        */
        for(int i = 0; i < combined_strings_p2.size(); i++)
        {

            string eoc = " ";
            if(i > 0)
            {
                eoc = ofToString((log(combined_strings_p2[i].second) - log(combined_strings_p2[i-1].second))/(log(combined_strings_p2[i].first) - log(combined_strings_p2[i-1].first)));
            }
            std::cout << ofToString(combined_strings_p2[i].first) + "& " + ofToString(combined_strings_p2[i].second) + "& "
                      + eoc + "\\\\ " << std::endl;
        }
        std::cout << "\\hline\\end{tabular}}\\caption{Convergence}\\end{table}\\normalsize" << std::endl;

    }


    //setNearClip(...)

}

void ofApp::newStructuredMeshes(size_t n)
{
    size_t N = n;
    // Build background mesh
    double a = 3.0;
    double b = -3.0;
    mesh->clear();
    collidingMesh->clear();
    std::shared_ptr<dolfin::Mesh> mesh0(
        new dolfin::RectangleMesh(a, a, b, b, 3 * N, 3 * N));
    mesh0->init(mesh0->topology().dim() - 1, mesh0->topology().dim());

    // Create sub mesh
    std::shared_ptr<dolfin::CellFunction<std::size_t> > domain_marker_0
    (new dolfin::CellFunction<std::size_t>(*mesh0, 0));

    InnerDomain().mark(*domain_marker_0, 2);
    std::shared_ptr<dolfin::SubMesh> mesh1(
        new dolfin::SubMesh(*mesh0, *domain_marker_0, 2));

    double scalefactor = 1.5;
    dolfin::MeshGeometry& geometry = mesh1->geometry();
    for (std::size_t i = 0; i < geometry.size(); i++)
    {
        // Get coordinate
        double* x = geometry.x(i);
        // Scale
        const double x0 = x[0]*scalefactor;
        const double x1 = x[1]*scalefactor;
        // Store coordinate
        x[0] = x0;
        x[1] = x1;
    }
    mesh = new dolfin::Mesh(*mesh0);
    collidingMesh = new dolfin::Mesh(*mesh1);

    collidingMesh->rotate(40,2);
    mesh->init();
    collidingMesh->init();

    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes_new(new cutfem::TriangulatedOverlappingMeshes);
    overlapping_meshes = overlapping_meshes_new;

    std::shared_ptr<dolfin::Mesh> mesh_mesh(new dolfin::Mesh(*mesh));
    std::shared_ptr<dolfin::Mesh> mesh_collidingMesh(new dolfin::Mesh(*collidingMesh));
    overlapping_meshes->add(mesh_mesh);
    overlapping_meshes->add(mesh_collidingMesh);

    // Computing intersections
    overlapping_meshes->compute_collisions();
    overlapping_meshes->compute_intersections();

}


void ofApp::refineMeshes()
{
    dolfin::CellFunction<bool> cell_markers(*mesh);
    cell_markers.set_all(true);
    dolfin::Mesh refined = dolfin::refine(*mesh, cell_markers);
    refined.init();
    mesh = new dolfin::Mesh(refined);

    dolfin::CellFunction<bool> cell_markers_c(*collidingMesh);
    cell_markers_c.set_all(true);
    dolfin::Mesh refined_c = dolfin::refine(*collidingMesh, cell_markers_c);
    refined_c.init();
    collidingMesh = new dolfin::Mesh(refined_c);


    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes_new(new cutfem::TriangulatedOverlappingMeshes);
    overlapping_meshes = overlapping_meshes_new;

    std::shared_ptr<dolfin::Mesh> mesh_mesh(new dolfin::Mesh(*mesh));
    std::shared_ptr<dolfin::Mesh> mesh_collidingMesh(new dolfin::Mesh(*collidingMesh));
    overlapping_meshes->add(mesh_mesh);
    overlapping_meshes->add(mesh_collidingMesh);

    // Computing intersections
    dolfin::Timer timer_domain_decomp("Total time domain decomposition");
    timer_domain_decomp.start();
    overlapping_meshes->compute_collisions();
    overlapping_meshes->compute_intersections();
    timer_domain_decomp.stop();

}

void ofApp::testVolume()
{
    double volume_mesh1;
    for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
    {
        dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
        volume_mesh1 += pickcell.volume();
    }
    double volume_mesh2;
    for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(1)); !cellit.end(); ++cellit)
    {
        dolfin::Cell pickcell(*overlapping_meshes->mesh(1), cellit->global_index());
        volume_mesh2 += pickcell.volume();
    }


    double volume_t_0;
    for (dolfin::CellIterator cellit(*overlapping_meshes->mesh(0)); !cellit.end(); ++cellit)
    {

        if(overlapping_meshes->cell_marker(0).get()[0][cellit->global_index()] == 0)
        {
            dolfin::Cell pickcell(*overlapping_meshes->mesh(0), cellit->global_index());
            volume_t_0 += pickcell.volume();
        }
    }

    double volume_local_triangulation;
    if(overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first->size(2) > 1)
    {
        for (dolfin::CellIterator cellit(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first); !cellit.end(); ++cellit)
        {
            dolfin::Cell pickcell(*overlapping_meshes->cut_mesh_and_parent_mesh_ids(0).first, cellit->global_index());
            dolfin::VertexIterator v(pickcell);
            volume_local_triangulation += pickcell.volume();
        }
    }

    double sum_untriangulated = volume_mesh1 - volume_mesh2;
    cout << "Sum: mesh - collidingMesh" << endl;
    std::cout << std::setprecision (30) << sum_untriangulated << std::endl;

    cout << "Sum: submesh T_0 + local triangulation" << endl;
    double sum_triangulated = volume_t_0 + volume_local_triangulation;
    std::cout << std::setprecision (30) << sum_triangulated << std::endl;

    double sum_difference = sum_untriangulated - sum_triangulated;
    cout << "Sum: difference" << endl;
    cout << std::setprecision (30)<< std::fixed << sum_difference << endl;

    double precision = 5.96e-14;
    cout << "Float precision" << endl;
    cout << std::setprecision (30)<< std::fixed << precision << endl;

}

void ofApp::draw_or_save_string(float x, float y, string indexstring)
{
    if(indexstring != "")
    {
        if(!outputeps)
        {
            ofDrawBitmapString(indexstring, x, y);
        }
        if(outputtex)
        {
            string texstring = "$" + indexstring + "$";
            glRasterPos2d(x, y);
            const char *fonts[] =
            {
                "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic",
                "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique",
                "Courier", "Courier-Bold", "Courier-Oblique", "Courier-BoldOblique",
                "Symbol", "ZapfDingbats"
            };
            const char *cstr = texstring.c_str();
            gl2psText(cstr, fonts[0], 16);
        }
    }
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key)
{
    outputeps = false;
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y )
{

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button)
{

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button)
{

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button)
{

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h)
{

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg)
{

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo)
{

}


void ofApp::runPoisson()
{
    dolfin::Timer timer_poisson("Standard Poisson");
    timer_poisson.start();
    int startsolve = ofGetElapsedTimeMillis();
    using namespace cutfem;
    {
        using namespace dolfin;
// Source term (right-hand side)
        class Source : public Expression
        {
            void eval(Array<double>& values, const Array<double>& x) const
            {
                float k_0 = 1;
                float k_1 = 1;
                values[0] = DOLFIN_PI*DOLFIN_PI*(k_0*k_0*sin((DOLFIN_PI)*x[0])*sin(DOLFIN_PI*x[1]) + k_1*k_1*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
            }
        };
        // Sub domain for Dirichlet boundary condition
        class DirichletBoundary : public SubDomain
        {
            bool inside(const Array<double>& x, bool on_boundary) const
            {
                return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;

            }
        };

        // Create mesh and function space
        Poisson::FunctionSpace V(*mesh);

        // Define boundary condition
        Constant u0(0.0);

        DirichletBoundary boundary;
        DirichletBC bc(V, u0, boundary);

        // Define variational forms
        Poisson::BilinearForm a(V, V);
        Poisson::LinearForm L(V);

        Source f;
        L.f = f;
        Constant gamma(gamma_u);
        Constant g(0);

        // Compute solution
        Function u(V);
        solve(a == L, u, bc);
        timer_poisson.stop();

        // Plot solution
        std::vector<dolfin::la_index> dof_vertex_set = vertex_to_dof_map(*u.function_space());
        std::vector<double> values2(dof_vertex_set.size());
        u.vector()->get_local(values2.data(), dof_vertex_set.size(), dof_vertex_set.data());
        vert_values_1 = values2;
        u_min_1 = (u.vector()->min());
        u_max_2 = (u.vector()->max());

        draw_t_0_1 = false;
        draw_first_mesh_t_0_u = false;
        draw_triangulation_mesh_u = false;
        draw_second_mesh_u = false;
        draw_first_mesh_u = true;
        draw_first_mesh = true;
        draw_second_mesh = false;
        draw_triangulation_mesh = false;
        draw_indexes = false;
        draw_index_triangulation = false;
        draw_border_first = false;
    }



    double standard_poisson_t = dolfin::timing("Standard Poisson");
    double combined_size = overlapping_meshes->mesh(0)->num_cells() + overlapping_meshes->mesh(1)->num_cells();
    timing_standard_poisson.push_back(make_pair(combined_size, standard_poisson_t));

    dolfin::list_timings(true);
    int endsolve = ofGetElapsedTimeMillis() - startsolve;
    std::cout << "time to solve standard Poisson " + ofToString(endsolve) + "ms"<< std::endl;
}

void ofApp::runPoissonNitscheBC()
{
    /*
    gamma_u = gamma_u + 0.1;
    if(gamma_u > 10)
    {
        gamma_u = 0.5;
    }
    */
    int startsolve = ofGetElapsedTimeMillis();
    using namespace cutfem;
    {
        using namespace dolfin;

        // Source term (right-hand side)
        class Source : public Expression
        {
            void eval(Array<double>& values, const Array<double>& x) const
            {
                float k_0 = 1;
                float k_1 = 1;
                values[0] = DOLFIN_PI*DOLFIN_PI*(k_0*k_0*sin((DOLFIN_PI)*x[0])*sin(DOLFIN_PI*x[1]) + k_1*k_1*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
            }
        };

        // Create mesh and function space
        PoissonNitscheBC::FunctionSpace V(*mesh);

        // Define variational forms
        PoissonNitscheBC::BilinearForm a(V, V);
        PoissonNitscheBC::LinearForm L(V);

        Source f;
        L.f = f;
        Constant gamma(gamma_u);
        Constant g(0);
        L.g = g;
        L.gamma = gamma;
        a.gamma = gamma;

        // Compute solution
        Function u(V);

        solve(a == L, u);


        // Plot solution

        std::vector<dolfin::la_index> dof_vertex_set = vertex_to_dof_map(*u.function_space());
        std::vector<double> values2(dof_vertex_set.size());
        u.vector()->get_local(values2.data(), dof_vertex_set.size(), dof_vertex_set.data());
        vert_values_1 = values2;
        u_min_1 = (u.vector()->min());
        u_max_2 = (u.vector()->max());

        draw_t_0_1 = false;
        draw_first_mesh_t_0_u = false;
        draw_triangulation_mesh_u = false;
        draw_second_mesh_u = false;
        draw_first_mesh_u = true;
        draw_first_mesh = true;
        draw_second_mesh = false;
        draw_triangulation_mesh = false;
        draw_indexes = false;
        draw_index_triangulation = false;
        draw_border_first = false;

    }
    int endsolve = ofGetElapsedTimeMillis() - startsolve;
    std::cout << "time to solve nitscheBC Poisson " + ofToString(endsolve) + "ms"<< std::endl;
}

void ofApp::runPoissonOLM()
{
    int startsolve = ofGetElapsedTimeMillis();
    dolfin::Timer timer_total("Total time inside OLM call");
    timer_total.start();
    std::clock_t start;
    double duration;
    start = std::clock();

    overlapping_meshes->compute_intersections();

    dolfin::Timer timer_solveassembly("Time without computational geometry");
    timer_solveassembly.start();
    // Reuse mesh 0 from MatchingCompositeMesh
    dolfin::Matrix A_2;

    // Create function spaces on each mesh and composite function space
    PoissonOLM0::FunctionSpace V0(overlapping_meshes->mesh(0));
    PoissonOLM1::FunctionSpace V1(overlapping_meshes->mesh(1));

    cutfem::CompositeFunctionSpace V_c(overlapping_meshes);
    V_c.add(V0);
    V_c.add(V1);
    V_c.build();

    V_c.view(0);

    // Create empty CompositeForm
    cutfem::CompositeForm a_c(V_c, V_c);

    // Create standard bilinear forms on mesh 0
    PoissonOLM0::BilinearForm a0(V0, V0);
    dolfin::Constant gamma(gamma_u);
    a0.gamma = gamma;


    // Get domain marker for mesh 0;
    auto cell_marker_0 = overlapping_meshes->cell_marker(0);

    // Attach data to standard forms
    a0.set_cell_domains(cell_marker_0);

    // Compute constrained dofs for V0
    auto & constrained_dofs = V_c.constrained_dofs();

    // (Candidates are those ones which lives in elements marked with 2).
    cutfem::CutFEMTools::compute_constrained_dofs(constrained_dofs,
            *V_c.view(0),
            cell_marker_0.get(),
            2);

    std::size_t dof_counter = 0;
    for (std::size_t i = 0; i < constrained_dofs.size(); ++i)
        if (constrained_dofs[i])
        {
            ++dof_counter;
            //dolfin::info("Dof %d is constrained.", i);
        }
    dolfin::info("In total, %d are constrained", dof_counter);

    cutfem::CompositeFunction vis_constr_dofs(V_c);
    cutfem::CutFEMTools::visualise_constrained_dofs(vis_constr_dofs,
            constrained_dofs);

    dolfin::File("constrained_dofs.pvd") << vis_constr_dofs.part(0);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a0)));

    // Get the cut mesh describing cut elements in mesh 0
    auto cut_mesh_and_parent_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(0);
    auto cut_mesh_0 = cut_mesh_and_parent_ids.first;
    auto parent_mesh_ids_0 = cut_mesh_and_parent_ids.second;

    //DEBUG
    dolfin::File("cut_mesh_0.pvd") << *cut_mesh_0;

    dolfin::info("parent_mesh_ids_0 = %d", parent_mesh_ids_0.size());
    // Create quadrature rule for mesh 0
    std::size_t order = 1;
    std::shared_ptr<cutfem::Quadrature> quadrature_a0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    a_c.cut_form(0)->set_quadrature(1, quadrature_a0_domain_1);
    a_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    a_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::BilinearForm a1(V1, V1);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a1)));

    //Add interface form
    // FIXME: Passed function space is just dummy.
    PoissonOLMInterface::BilinearForm ai(V0, V0);
    ai.gamma = gamma;

    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(ai)));

    // Get interface cut mesh and related information, stored as the second
    // cut_mesh_and_parent_ids component in TomOverlappingMeshes
    auto interface_mesh_and_parent_mesh_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(1);
    auto interface_mesh = interface_mesh_and_parent_mesh_ids.first;
    auto interface_parent_mesh_ids = interface_mesh_and_parent_mesh_ids.second;

    //DEBUG
    dolfin::File("interface_mesh.pvd") << *interface_mesh;

    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_interface(
        new cutfem::Quadrature(interface_mesh->type().cell_type(),
                               interface_mesh->geometry().dim(), order));

    // Build normal field for quadrature points
    std::shared_ptr<cutfem::FacetNormals> interface_facet_normals(
        new cutfem::FacetNormals(interface_mesh->facet_normals(),
                                 interface_mesh->geometry().dim(),
                                 quadrature_interface->size()));

    a_c.cut_form(2)->set_quadrature(0, quadrature_interface);
    a_c.cut_form(2)->set_cut_mesh(0, interface_mesh);
    a_c.cut_form(2)->set_parent_mesh_ids(0, interface_parent_mesh_ids);
    a_c.cut_form(2)->set_facet_normals(0, interface_facet_normals);


    // Assemble composite form and compare.
    cutfem::CompositeFormAssembler assembler;
    dolfin::Timer timer_A2("Assembling A_2");
    timer_A2.start();
    assembler.assemble(A_2, a_c);
    timer_A2.stop();

    // Temporary assemble of rhs and solution.
    cutfem::CompositeForm L_c(V_c);

    //dolfin::Constant f(1);

    // Source term (right-hand side)
    class Source : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            float k_0 = 1;
            float k_1 = 1;
            values[0] = DOLFIN_PI*DOLFIN_PI*(k_0*k_0*sin((DOLFIN_PI)*x[0])*sin(DOLFIN_PI*x[1]) + k_1*k_1*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
        }
    };

    class DirichletBC : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = 0;
        }
    };

    PoissonOLM0::LinearForm L0(V0);

    // Attach data to standard forms
    Source f;
    L0.f = f;

    DirichletBC g;
    L0.g = g;

    L0.gamma = gamma;

    L0.set_cell_domains(cell_marker_0);
    //L0.f = f;

    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L0)));

    // Get the cut mesh describing cut elements in mesh 0

    // Create quadrature rule for mesh 0
    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_L0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    L_c.cut_form(0)->set_quadrature(1, quadrature_L0_domain_1);
    L_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    L_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::LinearForm L1(V1);
    //L1.f = f;
    class Source_overlapping : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            float k_0 = 1;
            float k_1 = 1;
            values[0] = DOLFIN_PI*DOLFIN_PI*(k_0*k_0*sin((DOLFIN_PI)*x[0])*sin(DOLFIN_PI*x[1]) + k_1*k_1*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
        }
    };
    Source_overlapping f_2;
    L1.f = f_2;
    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L1)));

    dolfin::Vector b_2;

    dolfin::Timer timer_b2("Assembling b_2");
    timer_b2.start();
    assembler.assemble(b_2, L_c);
    timer_b2.stop();


    cutfem::CompositeFunction u(V_c);

    dolfin::Timer timer_solve("Solving system");
    timer_solve.start();
    dolfin::solve(A_2, *u.vector(), b_2);
    timer_solve.stop();

    timer_solveassembly.stop();
    timer_total.stop();
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<< "clock time Total " << duration <<'\n';

    double bvh_t = dolfin::timing("Bounding volume hierarchy, construction and traversal");
    double tot_triang_t = dolfin::timing("Total triangulation of mesh outside interface");
    double total_wo_cg_t = dolfin::timing("Time without computational geometry");
    double total_t = dolfin::timing("Total time inside OLM call");
    double triang_t = tot_triang_t - bvh_t;

    std::cout << "bvh " + ofToString(bvh_t) << std::endl;
    std::cout << "triangulation " + ofToString(triang_t) << std::endl;
    std::cout << "assembly and solving " + ofToString(total_wo_cg_t) << std::endl;
    std::cout << "b + t + a " + ofToString(bvh_t + triang_t + total_wo_cg_t) << std::endl;
    std::cout << "total_t " + ofToString(total_t) << std::endl;

    double combined_size = overlapping_meshes->mesh(0)->num_cells() + overlapping_meshes->mesh(1)->num_cells();

    timing_bvh.push_back(make_pair(combined_size, bvh_t));
    timing_triang.push_back(make_pair(combined_size,triang_t));
    timing_wo_cg.push_back(make_pair(combined_size,total_wo_cg_t));
    timing_total.push_back(make_pair(combined_size,bvh_t + triang_t + total_wo_cg_t));
    size_mesh1.push_back(make_pair(combined_size, overlapping_meshes->mesh(0)->num_cells()));
    size_mesh2.push_back(make_pair(combined_size, overlapping_meshes->mesh(1)->num_cells()));
    size_triangulation.push_back(make_pair(combined_size, interface_mesh->num_cells()));
    size_interface.push_back(make_pair(combined_size, cut_mesh_0->num_cells()));

    dolfin::list_timings(true);

    std::cout<< "size, meshes " + ofToString(overlapping_meshes->mesh(0)->num_cells()) << std::endl;
    std::cout<< "size, meshes " + ofToString(overlapping_meshes->mesh(1)->num_cells()) << std::endl;
    std::cout << "size interface_mesh " + ofToString(interface_mesh->num_cells()) << std::endl;
    std::cout << "size cut_mesh " + ofToString(cut_mesh_0->num_cells()) << std::endl;

    std::vector<dolfin::la_index> dof_vertex_set_1 = vertex_to_dof_map(*u.part(0).function_space());
    std::vector<double> values_1(dof_vertex_set_1.size());
    u.part(0).vector()->get_local(values_1.data(), dof_vertex_set_1.size(), dof_vertex_set_1.data());
    vert_values_1 = values_1;
    u_min_1 = (u.part(0).vector()->min());
    u_max_1 = (u.part(0).vector()->max());


    std::vector<dolfin::la_index> dof_vertex_set_2 = vertex_to_dof_map(*u.part(1).function_space());
    std::vector<double> values_2(dof_vertex_set_2.size());
    u.part(1).vector()->get_local(values_2.data(), dof_vertex_set_2.size(), dof_vertex_set_2.data());
    vert_values_2 = values_2;
    u_min_2 = (u.part(1).vector()->min());
    u_max_2 = (u.part(1).vector()->max());

    draw_t_0_1 = true;
    draw_first_mesh_t_0_u = true;
    draw_triangulation_mesh_u = true;
    draw_second_mesh_u = true;
    draw_first_mesh_u = false;
    draw_first_mesh = false;
    draw_second_mesh = true;
    draw_triangulation_mesh = true;
    draw_indexes = false;
    draw_index_triangulation = false;
    draw_border_first = true;
    //dolfin::File("u_olm_0.pvd") << u.part(0);
    //dolfin::File("u_olm_1.pvd") << u.part(1);
    int endsolve = ofGetElapsedTimeMillis() - startsolve;
    std::cout << "time to solve overlapping meshes Poisson " + ofToString(endsolve) + "ms"<< std::endl;
}


void ofApp::runPoissonOLM_convergence_p1()
{
    int startsolve = ofGetElapsedTimeMillis();
    // Reuse mesh 0 from MatchingCompositeMesh
    dolfin::Matrix A_2;

    // Create function spaces on each mesh and composite function space
    // Create function spaces on each mesh and composite function space
    PoissonOLM0::FunctionSpace V0(overlapping_meshes->mesh(0));
    PoissonOLM1::FunctionSpace V1(overlapping_meshes->mesh(1));

    cutfem::CompositeFunctionSpace V_c(overlapping_meshes);
    V_c.add(V0);
    V_c.add(V1);
    V_c.build();

    V_c.view(0);

    // Create empty CompositeForm
    cutfem::CompositeForm a_c(V_c, V_c);

    // Create standard bilinear forms on mesh 0
    PoissonOLM0::BilinearForm a0(V0, V0);
    dolfin::Constant gamma(gamma_u);
    a0.gamma = gamma;

    // Get domain marker for mesh 0;
    auto cell_marker_0 = overlapping_meshes->cell_marker(0);

    // Attach data to standard forms
    a0.set_cell_domains(cell_marker_0);

    // Compute constrained dofs for V0
    auto & constrained_dofs = V_c.constrained_dofs();

    // (Candidates are those ones which lives in elements marked with 2).
    cutfem::CutFEMTools::compute_constrained_dofs(constrained_dofs,
            *V_c.view(0),
            cell_marker_0.get(),
            2);

    std::size_t dof_counter = 0;
    for (std::size_t i = 0; i < constrained_dofs.size(); ++i)
        if (constrained_dofs[i])
        {
            ++dof_counter;
            //dolfin::info("Dof %d is constrained.", i);
        }
    dolfin::info("In total, %d are constrained", dof_counter);

    cutfem::CompositeFunction vis_constr_dofs(V_c);
    cutfem::CutFEMTools::visualise_constrained_dofs(vis_constr_dofs,
            constrained_dofs);

    //dolfin::File("constrained_dofs.pvd") << vis_constr_dofs.part(0);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a0)));

    // Get the cut mesh describing cut elements in mesh 0
    auto cut_mesh_and_parent_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(0);
    auto cut_mesh_0 = cut_mesh_and_parent_ids.first;
    auto parent_mesh_ids_0 = cut_mesh_and_parent_ids.second;

    //DEBUG
    //dolfin::File("cut_mesh_0.pvd") << *cut_mesh_0;

    dolfin::info("parent_mesh_ids_0 = %d", parent_mesh_ids_0.size());
    // Create quadrature rule for mesh 0
    std::size_t order = 1;
    std::shared_ptr<cutfem::Quadrature> quadrature_a0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    a_c.cut_form(0)->set_quadrature(1, quadrature_a0_domain_1);
    a_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    a_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::BilinearForm a1(V1, V1);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a1)));

    //Add interface form
    // FIXME: Passed function space is just dummy.
    PoissonOLMInterface::BilinearForm ai(V0, V0);
    ai.gamma = gamma;
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(ai)));

    // Get interface cut mesh and related information, stored as the second
    // cut_mesh_and_parent_ids component in TomOverlappingMeshes
    auto interface_mesh_and_parent_mesh_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(1);
    auto interface_mesh = interface_mesh_and_parent_mesh_ids.first;
    auto interface_parent_mesh_ids = interface_mesh_and_parent_mesh_ids.second;

    //DEBUG
    //dolfin::File("interface_mesh.pvd") << *interface_mesh;

    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_interface(
        new cutfem::Quadrature(interface_mesh->type().cell_type(),
                               interface_mesh->geometry().dim(), order));

    // Build normal field for quadrature points
    std::shared_ptr<cutfem::FacetNormals> interface_facet_normals(
        new cutfem::FacetNormals(interface_mesh->facet_normals(),
                                 interface_mesh->geometry().dim(),
                                 quadrature_interface->size()));

    a_c.cut_form(2)->set_quadrature(0, quadrature_interface);
    a_c.cut_form(2)->set_cut_mesh(0, interface_mesh);
    a_c.cut_form(2)->set_parent_mesh_ids(0, interface_parent_mesh_ids);
    a_c.cut_form(2)->set_facet_normals(0, interface_facet_normals);


    // Assemble composite form and compare.
    cutfem::CompositeFormAssembler assembler;
    assembler.assemble(A_2, a_c);

    // Temporary assemble of rhs and solution.
    cutfem::CompositeForm L_c(V_c);

    //dolfin::Constant f(1);

    using namespace dolfin;

    // Source term (right-hand side)
    class Source : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            //float k_0 = 1;
            //float k_1 = 1;
            //values[0] = DOLFIN_PI*DOLFIN_PI*(k_0*k_0*sin((DOLFIN_PI)*x[0])*sin(DOLFIN_PI*x[1]) + k_1*k_1*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
            values[0] = -6*(x[0] + x[1]);
        }
    };

    class DirichletBC : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            //values[0] = 0;
            values[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
        }
    };

    PoissonOLM0::LinearForm L0(V0);

    // Attach data to standard forms
    Source f;
    L0.f = f;

    DirichletBC g;
    L0.g = g;

    L0.gamma = gamma;

    L0.set_cell_domains(cell_marker_0);
    //L0.f = f;

    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L0)));

    // Get the cut mesh describing cut elements in mesh 0

    // Create quadrature rule for mesh 0
    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_L0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    L_c.cut_form(0)->set_quadrature(1, quadrature_L0_domain_1);
    L_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    L_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::LinearForm L1(V1);
    //L1.f = f;
    L1.f = f;
    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L1)));

    dolfin::Vector b_2;
    assembler.assemble(b_2, L_c);

    cutfem::CompositeFunction u(V_c);
    dolfin::solve(A_2, *u.vector(), b_2);


    std::vector<dolfin::la_index> dof_vertex_set_1 = vertex_to_dof_map(*u.part(0).function_space());
    std::vector<double> values_1(dof_vertex_set_1.size());
    u.part(0).vector()->get_local(values_1.data(), dof_vertex_set_1.size(), dof_vertex_set_1.data());
    vert_values_1 = values_1;
    u_min_1 = (u.part(0).vector()->min());
    u_max_1 = (u.part(0).vector()->max());


    std::vector<dolfin::la_index> dof_vertex_set_2 = vertex_to_dof_map(*u.part(1).function_space());
    std::vector<double> values_2(dof_vertex_set_2.size());
    u.part(1).vector()->get_local(values_2.data(), dof_vertex_set_2.size(), dof_vertex_set_2.data());
    vert_values_2 = values_2;
    u_min_2 = (u.part(1).vector()->min());
    u_max_2 = (u.part(1).vector()->max());


    draw_t_0_1 = false;
    draw_first_mesh_t_0_u = false;
    draw_triangulation_mesh_u = false;
    draw_second_mesh_u = false;
    draw_first_mesh_u = false;
    draw_first_mesh = false;
    draw_second_mesh = false;
    draw_triangulation_mesh = false;
    draw_indexes = false;
    draw_index_triangulation = false;
    draw_border_first = false;
    //dolfin::File("u_olm_0.pvd") << u.part(0);
    //dolfin::File("u_olm_1.pvd") << u.part(1);
    int endsolve = ofGetElapsedTimeMillis() - startsolve;
    std::cout << "time to solve overlapping meshes Poisson " + ofToString(endsolve) + "ms"<< std::endl;

    using namespace dolfin;

    // Source term (right-hand side)
    class Solution : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            //float k_0 = 1;
            //float k_1 = 1;
            //values[0] = (sin(k_0*DOLFIN_PI*x[0])*sin(k_1*DOLFIN_PI*x[1]));
             values[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
        }
    };
    //Reference solution and its interpolation


    Solution u_sol;

    //H1_Norm_u::FunctionSpace P3(overlapping_meshes->mesh(0));

    // mesh
    P3::FunctionSpace P3_0(overlapping_meshes->mesh(0));
    dolfin::Function u_sol_inter_0(P3_0);
    u_sol_inter_0.interpolate(u_sol);
    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_0(P3_0);
    u_error_0.interpolate(u.part(0));

    *u_error_0.vector() -= *u_sol_inter_0.vector();
    H1_Norm_u_0_only::Functional u_h1_norm_0(overlapping_meshes->mesh(0),u_error_0);
    u_h1_norm_0.set_cell_domains(overlapping_meshes->cell_marker(0));
    double H1_Norm_u_0 = std::sqrt(assemble(u_h1_norm_0));



    // collidingMesh
    P3::FunctionSpace P3_1(overlapping_meshes->mesh(1));
    dolfin::Function u_sol_inter_1(P3_1);
    u_sol_inter_1.interpolate(u_sol);

    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_1(P3_1);
    u_error_1.interpolate(u.part(1));

    *u_error_1.vector() -= *u_sol_inter_1.vector();
    H1_Norm_u_1::Functional u_h1_norm_1(overlapping_meshes->mesh(1),u_error_1);
    double H1_Norm_u_1 = std::sqrt(assemble(u_h1_norm_1));


    // Print
    std::cout << "hmax t_0 + t_gamma : " + ofToString(overlapping_meshes->mesh(0).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_0) << std::endl;
    std::cout << "hmax t_2: " + ofToString(overlapping_meshes->mesh(1).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_1) << std::endl;
    std::cout << "HMAX COMBINED : " + ofToString(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2))) + " H1 norm: " + ofToString(std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2)))) << std::endl;

    std::pair<double, double> t_0 = std::make_pair(overlapping_meshes->mesh(0).get()->hmax(), H1_Norm_u_0);
    std::pair<double, double> t_2 = std::make_pair(overlapping_meshes->mesh(1).get()->hmax(), H1_Norm_u_1);
    std::pair<double, double> t_combined = std::make_pair(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2)),std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2))));

    t_0_strings.push_back(t_0);
    t_2_strings.push_back(t_2);
    combined_strings.push_back(t_combined);
}

void ofApp::runPoissonOLM_convergence_p2()
{
    int startsolve = ofGetElapsedTimeMillis();
    // Reuse mesh 0 from MatchingCompositeMesh
    dolfin::Matrix A_2;

    // Create function spaces on each mesh and composite function space
    // Create function spaces on each mesh and composite function space
    PoissonOLM0_p2::FunctionSpace V0(overlapping_meshes->mesh(0));
    PoissonOLM1_p2::FunctionSpace V1(overlapping_meshes->mesh(1));

    cutfem::CompositeFunctionSpace V_c(overlapping_meshes);
    V_c.add(V0);
    V_c.add(V1);
    V_c.build();

    V_c.view(0);

    // Create empty CompositeForm
    cutfem::CompositeForm a_c(V_c, V_c);

    // Create standard bilinear forms on mesh 0
    PoissonOLM0_p2::BilinearForm a0(V0, V0);
    dolfin::Constant gamma(gamma_u);
    a0.gamma = gamma;
    // Get domain marker for mesh 0;
    auto cell_marker_0 = overlapping_meshes->cell_marker(0);

    // Attach data to standard forms
    a0.set_cell_domains(cell_marker_0);

    // Compute constrained dofs for V0
    auto & constrained_dofs = V_c.constrained_dofs();

    // (Candidates are those ones which lives in elements marked with 2).
    cutfem::CutFEMTools::compute_constrained_dofs(constrained_dofs,
            *V_c.view(0),
            cell_marker_0.get(),
            2);

    std::size_t dof_counter = 0;
    for (std::size_t i = 0; i < constrained_dofs.size(); ++i)
        if (constrained_dofs[i])
        {
            ++dof_counter;
            //dolfin::info("Dof %d is constrained.", i);
        }

    cutfem::CompositeFunction vis_constr_dofs(V_c);
    cutfem::CutFEMTools::visualise_constrained_dofs(vis_constr_dofs,
            constrained_dofs);

    dolfin::File("constrained_dofs.pvd") << vis_constr_dofs.part(0);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a0)));

    // Get the cut mesh describing cut elements in mesh 0
    auto cut_mesh_and_parent_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(0);
    auto cut_mesh_0 = cut_mesh_and_parent_ids.first;
    auto parent_mesh_ids_0 = cut_mesh_and_parent_ids.second;

    //DEBUG
    dolfin::File("cut_mesh_0.pvd") << *cut_mesh_0;

    dolfin::info("parent_mesh_ids_0 = %d", parent_mesh_ids_0.size());
    // Create quadrature rule for mesh 0
    std::size_t order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_a0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0_p2 has domain_id 1!
    a_c.cut_form(0)->set_quadrature(1, quadrature_a0_domain_1);
    a_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    a_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);


    // Create standard bilinear forms on mesh 1
    PoissonOLM1_p2::BilinearForm a1(V1, V1);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a1)));


    //Add interface form
    // FIXME: Passed function space is just dummy.
    PoissonOLMInterface_p2::BilinearForm ai(V0, V0);
    ai.gamma = gamma;
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(ai)));

    // Get interface cut mesh and related information, stored as the second
    // cut_mesh_and_parent_ids component in TomOverlappingMeshes
    auto interface_mesh_and_parent_mesh_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(1);
    auto interface_mesh = interface_mesh_and_parent_mesh_ids.first;
    auto interface_parent_mesh_ids = interface_mesh_and_parent_mesh_ids.second;

    order = 4;
    std::shared_ptr<cutfem::Quadrature> quadrature_interface(
        new cutfem::Quadrature(interface_mesh->type().cell_type(),
                               interface_mesh->geometry().dim(), order));

    // Build normal field for quadrature points
    std::shared_ptr<cutfem::FacetNormals> interface_facet_normals(
        new cutfem::FacetNormals(interface_mesh->facet_normals(),
                                 interface_mesh->geometry().dim(),
                                 quadrature_interface->size()));

    a_c.cut_form(2)->set_quadrature(0, quadrature_interface);
    a_c.cut_form(2)->set_cut_mesh(0, interface_mesh);
    a_c.cut_form(2)->set_parent_mesh_ids(0, interface_parent_mesh_ids);
    a_c.cut_form(2)->set_facet_normals(0, interface_facet_normals);

    // Assemble composite form and compare.
    cutfem::CompositeFormAssembler assembler;
    assembler.assemble(A_2, a_c);

    // Temporary assemble of rhs and solution.
    cutfem::CompositeForm L_c(V_c);

    //dolfin::Constant f(1);

    using namespace dolfin;

    // Source term (right-hand side)
    class Source : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            //values[0] = 2*DOLFIN_PI*DOLFIN_PI*(sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]));
            values[0] = -6*(x[0] + x[1]);
        }
    };

    class DirichletBC : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
           values[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
           //values[0] = 0;
        }
    };

    PoissonOLM0_p2::LinearForm L0(V0);

    // Attach data to standard forms

    Source f;
    L0.f = f;
    DirichletBC g;
    L0.g = g;
    L0.gamma = gamma;

    L0.set_cell_domains(cell_marker_0);

    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L0)));

    // Get the cut mesh describing cut elements in mesh 0

    // Create quadrature rule for mesh 0
    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_L0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0_p2 has domain_id 1!
    L_c.cut_form(0)->set_quadrature(1, quadrature_L0_domain_1);
    L_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    L_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1_p2::LinearForm L1(V1);
    //L1.f = f;
    L1.f = f;
    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L1)));

    dolfin::Vector b_2;
    assembler.assemble(b_2, L_c);

    cutfem::CompositeFunction u(V_c);
    dolfin::solve(A_2, *u.vector(), b_2);
    /*
        std::vector<dolfin::la_index> dof_vertex_set_1 = vertex_to_dof_map(*u.part(0).function_space());
        std::vector<double> values_1(dof_vertex_set_1.size());
        u.part(0).vector()->get_local(values_1.data(), dof_vertex_set_1.size(), dof_vertex_set_1.data());
        vert_values_1 = values_1;
        u_min_1 = (u.part(0).vector()->min());
        u_max_1 = (u.part(0).vector()->max());
    std::cout << "gets here14"<< std::endl;

        std::vector<dolfin::la_index> dof_vertex_set_2 = vertex_to_dof_map(*u.part(1).function_space());
        std::vector<double> values_2(dof_vertex_set_2.size());
        u.part(1).vector()->get_local(values_2.data(), dof_vertex_set_2.size(), dof_vertex_set_2.data());
        vert_values_2 = values_2;
        u_min_2 = (u.part(1).vector()->min());
        u_max_2 = (u.part(1).vector()->max());
        std::cout << "gets here15"<< std::endl;


        draw_t_0_1 = false;
        draw_first_mesh_t_0_u = false;
        draw_triangulation_mesh_u = false;
        draw_second_mesh_u = false;
        draw_first_mesh_u = false;
        draw_first_mesh = false;
        draw_second_mesh = false;
        draw_triangulation_mesh = false;
        draw_indexes = false;
        draw_index_triangulation = false;
        draw_border_first = false;
        //dolfin::File("u_olm_0.pvd") << u.part(0);
        //dolfin::File("u_olm_1.pvd") << u.part(1);
        int endsolve = ofGetElapsedTimeMillis() - startsolve;
        std::cout << "time to solve overlapping meshes Poisson " + ofToString(endsolve) + "ms"<< std::endl;
    */

    using namespace dolfin;

    // Source term (right-hand side)
    class Solution : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            //values[0] = sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]);
            values[0] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
        }
    };
    //Reference solution and its interpolation


    Solution u_sol;

    // mesh
    P3::FunctionSpace P3_0(overlapping_meshes->mesh(0));
    dolfin::Function u_sol_inter_0(P3_0);
    u_sol_inter_0.interpolate(u_sol);
    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_0(P3_0);
    u_error_0.interpolate(u.part(0));

    *u_error_0.vector() -= *u_sol_inter_0.vector();
    H1_Norm_u_0::Functional u_h1_norm_0(overlapping_meshes->mesh(0),u_error_0);
    u_h1_norm_0.set_cell_domains(overlapping_meshes->cell_marker(0));
    double H1_Norm_u_0 = std::sqrt(assemble(u_h1_norm_0));

    // collidingMesh
    P3::FunctionSpace P3_1(overlapping_meshes->mesh(1));
    dolfin::Function u_sol_inter_1(P3_1);
    u_sol_inter_1.interpolate(u_sol);

    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_1(P3_1);
    u_error_1.interpolate(u.part(1));

    *u_error_1.vector() -= *u_sol_inter_1.vector();
    H1_Norm_u_1::Functional u_h1_norm_1(overlapping_meshes->mesh(1),u_error_1);
    double H1_Norm_u_1 = std::sqrt(assemble(u_h1_norm_1));

    dolfin::File("u_olm_0.pvd") << u.part(0);
    dolfin::File("u_olm_1.pvd") << u.part(1);

    // Print
    std::cout << "hmax t_0 + t_gamma : " + ofToString(overlapping_meshes->mesh(0).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_0) << std::endl;
    std::cout << "hmax t_2: " + ofToString(overlapping_meshes->mesh(1).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_1) << std::endl;
    std::cout << "HMAX COMBINED : " + ofToString(max(overlapping_meshes->mesh(1).get()->hmax(), overlapping_meshes->mesh(0).get()->hmax())) + " H1 norm: " + ofToString(std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2)))) << std::endl;

    std::pair<double, double> t_0 = std::make_pair(overlapping_meshes->mesh(0).get()->hmax(), H1_Norm_u_0);
    std::pair<double, double> t_2 = std::make_pair(overlapping_meshes->mesh(1).get()->hmax(), H1_Norm_u_1);
    std::pair<double, double> t_combined = std::make_pair(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2)),std::sqrt( ( pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2) ) ) );

    t_0_strings_p2.push_back(t_0);
    t_2_strings_p2.push_back(t_2);
    combined_strings_p2.push_back(t_combined);
}

void ofApp::runPoissonOLM_patch_test_p1()
{
    int startsolve = ofGetElapsedTimeMillis();
    // Reuse mesh 0 from MatchingCompositeMesh
    dolfin::Matrix A_2;

    // Create function spaces on each mesh and composite function space
    // Create function spaces on each mesh and composite function space
    PoissonOLM0::FunctionSpace V0(overlapping_meshes->mesh(0));
    PoissonOLM1::FunctionSpace V1(overlapping_meshes->mesh(1));

    cutfem::CompositeFunctionSpace V_c(overlapping_meshes);
    V_c.add(V0);
    V_c.add(V1);
    V_c.build();

    V_c.view(0);

    // Create empty CompositeForm
    cutfem::CompositeForm a_c(V_c, V_c);

    // Create standard bilinear forms on mesh 0
    PoissonOLM0::BilinearForm a0(V0, V0);
    dolfin::Constant gamma(gamma_u);
    a0.gamma = gamma;

    // Get domain marker for mesh 0;
    auto cell_marker_0 = overlapping_meshes->cell_marker(0);

    // Attach data to standard forms
    a0.set_cell_domains(cell_marker_0);

    // Compute constrained dofs for V0
    auto & constrained_dofs = V_c.constrained_dofs();

    // (Candidates are those ones which lives in elements marked with 2).
    cutfem::CutFEMTools::compute_constrained_dofs(constrained_dofs,
            *V_c.view(0),
            cell_marker_0.get(),
            2);

    std::size_t dof_counter = 0;
    for (std::size_t i = 0; i < constrained_dofs.size(); ++i)
        if (constrained_dofs[i])
        {
            ++dof_counter;
            //dolfin::info("Dof %d is constrained.", i);
        }
    dolfin::info("In total, %d are constrained", dof_counter);

    cutfem::CompositeFunction vis_constr_dofs(V_c);
    cutfem::CutFEMTools::visualise_constrained_dofs(vis_constr_dofs,
            constrained_dofs);

    //dolfin::File("constrained_dofs.pvd") << vis_constr_dofs.part(0);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a0)));

    // Get the cut mesh describing cut elements in mesh 0
    auto cut_mesh_and_parent_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(0);
    auto cut_mesh_0 = cut_mesh_and_parent_ids.first;
    auto parent_mesh_ids_0 = cut_mesh_and_parent_ids.second;

    //DEBUG
    //dolfin::File("cut_mesh_0.pvd") << *cut_mesh_0;

    dolfin::info("parent_mesh_ids_0 = %d", parent_mesh_ids_0.size());
    // Create quadrature rule for mesh 0
    std::size_t order = 1;
    std::shared_ptr<cutfem::Quadrature> quadrature_a0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    a_c.cut_form(0)->set_quadrature(1, quadrature_a0_domain_1);
    a_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    a_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::BilinearForm a1(V1, V1);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a1)));

    //Add interface form
    // FIXME: Passed function space is just dummy.
    PoissonOLMInterface::BilinearForm ai(V0, V0);
    ai.gamma = gamma;
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(ai)));

    // Get interface cut mesh and related information, stored as the second
    // cut_mesh_and_parent_ids component in TomOverlappingMeshes
    auto interface_mesh_and_parent_mesh_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(1);
    auto interface_mesh = interface_mesh_and_parent_mesh_ids.first;
    auto interface_parent_mesh_ids = interface_mesh_and_parent_mesh_ids.second;

    //DEBUG
    //dolfin::File("interface_mesh.pvd") << *interface_mesh;

    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_interface(
        new cutfem::Quadrature(interface_mesh->type().cell_type(),
                               interface_mesh->geometry().dim(), order));

    // Build normal field for quadrature points
    std::shared_ptr<cutfem::FacetNormals> interface_facet_normals(
        new cutfem::FacetNormals(interface_mesh->facet_normals(),
                                 interface_mesh->geometry().dim(),
                                 quadrature_interface->size()));

    a_c.cut_form(2)->set_quadrature(0, quadrature_interface);
    a_c.cut_form(2)->set_cut_mesh(0, interface_mesh);
    a_c.cut_form(2)->set_parent_mesh_ids(0, interface_parent_mesh_ids);
    a_c.cut_form(2)->set_facet_normals(0, interface_facet_normals);


    // Assemble composite form and compare.
    cutfem::CompositeFormAssembler assembler;
    assembler.assemble(A_2, a_c);

    // Temporary assemble of rhs and solution.
    cutfem::CompositeForm L_c(V_c);

    //dolfin::Constant f(1);

    using namespace dolfin;

    // Source term (right-hand side)
    class Source : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = Constant(0);
        }
    };

    class DirichletBC : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = 1 + x[0] + x[1];
        }
    };

    PoissonOLM0::LinearForm L0(V0);

    // Attach data to standard forms

    Source f;
    L0.f = f;
    DirichletBC g;
    L0.g = g;
    L0.gamma = gamma;

    L0.set_cell_domains(cell_marker_0);
    //L0.f = f;

    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L0)));

    // Get the cut mesh describing cut elements in mesh 0

    // Create quadrature rule for mesh 0
    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_L0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0 has domain_id 1!
    L_c.cut_form(0)->set_quadrature(1, quadrature_L0_domain_1);
    L_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    L_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1::LinearForm L1(V1);
    //L1.f = f;
    L1.f = f;
    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L1)));

    dolfin::Vector b_2;
    assembler.assemble(b_2, L_c);

    cutfem::CompositeFunction u(V_c);
    dolfin::solve(A_2, *u.vector(), b_2);


    std::vector<dolfin::la_index> dof_vertex_set_1 = vertex_to_dof_map(*u.part(0).function_space());
    std::vector<double> values_1(dof_vertex_set_1.size());
    u.part(0).vector()->get_local(values_1.data(), dof_vertex_set_1.size(), dof_vertex_set_1.data());
    vert_values_1 = values_1;
    u_min_1 = (u.part(0).vector()->min());
    u_max_1 = (u.part(0).vector()->max());


    std::vector<dolfin::la_index> dof_vertex_set_2 = vertex_to_dof_map(*u.part(1).function_space());
    std::vector<double> values_2(dof_vertex_set_2.size());
    u.part(1).vector()->get_local(values_2.data(), dof_vertex_set_2.size(), dof_vertex_set_2.data());
    vert_values_2 = values_2;
    u_min_2 = (u.part(1).vector()->min());
    u_max_2 = (u.part(1).vector()->max());


    draw_t_0_1 = false;
    draw_first_mesh_t_0_u = false;
    draw_triangulation_mesh_u = false;
    draw_second_mesh_u = false;
    draw_first_mesh_u = false;
    draw_first_mesh = false;
    draw_second_mesh = false;
    draw_triangulation_mesh = false;
    draw_indexes = false;
    draw_index_triangulation = false;
    draw_border_first = false;
    //dolfin::File("u_olm_0.pvd") << u.part(0);
    //dolfin::File("u_olm_1.pvd") << u.part(1);
    int endsolve = ofGetElapsedTimeMillis() - startsolve;
    std::cout << "time to solve overlapping meshes Poisson " + ofToString(endsolve) + "ms"<< std::endl;

    using namespace dolfin;

    // Source term (right-hand side)
    class Solution : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = 1 + x[0] + x[1];
        }
    };
    //Reference solution and its interpolation


    Solution u_sol;

    //H1_Norm_u::FunctionSpace P3(overlapping_meshes->mesh(0));

    // mesh
    P3::FunctionSpace P3_0(overlapping_meshes->mesh(0));
    dolfin::Function u_sol_inter_0(P3_0);
    u_sol_inter_0.interpolate(u_sol);
    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_0(P3_0);
    u_error_0.interpolate(u.part(0));

    *u_error_0.vector() -= *u_sol_inter_0.vector();
    H1_Norm_u_0::Functional u_h1_norm_0(overlapping_meshes->mesh(0),u_error_0);
    u_h1_norm_0.set_cell_domains(overlapping_meshes->cell_marker(0));
    double H1_Norm_u_0 = std::sqrt(assemble(u_h1_norm_0));



    // collidingMesh
    P3::FunctionSpace P3_1(overlapping_meshes->mesh(1));
    dolfin::Function u_sol_inter_1(P3_1);
    u_sol_inter_1.interpolate(u_sol);

    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_1(P3_1);
    u_error_1.interpolate(u.part(1));

    *u_error_1.vector() -= *u_sol_inter_1.vector();
    H1_Norm_u_1::Functional u_h1_norm_1(overlapping_meshes->mesh(1),u_error_1);
    double H1_Norm_u_1 = std::sqrt(assemble(u_h1_norm_1));


    // Print
    std::cout << "hmax t_0 + t_gamma : " + ofToString(overlapping_meshes->mesh(0).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_0) << std::endl;
    std::cout << "hmax t_2: " + ofToString(overlapping_meshes->mesh(1).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_1) << std::endl;
    std::cout << "HMAX COMBINED : " + ofToString(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2))) + " H1 norm: " + ofToString(std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2)))) << std::endl;

    std::pair<double, double> t_0 = std::make_pair(overlapping_meshes->mesh(0).get()->hmax(), H1_Norm_u_0);
    std::pair<double, double> t_2 = std::make_pair(overlapping_meshes->mesh(1).get()->hmax(), H1_Norm_u_1);
    std::pair<double, double> t_combined = std::make_pair(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2)),std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2))));
    t_0_strings.push_back(t_0);
    t_2_strings.push_back(t_2);
    combined_strings.push_back(t_combined);

}

void ofApp::runPoissonOLM_patch_test_p2()
{
    int startsolve = ofGetElapsedTimeMillis();
    // Reuse mesh 0 from MatchingCompositeMesh
    dolfin::Matrix A_2;

    // Create function spaces on each mesh and composite function space
    // Create function spaces on each mesh and composite function space
    PoissonOLM0_p2::FunctionSpace V0(overlapping_meshes->mesh(0));
    PoissonOLM1_p2::FunctionSpace V1(overlapping_meshes->mesh(1));

    cutfem::CompositeFunctionSpace V_c(overlapping_meshes);
    V_c.add(V0);
    V_c.add(V1);
    V_c.build();

    V_c.view(0);

    // Create empty CompositeForm
    cutfem::CompositeForm a_c(V_c, V_c);

    // Create standard bilinear forms on mesh 0
    PoissonOLM0_p2::BilinearForm a0(V0, V0);
    dolfin::Constant gamma(gamma_u);
    a0.gamma = gamma;
    // Get domain marker for mesh 0;
    auto cell_marker_0 = overlapping_meshes->cell_marker(0);

    // Attach data to standard forms
    a0.set_cell_domains(cell_marker_0);

    // Compute constrained dofs for V0
    auto & constrained_dofs = V_c.constrained_dofs();

    // (Candidates are those ones which lives in elements marked with 2).
    cutfem::CutFEMTools::compute_constrained_dofs(constrained_dofs,
            *V_c.view(0),
            cell_marker_0.get(),
            2);

    std::size_t dof_counter = 0;
    for (std::size_t i = 0; i < constrained_dofs.size(); ++i)
        if (constrained_dofs[i])
        {
            ++dof_counter;
            //dolfin::info("Dof %d is constrained.", i);
        }
    cutfem::CompositeFunction vis_constr_dofs(V_c);
    cutfem::CutFEMTools::visualise_constrained_dofs(vis_constr_dofs,
            constrained_dofs);

    //dolfin::File("constrained_dofs.pvd") << vis_constr_dofs.part(0);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a0)));

    // Get the cut mesh describing cut elements in mesh 0
    auto cut_mesh_and_parent_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(0);
    auto cut_mesh_0 = cut_mesh_and_parent_ids.first;
    auto parent_mesh_ids_0 = cut_mesh_and_parent_ids.second;

    //DEBUG
    //dolfin::File("cut_mesh_0.pvd") << *cut_mesh_0;

    dolfin::info("parent_mesh_ids_0 = %d", parent_mesh_ids_0.size());
    // Create quadrature rule for mesh 0
    std::size_t order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_a0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0_p2 has domain_id 1!
    a_c.cut_form(0)->set_quadrature(1, quadrature_a0_domain_1);
    a_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    a_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);


    // Create standard bilinear forms on mesh 1
    PoissonOLM1_p2::BilinearForm a1(V1, V1);

    // Add them as CutForm to CompositeForm
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(a1)));


    //Add interface form
    // FIXME: Passed function space is just dummy.
    PoissonOLMInterface_p2::BilinearForm ai(V0, V0);
    ai.gamma = gamma;
    a_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(ai)));

    // Get interface cut mesh and related information, stored as the second
    // cut_mesh_and_parent_ids component in TomOverlappingMeshes
    auto interface_mesh_and_parent_mesh_ids =
        overlapping_meshes->cut_mesh_and_parent_mesh_ids(1);
    auto interface_mesh = interface_mesh_and_parent_mesh_ids.first;
    auto interface_parent_mesh_ids = interface_mesh_and_parent_mesh_ids.second;

    order = 4;
    std::shared_ptr<cutfem::Quadrature> quadrature_interface(
        new cutfem::Quadrature(interface_mesh->type().cell_type(),
                               interface_mesh->geometry().dim(), order));

    // Build normal field for quadrature points
    std::shared_ptr<cutfem::FacetNormals> interface_facet_normals(
        new cutfem::FacetNormals(interface_mesh->facet_normals(),
                                 interface_mesh->geometry().dim(),
                                 quadrature_interface->size()));

    a_c.cut_form(2)->set_quadrature(0, quadrature_interface);
    a_c.cut_form(2)->set_cut_mesh(0, interface_mesh);
    a_c.cut_form(2)->set_parent_mesh_ids(0, interface_parent_mesh_ids);
    a_c.cut_form(2)->set_facet_normals(0, interface_facet_normals);

    // Assemble composite form and compare.
    cutfem::CompositeFormAssembler assembler;
    assembler.assemble(A_2, a_c);

    // Temporary assemble of rhs and solution.
    cutfem::CompositeForm L_c(V_c);

    //dolfin::Constant f(1);

    using namespace dolfin;

    // Source term (right-hand side)
    class Source : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = -4;
        }
    };

    class BoundaryCondition : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = x[0]*x[0] + x[1]*x[1] + x[0]*x[1] + x[0] + x[1] + 1;
        }
    };

    PoissonOLM0_p2::LinearForm L0(V0);

    // Attach data to standard forms

    Source f;
    L0.f = f;
    BoundaryCondition g;
    L0.g = g;
    L0.gamma = gamma;

    L0.set_cell_domains(cell_marker_0);

    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L0)));

    // Get the cut mesh describing cut elements in mesh 0

    // Create quadrature rule for mesh 0
    order = 2;
    std::shared_ptr<cutfem::Quadrature> quadrature_L0_domain_1(
        new cutfem::Quadrature(cut_mesh_0->type().cell_type(),
                               cut_mesh_0->geometry().dim(), order));

    // Remember that dxq in PoissonOLM0_p2 has domain_id 1!
    L_c.cut_form(0)->set_quadrature(1, quadrature_L0_domain_1);
    L_c.cut_form(0)->set_cut_mesh(1, cut_mesh_0);
    L_c.cut_form(0)->set_parent_mesh_ids(1, parent_mesh_ids_0);

    // Create standard bilinear forms on mesh 1
    PoissonOLM1_p2::LinearForm L1(V1);
    //L1.f = f;
    L1.f = f;
    // Add them as CutForm to CompositeForm
    L_c.add(std::shared_ptr<cutfem::CutForm>(new cutfem::CutForm(L1)));

    dolfin::Vector b_2;
    assembler.assemble(b_2, L_c);

    cutfem::CompositeFunction u(V_c);
    dolfin::solve(A_2, *u.vector(), b_2);
    /*
        std::vector<dolfin::la_index> dof_vertex_set_1 = vertex_to_dof_map(*u.part(0).function_space());
        std::vector<double> values_1(dof_vertex_set_1.size());
        u.part(0).vector()->get_local(values_1.data(), dof_vertex_set_1.size(), dof_vertex_set_1.data());
        vert_values_1 = values_1;
        u_min_1 = (u.part(0).vector()->min());
        u_max_1 = (u.part(0).vector()->max());
    std::cout << "gets here14"<< std::endl;

        std::vector<dolfin::la_index> dof_vertex_set_2 = vertex_to_dof_map(*u.part(1).function_space());
        std::vector<double> values_2(dof_vertex_set_2.size());
        u.part(1).vector()->get_local(values_2.data(), dof_vertex_set_2.size(), dof_vertex_set_2.data());
        vert_values_2 = values_2;
        u_min_2 = (u.part(1).vector()->min());
        u_max_2 = (u.part(1).vector()->max());
        std::cout << "gets here15"<< std::endl;


        draw_t_0_1 = false;
        draw_first_mesh_t_0_u = false;
        draw_triangulation_mesh_u = false;
        draw_second_mesh_u = false;
        draw_first_mesh_u = false;
        draw_first_mesh = false;
        draw_second_mesh = false;
        draw_triangulation_mesh = false;
        draw_indexes = false;
        draw_index_triangulation = false;
        draw_border_first = false;
        //dolfin::File("u_olm_0.pvd") << u.part(0);
        //dolfin::File("u_olm_1.pvd") << u.part(1);
        int endsolve = ofGetElapsedTimeMillis() - startsolve;
        std::cout << "time to solve overlapping meshes Poisson " + ofToString(endsolve) + "ms"<< std::endl;
    */
    using namespace dolfin;

    // Source term (right-hand side)
    class Solution : public Expression
    {
        void eval(Array<double>& values, const Array<double>& x) const
        {
            values[0] = x[0]*x[0] + x[1]*x[1] + x[0]*x[1] + x[0] + x[1] + 1;
        }
    };
    //Reference solution and its interpolation


    Solution u_sol;

    //H1_Norm_u::FunctionSpace P3(overlapping_meshes->mesh(0));

    // mesh
    P3::FunctionSpace P3_0(overlapping_meshes->mesh(0));
    dolfin::Function u_sol_inter_0(P3_0);
    u_sol_inter_0.interpolate(u_sol);
    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_0(P3_0);
    u_error_0.interpolate(u.part(0));

    *u_error_0.vector() -= *u_sol_inter_0.vector();
    H1_Norm_u_0::Functional u_h1_norm_0(overlapping_meshes->mesh(0),u_error_0);
    u_h1_norm_0.set_cell_domains(overlapping_meshes->cell_marker(0));
    double H1_Norm_u_0 = std::sqrt(assemble(u_h1_norm_0));



    // collidingMesh
    P3::FunctionSpace P3_1(overlapping_meshes->mesh(1));
    dolfin::Function u_sol_inter_1(P3_1);
    u_sol_inter_1.interpolate(u_sol);

    // Compute error on the background mesh cells mark with either 0 or 1 FEM
    // solution "extension" and analytical solution should match on that part
    Function u_error_1(P3_1);
    u_error_1.interpolate(u.part(1));

    *u_error_1.vector() -= *u_sol_inter_1.vector();
    H1_Norm_u_1::Functional u_h1_norm_1(overlapping_meshes->mesh(1),u_error_1);
    double H1_Norm_u_1 = std::sqrt(assemble(u_h1_norm_1));


    // Print
    std::cout << "hmax t_0 + t_gamma : " + ofToString(overlapping_meshes->mesh(0).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_0) << std::endl;
    std::cout << "hmax t_2: " + ofToString(overlapping_meshes->mesh(1).get()->hmax()) + " H1 norm: " + ofToString(H1_Norm_u_1) << std::endl;
    std::cout << "HMAX COMBINED : " + ofToString(std::sqrt(pow(overlapping_meshes->mesh(1).get()->hmax(),2) + pow(overlapping_meshes->mesh(0).get()->hmax(), 2))) + " H1 norm: " + ofToString(std::sqrt((pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2)))) << std::endl;


    std::pair<double, double> t_0 = std::make_pair(overlapping_meshes->mesh(0).get()->hmax(), H1_Norm_u_0);
    std::pair<double, double> t_2 = std::make_pair(overlapping_meshes->mesh(1).get()->hmax(), H1_Norm_u_1);
    std::pair<double, double> t_combined = std::make_pair(max(overlapping_meshes->mesh(1).get()->hmax(), overlapping_meshes->mesh(0).get()->hmax()),std::sqrt( ( pow(H1_Norm_u_1,2) + pow(H1_Norm_u_0,2) ) ) );
    t_0_strings_p2.push_back(t_0);
    t_2_strings_p2.push_back(t_2);
    combined_strings_p2.push_back(t_combined);

}


/*
// Source term (right-hand side)
class Source : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double dx = x[0] - 0.5;
        double dy = x[1] - 0.5;
        values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
    }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0] = sin(5*x[0]);
    }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS;

    }
};

// Create mesh and function space
Poisson::FunctionSpace V(*mesh);

// Define boundary condition
Constant u0(0.0);
DirichletBoundary boundary;
DirichletBC bc(V, u0, boundary);

// Define variational forms
Poisson::BilinearForm a(V, V);
Poisson::LinearForm L(V);

L.f = f;
//L.g = gs;

// Compute solution
Function u(V);
solve(a == L, u, bc);



Function u(V0);
solve(a0 == L0, u, bc1);


// Assemble system
Matrix A;
Vector b;

OLMAssembler::assemble(A,a_olm);
OLMAssembler::assemble(b,L_olm);

// Apply strong bcs.
bc1.apply(A,b);

// CC: The last operation works only for bcs on mesh0. For mesh1
// the shift in the dofs is not taken into account. Use OLMAssembler::assemble_system
// instead.

//Solve linear algebra system
Vector x(A.size(1));
solve(A,x,b);

// Split up to solution to original function spaces
Function u0(V0);
Function u1(V1);
a_olm.distribute_solution(x,u0,u1);

std::cout << "gets here 12" << std::endl;
std::vector<dolfin::la_index> dof_vertex_set = vertex_to_dof_map(*u.function_space());
std::vector<double> values2(dof_vertex_set.size());
u.vector()->get_local(values2.data(), dof_vertex_set.size(), dof_vertex_set.data());
vert_values = values2;
u_min = (u.vector()->min());
u_max= (u.vector()->max());


}

}

*/

bool ofApp::point_on_line(dolfin::Point linePointA, dolfin::Point linePointB, dolfin::Point point)
{
    const double EPSILON = 0.000001f;

    double a = (linePointB.y() - linePointA.y()) / (linePointB.x() - linePointA.x());

    double b = linePointA.y() - a * linePointA.x();

    double maxx = std::max(linePointA.x(), linePointB.x());
    double minx = std::min(linePointA.x(), linePointB.x());
    double maxy = std::max(linePointA.y(), linePointB.y());
    double miny = std::min(linePointA.y(), linePointB.y());

    if ( abs(point.y() - (a*point.x()+b)) < EPSILON)
    {
        if(point.y() <= maxy + EPSILON & point.y() >= miny - EPSILON & point.x() <= maxx + EPSILON & point.x() >= minx - EPSILON)
        {
            return true;
        }
    }
    return false;
}

double ofApp::interpolate_triangle(dolfin::Point p0,
                                   double w1,
                                   dolfin::Point p1,
                                   double w2,
                                   dolfin::Point p2,
                                   double w3,
                                   dolfin::Point pi)
{

    double x1 = p0.x();
    double y1 = p0.y();
    double x2 = p1.x();
    double y2 = p1.y();
    double x3 = p2.x();
    double y3 = p2.y();
    double x = pi.x();
    double y = pi.y();


    double DET= x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3;
    double A 	= ((y2-y3)*w1+(y3-y1)*w2+(y1-y2)*w3) / DET;
    double B 	= ((x3-x2)*w1+(x1-x3)*w2+(x2-x1)*w3) / DET;
    double C 	= ((x2*y3-x3*y2)*w1+(x3*y1-x1*y3)*w2+(x1*y2-x2*y1)*w3) / DET;
    double w 	= A*x+B*y+C;

    return w;
}

bool ofApp::do_intersect_two_lines(dolfin::Point p1, dolfin::Point p2, dolfin::Point p3, dolfin::Point p4)
{
// Store the values for fast access and easy
// equations-to-code conversion
    double x1 = p1.x(), x2 = p2.x(), x3 = p3.x(), x4 = p4.x();
    double y1 = p1.y(), y2 = p2.y(), y3 = p3.y(), y4 = p4.y();

    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
// If d is zero, there is no intersection
    if (d == 0) return false;

// Get the x and y
    double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
    double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;

    // Check if the x and y coordinates are within both lines
    if ( x < std::min(x1, x2) - 3.0e-13 || x > std::max(x1, x2) + 3.0e-13 ||
            x < std::min(x3, x4) - 3.0e-13 || x > std::max(x3, x4) + 3.0e-13 )
    {
        return false;
    }

    if ( y < std::min(y1, y2) - 3.0e-13 || y > std::max(y1, y2) + 3.0e-13||
            y < std::min(y3, y4) - 3.0e-13 || y > std::max(y3, y4) + 3.0e-13)
    {
        return false;
    }
    // Return the point of intersection
    return true;
}

void ofApp::runPoissonFiDo()
{
    /*
        using namespace cutfem;
        using namespace dolfin;

        class LevelSet : public Expression
        {
        public:
            void eval(Array<double>& values, const Array<double>& x) const
            {
                double f1 = (x[0] - 0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5) - 0.1;
                values[0] = f1;

            }
        };

        set_log_level(DBG);
        size_t N = 50;
        //boost::shared_ptr<Mesh> mesh(new RectangleMesh(-1.5,-1.5,1.5,1.5,N,N));
        std::shared_ptr<GenericFunction> level_set(new LevelSet);

        std::shared_ptr<Mesh> poisson_mesh(new Mesh(*mesh));
        //auto mesh_poisson = boost::make_shared<const Mesh>(mesh);
        MeshCutter mesh_cutter(poisson_mesh, level_set);
        mesh_cutter.compute_cut_mesh_topology();
        mesh_cutter.compute_cut_mesh_geometry();
        mesh_cutter.mark_ghost_penalty_facets();

        File("interior_facet_marker.pvd") << *(mesh_cutter.interior_facet_marker()[0]);
        File("domain_marker.pvd") << *(mesh_cutter.domain_marker());
        File("interior_cut_cells.pvd") << *(mesh_cutter.cut_cells()[0]);
        File("interface_cells.pvd") << *(mesh_cutter.interface_cells());

        ////////////////////////////////////////////////////////////////////////////////////////////
        // Compute Normal                                                                         //
        ////////////////////////////////////////////////////////////////////////////////////////////
        Normal::CoefficientSpace_ls LS(poisson_mesh);
        Function level_set_inter(LS);
        level_set_inter.interpolate(*level_set);
        File("level_set_inter.pvd") << level_set_inter;

        Normal::FunctionSpace NV(poisson_mesh);
        Normal::BilinearForm na(NV,NV);
        Normal::LinearForm nL(NV);
        na.ls = level_set_inter;
        nL.ls = level_set_inter;

        Matrix An;
        Vector bn;

        assemble(An,na);
        assemble(bn,nL);

        Function normal(NV);

        LinearSolver solver("cg","hypre_amg");
        solver.parameters["monitor_convergence"] = true;
        solver.solve(An,*normal.vector(),bn);

        //Write out fido solution
        std::string n_fido_fname = "normal_fido.pvd";
        File n_fido_file(n_fido_fname.c_str());
        n_fido_file << normal;

        ////////////////////////////////////////////////////////////////////////////////////////////
        // Compute Poisson problem on FiDo
        ////////////////////////////////////////////////////////////////////////////////////////////

        //Source function
        Constant f(1.0);
        Constant g(0.0);
        Constant gamma(10.0);
        Constant beta(0.5);

        // Create standard forms
        Poisson::FunctionSpace V(poisson_mesh);
        Poisson::BilinearForm a0(V,V);
        a0.beta =  beta;
        a0.set_cell_domains(mesh_cutter.domain_marker());
        // Have to use facet marker 0 here.
        auto interior_facet_domains = mesh_cutter.interior_facet_marker()[0];
        a0.set_interior_facet_domains(interior_facet_domains);

        Poisson::LinearForm L0(V);
        L0.f = f;
        L0.set_cell_domains(mesh_cutter.domain_marker());

        // Create FidoForm starting with the standard form.
        FidoForm a_fd(a0);
        FidoForm L_fd(L0);

        // Create forms on cut cells
        PoissonCut::BilinearForm a_cut(V,V);
        PoissonCut::LinearForm L_cut(V);
        L_cut.f = f;

        // Get cut mesh
        std::size_t material_id = 0;
        auto cut_cells = mesh_cutter.cut_cells()[material_id];

        // Create quadrature rule
        std::size_t order = 1;
        Quadrature quadrature_a_cut(cut_cells->type().cell_type(),
                                    cut_cells->geometry().dim(),
                                    order);

        // Create cut cells forms
        a_fd.add_cut_form(boost::shared_ptr<CutForm>(new CutForm(a_cut,
                          *cut_cells,
                          quadrature_a_cut)));
        L_fd.add_cut_form(boost::shared_ptr<CutForm>(new CutForm(L_cut,
                          *cut_cells,
                          quadrature_a_cut)));

        // Create surface form
        PoissonSurface::BilinearForm a_s(V,V);
        PoissonSurface::LinearForm L_s(V);

        a_s.gamma = gamma;
        a_s.n = normal;

        L_s.g = g;
        L_s.gamma = gamma;
        L_s.n = normal;

        // Get cut mesh
        auto interface_cells = mesh_cutter.interface_cells();

        // Create quadrature rule
        Quadrature quadrature_a_s(interface_cells->type().cell_type(),
                                  interface_cells->geometry().dim(),
                                  order);

        // Create cut cells forms
        boost::shared_ptr<CutForm> cut_form_a_s(new CutForm(a_s,
                                                *interface_cells,
                                                quadrature_a_s));

        boost::shared_ptr<CutForm> cut_form_L_s(new CutForm(L_s,
                                                *interface_cells,
                                                quadrature_a_s));
        a_fd.add_cut_form(cut_form_a_s);
        L_fd.add_cut_form(cut_form_L_s);

        // Create constraints
        std::size_t domain_number = 2;

        a_fd.constrained_dofs() =
            // FIXME: Too many copies of constrained_dofs right now

        FidoTools::compute_constrained_dofs(V,
                                            mesh_cutter.domain_marker().get(),
                                            domain_number);


        //Assemble Nitsche
        Matrix A;
        Vector b;

        FidoAssembler fidoassembler;
        fidoassembler.assemble(b, L_fd);
        fidoassembler.assemble(A, a_fd);

        // Solve
        Function u(V);
        LinearSolver solver2("cg","hypre_amg");
        solver2.parameters["monitor_convergence"] = true;
        std::cout << "solving::" << std::endl;
        std::cout << A.str(true) << std::endl;
        std::cout << b.str(true) << std::endl;
        solver2.solve(A,*u.vector(),b);

        //Make plotting values of fido out fido solution
        std::vector<dolfin::la_index> dof_vertex_set = vertex_to_dof_map(*u.function_space());
        std::vector<double> values2(dof_vertex_set.size());
        u.vector()->get_local(values2.data(), dof_vertex_set.size(), dof_vertex_set.data());
        vert_values = values2;
        u_min = (u.vector()->min());
        u_max= (u.vector()->max());
    */
}
