#include "ofMain.h"
#include "ofxGrabCam.h"
#include "ofxGui.h"
#include "gl2ps.h"

// Various undefs and defs
#undef restrict
#undef Modifiers
#undef SHIFT
#undef Status
#undef Success

#define HAS_CGAL
#define HAS_MPI

// Dolfin stuff
#include <dolfin.h>
#include "/home/tom/Desktop/FEniCS-cutfem/cutfem/local.naerland.overlappingmeshes/include/cutfem.h"

#include "poisson/PoissonOLM0.h"
#include "poisson/PoissonOLM1.h"
#include "poisson/PoissonOLMInterface.h"
#include "poisson/PoissonOLM0_p2.h"
#include "poisson/PoissonOLM1_p2.h"
#include "poisson/PoissonOLMInterface_p2.h"
#include "poisson/Poisson.h"
#include "poisson/PoissonNitscheBC.h"
#include "poisson/H1_Norm_u_0.h"
#include "poisson/H1_Norm_u_0_only.h"
#include "poisson/H1_Norm_u_1.h"
#include "poisson/P3.h"


class ofApp : public ofBaseApp
{

public:
    void setup();
    void update();
    void draw();

    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);



    dolfin::Mesh * mesh;
    dolfin::Mesh * collidingMesh;

    void testVolume();

    void draw_or_save_string(float x, float y, string indexstring);
    std::vector<size_t> draw_indexes_mesh_vec;
    std::vector<size_t> draw_indexes_collidingMesh_vec;
    std::vector<size_t> draw_indexes_subtriang_mesh_vec;
    std::vector<size_t> draw_indexes_facettriang_collidingMesh_vec;

    void clearsets();

    void drawTet(dolfin::Point v0,dolfin::Point v1,dolfin::Point v2, dolfin::Point v3);
    void drawTetLines3D(dolfin::Point v0,dolfin::Point v1,dolfin::Point v2, dolfin::Point v3);
    void drawTetLines2D(dolfin::Point v0,dolfin::Point v1,dolfin::Point v2);
    void drawLines2D(dolfin::Point v0,dolfin::Point v1);
    void drawBoxPoint2D(dolfin::Point v0);
    void drawTriangleColor2D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2, ofColor c0, ofColor c1, ofColor c2);
    void drawTriangleLinesColor2D(dolfin::Point v0,dolfin::Point v1, dolfin::Point v2, ofColor c0, ofColor c1, ofColor c2);

    // DRAW HELPERS
    ofxGrabCam cam;
    bool drawsomething;
    bool outputeps;
    bool outputtex;
    ofMatrix4x4 savedPose;

    // GUI
    ofxPanel gui;
    ofxPanel gui2;
    ofxToggle draw_indexes;
    ofxToggle draw_first_mesh;
    ofxToggle draw_second_mesh;
    ofxToggle draw_t_0_1;
    ofxToggle draw_t_0_gamma;
    ofxToggle draw_t_0_2;
    ofxToggle draw_border_first;
    ofxToggle draw_triangulation_mesh;
    ofxToggle draw_index_triangulation;
    ofxToggle draw_first_mesh_u;
    ofxToggle draw_first_mesh_t_0_u;
    ofxToggle draw_first_mesh_t_0_t_1_u;
    ofxToggle draw_second_mesh_u;
    ofxToggle draw_triangulation_mesh_u;

    ofxButton mesh2dbutton;
    ofxButton mesh2dgenbutton;

    int lastmeshu_r;
    ofxSlider<int> meshu_r;
    ofxSlider<int> meshx_r;
    ofxSlider<int> meshy_r;
    ofxSlider<int> meshz_r;

    int lastcmeshu_r;
    ofxSlider<int> cmeshu_r;
    ofxSlider<int> cmeshx_r;
    ofxSlider<int> cmeshy_r;
    ofxSlider<int> cmeshz_r;


    ofxSlider<float> zaxisscaler;
    ofxSlider<int> transparency;
    ofxSlider<float> svg_n_slider;

    ofxSlider<float> allscaler;
    ofxSlider<float> rotatez;
    ofxSlider<float> trans_u_z;
    ofxSlider<float> trans_first_mesh;
    ofxSlider<float> trans_second_mesh;
    ofxSlider<float> trans_t_0_1;
    ofxSlider<float> trans_t_0_gamma;
    ofxSlider<float> trans_t_0_2;
    ofxSlider<float> trans_triangulation_mesh;
    ofxSlider<float> trans_border_first;

    ofxSlider<float> gamma_u;

    ofParameter<ofColor> firstcolor;
    ofParameter<ofColor> secondcolor;

    void meshuniform();
    void cmeshuniform();

    void mesh2dbuttonPressed();
    void mesh2dgenbuttonPressed();

    bool checksets;
    float magicnumber;

    void refineMeshes();
    void newStructuredMeshes(size_t n);

    // POISSON
    void runPoisson();
    void runPoissonNitscheBC();

    void runPoissonFiDo();
    void runPoissonOLM();
    void runPoissonOLM_convergence_p1();
    void runPoissonOLM_convergence_p2();
    void runPoissonOLM_patch_test_p1();
    void runPoissonOLM_patch_test_p2();

    std::shared_ptr<cutfem::CompositeMesh> overlapping_meshes;
    void plotOnMesh(dolfin::Function u);

    std::vector<std::vector<double> > cell_values;
    std::vector<double> vert_values_1;
    std::vector<double> vert_values_2;
    std::vector<pair<double, double> > t_0_strings;
    std::vector<pair<double, double> > t_2_strings;
    std::vector<pair<double, double> > combined_strings;

    std::vector<pair<double, double> > t_0_strings_p2;
    std::vector<pair<double, double> > t_2_strings_p2;
    std::vector<pair<double, double> > combined_strings_p2;

    std::vector<pair<double, double> > timing_bvh;
    std::vector<pair<double, double> > timing_triang;
    std::vector<pair<double, double> > timing_wo_cg;
    std::vector<pair<double, double> > timing_total;
    std::vector<pair<double, double> > timing_standard_poisson;
    std::vector<pair<double, double> > size_mesh1;
    std::vector<pair<double, double> > size_mesh2;
    std::vector<pair<double, double> > size_triangulation;
    std::vector<pair<double, double> > size_interface;

    int n_structured = 0;
    double u_min_1;
    double u_max_1;

    double u_min_2;
    double u_max_2;

    // Interpolation for triangle
double interpolate_triangle(dolfin::Point p0,
                                   double w1,
                                   dolfin::Point p1,
                                   double w2,
                                   dolfin::Point p2,
                                   double w3,
                                   dolfin::Point pi);

bool do_intersect_two_lines(dolfin::Point p1, dolfin::Point p2, dolfin::Point p3, dolfin::Point p4);

bool point_on_line(dolfin::Point linePointA, dolfin::Point linePointB, dolfin::Point point);

};
