#include "ofApp.h"
#define STRINGIFY(X) #X

//--------------------------------------------------------------
void ofApp::setup(){
    gui2.setup("Draw u");
    gui2.add(tree_cwidth.setup("tree_cwidth", 0.01, 0, 10000 ));
    gui2.add(tree_cheight.setup("tree_cheight", 0.01, 0, 1000 ));
    gui2.add(tree_swidth.setup("tree_swidth", 0.01, 0, 1000 ));
    gui2.add(tree_sheight.setup("tree_sheight", 0.01, 0, 1000 ));
    gui2.add(tree_center.setup("tree_center", 0.01, 0, 1000 ));
    gui2.add(svg_n_slider.setup( "svg_n_slider", 0.5, 0, 10 ));
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
ofBackground(255);

 gl2psLineWidth(0.5);
gl2psEnable(GL2PS_BLEND);
    const std::string treestring =
        STRINGIFY(
2llllllrl1
2llllllrr1
2lllllll1
2lllllrll1
2lllllrlr1
2lllllrrl1
2lllllrrr1
2llllrlll1
2llllrllr1
2llllrlrl1
2llllrlrr1
2llllrrll1
2llllrrlr1
2llllrrrl1
2llllrrrr1
2lllrllrl1
2lllrllrr1
2lllrlll1
2lllrlrll1
2lllrlrlr1
2lllrlrrl1
2lllrlrrr1
2lllrrlll1
2lllrrllr1
2lllrrlrl1
2lllrrlrr1
2lllrrrll1
2lllrrrlr1
2lllrrrrl1
2lllrrrrr1
2llrlllrl1
2llrlllrr1
2llrllll1
2llrllrll1
2llrllrlr1
2llrllrrl1
2llrllrrr1
2llrlrlll1
2llrlrllr1
2llrlrlrl1
2llrlrlrr1
2llrlrrll1
2llrlrrlr1
2llrlrrrl1
2llrlrrrr1
2llrrllrl1
2llrrllrr1
2llrrlll1
2llrrlrll1
2llrrlrlr1
2llrrlrrl1
2llrrlrrr1
2llrrrlll1
2llrrrllr1
2llrrrlrl1
2llrrrlrr1
2llrrrrll1
2llrrrrlr1
2llrrrrrl1
2llrrrrrr1
2lrllllrl1
2lrllllrr1
2lrlllll1
2lrlllrll1
2lrlllrlr1
2lrlllrrl1
2lrlllrrr1
2lrllrlll1
2lrllrllr1
2lrllrlrl1
2lrllrlrr1
2lrllrrll1
2lrllrrlr1
2lrllrrrl1
2lrllrrrr1
2lrlrllrl1
2lrlrllrr1
2lrlrlll1
2lrlrlrll1
2lrlrlrlr1
2lrlrlrrl1
2lrlrlrrr1
2lrlrrlll1
2lrlrrllr1
2lrlrrlrl1
2lrlrrlrr1
2lrlrrrll1
2lrlrrrlr1
2lrlrrrrl1
2lrlrrrrr1
2lrrlllrl1
2lrrlllrr1
2lrrllll1
2lrrllrll1
2lrrllrlr1
2lrrllrrl1
2lrrllrrr1
2lrrlrlll1
2lrrlrllr1
2lrrlrlrl1
2lrrlrlrr1
2lrrlrrll1
2lrrlrrlr1
2lrrlrrrl1
2lrrlrrrr1
2lrrrllll1
2lrrrlllr1
2lrrrllrl1
2lrrrllrr1
2lrrrlrll1
2lrrrlrlr1
2lrrrlrrl1
2lrrrlrrr1
2lrrrrlll1
2lrrrrllr1
2lrrrrlrl1
2lrrrrlrr1
2lrrrrrll1
2lrrrrrlr1
2lrrrrrrl1
2lrrrrrrr1
2rlllllrl1
2rlllllrr1
2rllllll1
2rllllrll1
2rllllrlr1
2rllllrrl1
2rllllrrr1
2rlllrlll1
2rlllrllr1
2rlllrlrl1
2rlllrlrr1
2rlllrrll1
2rlllrrlr1
2rlllrrrl1
2rlllrrrr1
2rllrllrl1
2rllrllrr1
2rllrlll1
2rllrlrll1
2rllrlrlr1
2rllrlrrl1
2rllrlrrr1
2rllrrlll1
2rllrrllr1
2rllrrlrl1
2rllrrlrr1
2rllrrrll1
2rllrrrlr1
2rllrrrrl1
2rllrrrrr1
2rlrlllrl1
2rlrlllrr1
2rlrllll1
2rlrllrll1
2rlrllrlr1
2rlrllrrl1
2rlrllrrr1
2rlrlrlll1
2rlrlrllr1
2rlrlrlrl1
2rlrlrlrr1
2rlrlrrll1
2rlrlrrlr1
2rlrlrrrl1
2rlrlrrrr1
2rlrrllll1
2rlrrlllr1
2rlrrllrl1
2rlrrllrr1
2rlrrlrll1
2rlrrlrlr1
2rlrrlrrl1
2rlrrlrrr1
2rlrrrlll1
2rlrrrllr1
2rlrrrlrl1
2rlrrrlrr1
2rlrrrrll1
2rlrrrrlr1
2rlrrrrrl1
2rlrrrrrr1
2rrllllrl1
2rrllllrr1
2rrlllll1
2rrlllrll1
2rrlllrlr1
2rrlllrrl1
2rrlllrrr1
2rrllrlll1
2rrllrllr1
2rrllrlrl1
2rrllrlrr1
2rrllrrll1
2rrllrrlr1
2rrllrrrl1
2rrllrrrr1
2rrlrllrl1
2rrlrllrr1
2rrlrlll1
2rrlrlrll1
2rrlrlrlr1
2rrlrlrrl1
2rrlrlrrr1
2rrlrrlll1
2rrlrrllr1
2rrlrrlrl1
2rrlrrlrr1
2rrlrrrll1
2rrlrrrlr1
2rrlrrrrl1
2rrlrrrrr1
2rrrlllrl1
2rrrlllrr1
2rrrllll1
2rrrllrll1
2rrrllrlr1
2rrrllrrl1
2rrrllrrr1
2rrrlrlll1
2rrrlrllr1
2rrrlrlrl1
2rrrlrlrr1
2rrrlrrll1
2rrrlrrlr1
2rrrlrrrl1
2rrrlrrrr1
2rrrrllll1
2rrrrlllr1
2rrrrllrl1
2rrrrllrr1
2rrrrlrll1
2rrrrlrlr1
2rrrrlrrl1
2rrrrlrrr1
2rrrrrlll1
2rrrrrllr1
2rrrrrlrl1
2rrrrrlrr1
2rrrrrrll1
2rrrrrrlr1
2rrrrrrrl1
2rrrrrrrr1
                 );

    float swidth = tree_swidth;
    float sheight = tree_sheight;
    float center = tree_center;
    float cheight = tree_cheight;
    float cwidth = tree_cwidth;
    float depth = 0;

    ofSetColor(0);
    int counter_total = 0;
    for(int i = 0; i < treestring.size(); i++)
    {
        if(treestring[i] == '2') {
         cheight = tree_cheight;
         depth = 1;
         }

        if(treestring[i] == 'r') {
        ofLine(cwidth, cheight, cwidth + (1.0/(pow(2,depth)))*swidth, cheight + sheight);
        cwidth += (1.0/(pow(2,depth)))*swidth;
        cheight += sheight;
        depth +=1;
        }
        if(treestring[i] == 'l') {
        ofLine(cwidth, cheight, cwidth -  (1.0/(pow(2,depth)))*swidth, cheight +  sheight);
        cwidth -= (1.0/(pow(2,depth)))*swidth;
        cheight += sheight;
        depth +=1;
        }

        if(treestring[i] == '1')
        {
            counter_total++;
            ofRect(cwidth - 1, cheight -1, 2, 2);
            cheight = tree_cheight;
            cwidth = tree_cwidth;
            depth = 0;
        }

    }

    if(!outputeps) gui2.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
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
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){

}
