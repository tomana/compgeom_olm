#pragma once

#include "ofMain.h"
#include "ofxGui.h"
#include "gl2ps.h"

class ofApp : public ofBaseApp{

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


		ofxPanel gui2;

        ofxSlider<float> tree_cwidth;
        ofxSlider<float> tree_cheight;
        ofxSlider<float> tree_swidth;
        ofxSlider<float> tree_sheight;
        ofxSlider<float> tree_center;
        ofxSlider<float> svg_n_slider;

        bool outputeps = false;
        bool outputtex = false;

};
