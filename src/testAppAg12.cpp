#include "testApp.h"



// A vector with all the frames.
vector <cv::Mat> TheFramesInput;
vector<Grain> GrainList; 
const int Nlast = 1065;
ofVec3f BlockDims;
float Op;//  overlaping porcentage
cv::Mat CurrentOutput;
int fc=0; // frame counter
//--------------------------------------------------------------
void testApp::setup(){
   
    // loading all the Paris data
    // To create the input video block
    for (int k = 1; k<= Nlast; k++){
    
        ofImage auxImage;
        auxImage.loadImage("ParisFrames/paris"+ofToString(k)+".jpg");
        BlockDims.x = auxImage.width;
        BlockDims.y = auxImage.height;
        BlockDims.z = Nlast;
        
        ofxCvColorImage AuxCVImage;
        AuxCVImage.allocate(auxImage.width,auxImage.height);
        AuxCVImage.setFromPixels(auxImage.getPixels(),auxImage.width,auxImage.height);
        cv::Mat AuxMat(AuxCVImage.getCvImage());
        
        TheFramesInput.push_back(AuxMat.clone());
    
    }
    
 //   vidGrabber.setVerbose(true);
 //   vidGrabber.initGrabber(Nx,Ny);
    
    ofSetFrameRate(15);
   // FillGrainList(15);
    Op =50;
}

//--------------------------------------------------------------
void testApp::update(){
    
    CurrentOutput.create(BlockDims.y, BlockDims.x, CV_8UC3);
    
    // calculating the current Output
    unsigned char *input = (unsigned char*)(TheFramesInput[fc].data);
    unsigned char *output = (unsigned char*)(CurrentOutput.data);
    for (int x = 0; x < BlockDims.x; x++) {
        for (int y = 0; y < BlockDims.y; y++) {
//    for (int x = 0; x < 110; x++) {
//        for (int y = 0; y < 110; y++) {
            float Scale;
            vector<Grain> GC = ListGrainsThatContribute(x,y,fc,21);
            
//                for (int s = 0; s < GrainList.size(); s++) {
//                    Scale+=GrainList[s].CalculateContrib(x, y, fc);
//                 }
            Scale = 0.0;
            for (int s= 0; s< GC.size(); s++) {
                Scale += GC[s].CalculateContrib(x, y, fc);
//                cout<<"x =" <<x<<"y = "<<y<<" || IndX = "<<GC[s].GrInd.x<<", IndY ="<<GC[s].GrInd.y<<" IndZ ="<<GC[s].GrInd.z<<endl;
            }
           
                 float r,g,b;
//                    b = input[(int)BlockDims.x * y + x ] ;
//                    g = input[(int)BlockDims.x * y + x + 1];
//                    r = input[(int)BlockDims.x * y + x + 2];
//            
//                    output[(int)BlockDims.x * y + x ] = (uchar)(Scale*b);
//                    output[(int)BlockDims.x * y + x +1] = (uchar)(Scale*g);
//                    output[(int)BlockDims.x * y + x +2] = (uchar)(Scale*r);
            b = input[3*TheFramesInput[fc].cols * y + 3*x ] ;
            g = input[3*TheFramesInput[fc].cols * y + 3*x + 1];
            r = input[3*TheFramesInput[fc].cols * y + 3*x + 2];
            
            output[3*CurrentOutput.cols * y + 3*x ] = (uchar)(Scale*b);
            output[3*CurrentOutput.cols * y + 3*x +1] = (uchar)(Scale*g);
            output[3*CurrentOutput.cols * y + 3*x +2] = (uchar)(Scale*r);
           
            GC.clear();
        }
      }
    
    // bool bNewFrame = false;
    
    
    //vidGrabber.update();
    //bNewFrame = vidGrabber.isFrameNew();
    
	//if (bNewFrame){
    // }
    
}




//--------------------------------------------------------------
void testApp::draw(){
	ofSetHexColor(0xffffff);
    ofxCvColorImage AuxDrawImage;

    AuxDrawImage.allocate(TheFramesInput[fc].cols, TheFramesInput[fc].rows);
    AuxDrawImage = TheFramesInput[fc].data;
    AuxDrawImage.draw(0, 0);
    
    ofxCvColorImage AuxDrawImage2;
    
    AuxDrawImage2.allocate(CurrentOutput.cols, CurrentOutput.rows);
    AuxDrawImage2 = CurrentOutput.data;
    AuxDrawImage2.draw(AuxDrawImage.width, 0);
   
    
    fc++;
    if(fc>Nlast-1){fc=0;}
    
}


void testApp::FillGrainList(float GrainSize){

    int Osamp = (int)(ceil(GrainSize*Op/100.0));
    int jp = GrainSize - Osamp;
    
    for (int x=0; x< BlockDims.x; x+=jp) {
        for (int y = 0; y < BlockDims.y; y+=jp) {
            for (int t=0; t < BlockDims.z; t+=jp) {
              GrainList.push_back(Grain(x,y,t,GrainSize));
            }
        }
    }
}


vector<Grain> testApp::ListGrainsThatContribute(float x, float y, float t, float GrainSize){

    vector<Grain> TobeReturned;
    TobeReturned.clear();
    int Osamp = (int)(ceil(GrainSize*Op/100.0));
    int jp = GrainSize - Osamp;
    vector<int> Indx;
    int TemInd = floor(x/(float)jp);
    Indx.push_back(jp*TemInd);
    if (((x-jp*TemInd) < Osamp)&&(TemInd!=0)) {
        Indx.push_back(jp*(TemInd-1));
    }
    
    vector<int> Indy;
    TemInd = floor(y/(float)jp);
    Indy.push_back(jp*TemInd);
    if (((y-jp*TemInd) < Osamp)&&(TemInd!=0)) {
        Indy.push_back(jp*(TemInd-1));
    }
    
    vector<int> Indt;
    TemInd = floor(t/(float)jp);
    Indt.push_back(jp*TemInd);
    if (((t-jp*TemInd) < Osamp)&&(TemInd!=0)) {
        Indt.push_back(jp*(TemInd-1));
    }

    for (int k=0; k<Indx.size(); k++) {
        for (int q=0; q<Indy.size(); q++) {
            for (int r=0; r < Indt.size(); r++) {
                TobeReturned.push_back(Grain(Indx[k],Indy[q],Indt[r],GrainSize));
            }
        }
    }
    
    return TobeReturned;
    
    
    
}









//--------------------------------------------------------------
void testApp::keyPressed(int key){
    
    switch (key) {
        case 'm':
            

            break;
            
        case 's':
      
            break;
        case OF_KEY_RIGHT:
            break;
        case OF_KEY_LEFT:
            
            break;
        case 'e':
            break;
            
        case 'g':
            break;
        case 'b':
            break;
        default:
            break;
    }
    
    
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y){
    
}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
    
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
    
}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 
    
}