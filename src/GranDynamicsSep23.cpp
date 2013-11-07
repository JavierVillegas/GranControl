#include "testApp.h"

const int Nx = 352;
const int Ny = 288;




// A vector with all the frames.
vector <cv::Mat> TheFramesInput;

//vector<Grain> GrainList;
const int Nlast = 4;
ofVec3f BlockDims;

float inOp;//  overlaping porcentage
int Gs = 25; // for fixed grain size in all dimensions
cv::Mat CurrentOutput;
int fc=0; // frame counter
int ofc =0; // output frame counter

vector<ofVec2f > GrainsIN;
vector<ofVec4f > GrainsOut;
vector<ofVec2f > GrainsVelOut;
vector <float> GrainMeasurements;
vector <float> GrainRandAngs;

void UpdateMeasurements(int);
void GrainInit(void);

int OutFramesN;
bool UpdateBool;
float CalculateScalecorrect(void);
float ScaleCorrectFactor=1.0;
ofVec3f GrainCamera;
ofVec3f GrainUpVector;
int outYsize = Ny;
int outXsize = Nx;

int viewerDir =0;
int circuIndex = 0;




float dt = 0.9;
float Kback = 0.1097;
float Kdam =0.0567;
float Krand = 0.042;//
float MTh = 1.3; // motion threshold

// TODO map as globals to find the sweet operational point

//--------------------------------------------------------------
void testApp::setup(){
   
    // First creating an empty array
   
    BlockDims.x = Nx;
    BlockDims.y = Ny;
    BlockDims.z = Nlast;
    cv::Mat AuxMat(Ny,Nx,CV_8UC3,cv::Scalar(0,0,0));
    for (int k = 0; k< Nlast; k++){
          TheFramesInput.push_back(AuxMat.clone());
    }
    

     
     //initial value
    OutFramesN = BlockDims.z;
    
 //   vidGrabber.setVerbose(true);
 //   vidGrabber.initGrabber(Nx,Ny);
    
    ofSetFrameRate(15);
    inOp =50;
  // the input grain list 
    GrainInit();

     ScaleCorrectFactor =CalculateScalecorrect();
    UpdateBool = true;
    vidGrabber.initGrabber(Nx,Ny);
    colorImg.allocate(Nx,Ny);
	grayImage.allocate(Nx,Ny);
}

//--------------------------------------------------------------
void testApp::update(){
    
   
    bool bNewFrame = false;
    
    vidGrabber.update();
    bNewFrame = vidGrabber.isFrameNew();
    //
    if (bNewFrame){
        //
        //
        colorImg.setFromPixels(vidGrabber.getPixels(), Nx,Ny);
        cv::Mat AuxMat;
        AuxMat = colorImg.getCvImage();
        AuxMat.copyTo(TheFramesInput[circuIndex]);
        circuIndex++;
        if (circuIndex==Nlast){
            circuIndex = 0;
        }
    
    }
    
    CurrentOutput =cv::Mat::zeros(outYsize,outXsize, CV_8UC3);

    unsigned char *output = (unsigned char*)(CurrentOutput.data);
    
    

    
    
    
    // input/output jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
     
//   if (((fc%jp) == jp-1 )&&(fc < (Nlast-Gs))) {
    if ((fc < (Nlast-1.0))) {

        UpdateMeasurements(fc);
        
    }
    
    
    
    for (int g =0; g < GrainsIN.size(); g++) {
//        if (((fc%jp) == jp-1 )&&(fc < (Nlast-Gs))) {
 //      if (((fc%60) == 0 )&&(fc < (Nlast-Gs))) {
        
   

//            GrainsVelOut[g].x=0;
//            GrainsVelOut[g].y = 0;
 //      }

            
             unsigned char *input = (unsigned char*)(TheFramesInput[fc].data);
        
        
        // for the 2 overlaping windows in time 
        for (int m = 1; m<2; m++) {
        
            for (int x =0; x<Gs; x++) {
                for (int y=0; y<Gs; y++) {
                    
                    float rIn,gIn,bIn;
                    float rOut,gOut,bOut;
                    float Scale;
                    Scale = ((0.5*0.5*0.5)*(1.0 -cosf(2*PI*x/(float)(Gs-1)))*
                             (1.0 -cosf(2*PI*y/(float)(Gs-1)))  );//*
                             //(1.0 -cosf(2*PI*(m*(Gs-1)/2.0 +(fc%jp))/(float)(Gs-1))));
                    
                    Scale/=(1.1*ScaleCorrectFactor*ScaleCorrectFactor);
                    // pixels at the input
                    int yindIn = y + GrainsIN[g].y;
                    int xindIn = x + GrainsIN[g].x;
                    
                  
                  
                    bIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn)) ] ;
                    gIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 1)];
                    rIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 2)];
                    //pixels at the output
                    // output indexes
                    int xindOut,yindOut;
                    if (m==0){
                        xindOut = x + GrainsOut[g].x;
                        yindOut = y + GrainsOut[g].y;
                    }
                    else{
                        xindOut = x + GrainsOut[g].z;
                        yindOut = y + GrainsOut[g].w;
                    }
                    
                    xindOut = (xindOut<0)?0:xindOut;
                    xindOut = (xindOut>outXsize-1)?outXsize-1:xindOut;
                    
                    yindOut = (yindOut<0)?0:yindOut;
                    yindOut = (yindOut>outYsize-1)?outYsize-1:yindOut;
                    
                    bOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] ;
                    gOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1];
                    rOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2];
                    

                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] = (uchar)(bOut + Scale*bIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1]=(uchar)(gOut +Scale*gIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2]=(uchar)(rOut + Scale*rIn);
//                    float sca2 = GrainMeasurements[g]/100.0;
//                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] = (uchar)(bOut + Scale*bIn*sca2);
//                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1]=(uchar)(gOut +Scale*gIn*sca2);
//                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2]=(uchar)(rOut + Scale*rIn*sca2);
                    
                    
                    
                
                
                }// for y
            } // for x
        } // for of the 2 halves
        
        // Updating the halves
//        if ((fc%jp) == jp-1 ) {
//            GrainsOut[g].x = GrainsOut[g].z;
//            GrainsOut[g].y = GrainsOut[g].w;
//        }
        
            // sum of forces, acceleration
            
      

        
        
        
        ofVec2f NormaDir(GrainsIN[g].x - GrainsOut[g].z,GrainsIN[g].y - GrainsOut[g].w);
        NormaDir.normalize();
        
//            float accX = Kback*(GrainsIN[g].x - GrainsOut[g].z) + Krand*(cos(GrainRandAngs[g]))
//            - Kdam*GrainsVelOut[g].x;
//            float accY = Kback*(GrainsIN[g].y - GrainsOut[g].w) + Krand*(sin(GrainRandAngs[g]))
//            - Kdam*GrainsVelOut[g].y;
        
        
        float disturvance;
        if(GrainMeasurements[g]>MTh){
        
            disturvance = GrainMeasurements[g];
            
        }
        else{
            disturvance =0;

        }
        
        float accX = Kback*(NormaDir.x) + Krand*disturvance*(cos(GrainRandAngs[g]))
        - Kdam*GrainsVelOut[g].x;
        float accY = Kback*(NormaDir.y) + Krand*disturvance*(sin(GrainRandAngs[g]))
        - Kdam*GrainsVelOut[g].y;
            // velocity;
        

            GrainsVelOut[g].x += accX*dt;
            GrainsVelOut[g].y += accY*dt;
            
            // Position
            
            GrainsOut[g].z  += GrainsVelOut[g].x*dt;
            GrainsOut[g].w  += GrainsVelOut[g].y*dt;

        
        // if is close enought fixed it!
        
        float distagrain = (GrainsOut[g].z-GrainsIN[g].x)*(GrainsOut[g].z-GrainsIN[g].x) +
        (GrainsOut[g].w-GrainsIN[g].y)*(GrainsOut[g].w-GrainsIN[g].y);
        if((distagrain < 1.0)&&(disturvance==0)){
            GrainsOut[g].z = GrainsIN[g].x;
            GrainsOut[g].w = GrainsIN[g].y;
        }
        
//            GrainsOut[g].z   += (int)(cos(GrainRandAngs[g])*sca2);
//            GrainsOut[g].w  += (int)(sin(GrainRandAngs[g])*sca2);
       
        
    }// for g
    

    
}




//--------------------------------------------------------------
void testApp::draw(){
	ofSetHexColor(0xffffff);
    ofxCvColorImage AuxDrawImage;

    AuxDrawImage.allocate(TheFramesInput[fc].cols, TheFramesInput[fc].rows);
    AuxDrawImage = TheFramesInput[fc].data;
    AuxDrawImage.draw(0, 0,640,480);
    
    ofxCvColorImage AuxDrawImage2;
    
    AuxDrawImage2.allocate(CurrentOutput.cols, CurrentOutput.rows);
    AuxDrawImage2 = CurrentOutput.data;
    AuxDrawImage2.draw(0, 480,640,480);
   

    
    
    ofc++;
    
    if(ofc>OutFramesN){ofc = 0;}
    fc++;
    if(fc>BlockDims.z-1){fc=0;}
    
}




void GrainInit(){
// Initialices the grain scheduler for the input, and the output positions and velocities
    GrainsIN.clear();
    GrainsOut.clear();
    GrainsVelOut.clear();
    GrainMeasurements.clear();
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
       
    
    for (int y = 0; y < (BlockDims.y - Gs); y+=jp) {
        for (int x=0; x< (BlockDims.x - Gs); x+=jp) {
                ofVec2f tempStorage;
                tempStorage.x =x;
                tempStorage.y = y;
                GrainsIN.push_back(tempStorage);
                ofVec4f temvecout;
                temvecout.x = x;
                temvecout.y = y;
                temvecout.z = x;
                temvecout.w = y;
                GrainsOut.push_back(temvecout);
                GrainsVelOut.push_back(ofVec2f(0.0,0.0));
            GrainMeasurements.push_back(0.0);
            
            float randa = ofRandom(2*PI);
            GrainRandAngs.push_back(randa);
            
   
            }
        }
   

}



void UpdateMeasurements(int initframe){

    
    float rIn,gIn,bIn;
    
    float Ns =2.0; //= ((float)(Gs));
    for (int g =0; g < GrainsIN.size(); g++) {
    

            float StdAvg =0;
            for (int x =0; x<Gs; x++) {
                for (int y=0; y<Gs; y++) {
                
                    int yindIn = y + GrainsIN[g].y;
                    int xindIn = x + GrainsIN[g].x;
                    float Theavg=0;
                    float TheMSE =0;
                    for (int t = 0; t < Ns; t++) {
                         unsigned char *input = (unsigned char*)(TheFramesInput[initframe + t].data);
                         bIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn)) ] ;
                         gIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 1)];
                         rIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 2)];
                        // gray value pixel
                        float val = 0.229*rIn + 0.587*gIn + 0.114*bIn;
                        Theavg+= val;
                        TheMSE += val*val;
                    } // for t
                   // Variance for this x,y

                    StdAvg+= sqrtf(abs(TheMSE/(float)Ns -(Theavg/(float)Ns)*(Theavg/(float)Ns)));
                    
             }// for y
        } // for x
  
        
        GrainMeasurements[g] = StdAvg/((float)Gs*(float)Gs);
      }
}




float CalculateScalecorrect(void){

    vector<float> Tempovec(10*Gs);
    // initializing
    for (int k=0; k<(10*Gs); k++) {
        Tempovec[k]=0.0;
    }
    
    
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
    
    // overlaping
    for (int k=0; k<(9*Gs); k+=jp) {
        for (int n=0; n<Gs; n++) {
            Tempovec[k+n]+= 0.5*(1.0 -cosf(2*PI*n/(float)(Gs-1)));
        }
    }

    // finding the max;
    float TheMax =0;
    for (int k=0; k<(10*Gs); k++) {
        if(Tempovec[k]>TheMax){
            TheMax = Tempovec[k];
        }
    }

    return TheMax;

}













//--------------------------------------------------------------
void testApp::keyPressed(int key){
    
    switch (key) {
        case 't':
            MTh +=0.1;
            cout<<"MTH:"<<MTh<<endl;
            break;
        case 'g':
            MTh -=0.1;
            cout<<"MTH:"<<MTh<<endl;
            break;
        case 'r':
            GrainInit();
            break;
        case OF_KEY_RIGHT:
            Kback+=.0001;
            cout<<"Kback:"<<Kback<<endl;
            break;
        case OF_KEY_LEFT:
            Kback-=.0001;
             cout<<"Kback:"<<Kback<<endl;
            break;
        case OF_KEY_UP:
            Krand+=.001;
             cout<<"Krand:"<<Krand<<endl;
            break;
        case OF_KEY_DOWN:
            Krand-=.001;
              cout<<"Krand:"<<Krand<<endl;
            break;
        
        
        case 'e':
            break;
            
        case 'a':
            Kdam +=0.0001;
            cout<<"Kdam:"<<Kdam<<endl;
            break;
        case 'z':
            Kdam -=0.0001;
            cout<<"Kdam:"<<Kdam<<endl;
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