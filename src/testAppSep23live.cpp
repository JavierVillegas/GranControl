#include "testApp.h"


const int Nx = 352;
const int Ny = 288;
// A vector with all the frames.
vector <cv::Mat> TheFramesInput;

//vector<Grain> GrainList;
const int Nlast = 700;
ofVec3f BlockDims;
float inOp;//  overlaping porcentage
float outOp;//  overlaping porcentage
float NewoutOp;
int Gs = 51; // for fixed grain size in all dimensions
cv::Mat CurrentOutput;
int fc=0; // frame counter
int ofc =0; // output frame counter

vector<vector<int> > GrainScheduler;
void UpdateScheduler(void);
void UpdateSchedulerUpDown(int Upf, int Downf);
int OutFramesN;
bool UpdateBool;
float CalculateScalecorrect(void);
float ScaleCorrectFactor=1.0;
int TheUf =1;

bool doRotation = false;
void GrainRandomizer();
bool boolshufflegrains = false;
bool boolfrozesomegrains = false;
int circuIndex = Nlast-1;
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
    outOp =50;
    NewoutOp = outOp;

   // UpdateScheduler();
    UpdateSchedulerUpDown(TheUf,1);
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
        fc = circuIndex;
        circuIndex--;
        if (circuIndex<0){
            circuIndex = Nlast-1;
        }
        
    

    

    if (UpdateBool){
        outOp = NewoutOp;
      //UpdateScheduler();
        UpdateSchedulerUpDown(TheUf,1);
        if (boolshufflegrains==true){
            GrainRandomizer();}
        ScaleCorrectFactor =CalculateScalecorrect();
        UpdateBool = false;
    }
    // calculate new output size:
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
 
    // input jump size
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;

    
    int outYsize = Gs + TheUf*ojp*((BlockDims.y-Gs)/jp);
    int outXsize = Gs + TheUf*ojp*((BlockDims.x-Gs)/jp);
    // calculate the number of output frames.
    OutFramesN = Gs + TheUf*ojp*((BlockDims.z-Gs)/jp);
    
   // creating an empty frame for the output

    CurrentOutput =cv::Mat::zeros(outYsize,outXsize, CV_8UC3);
 
    // pointer to output data:
    
    unsigned char *output = (unsigned char*)(CurrentOutput.data);
    
    // running throught the list of grains
    
    for (int g =0; g < GrainScheduler.size(); g++) {
        
        // check if in the current frame the current grain needs to be ploted
        
        int indT = ofc - GrainScheduler[g][5];
        if ((indT >=0)&&(indT<Gs)) {
          // copy the grain from the input
          // getting a pointer to the input frame
//            cout<<fc<<endl;
//            cout<<GrainScheduler[g][2]<<endl;
//            cout<<indT<<endl;
//            cout<<"mod:"<<((fc + GrainScheduler[g][2])%Nlast)<<endl;
            //((fc + GrainScheduler[g][2]+indT)%Nlast)
             unsigned char *input = (unsigned char*)(TheFramesInput[GrainScheduler[g][2] + indT].data);
            for (int x =0; x<Gs; x++) {
                for (int y=0; y<Gs; y++) {
                    
                    float rIn,gIn,bIn;
                    float rOut,gOut,bOut;
                    float Scale;
                    Scale = ((0.5*0.5*0.5)*(1.0 -cosf(2*PI*x/(float)(Gs-1)))*
                             (1.0 -cosf(2*PI*y/(float)(Gs-1)))*
                             (1.0 -cosf(2*PI*indT/(float)(Gs-1))));
                    
                    Scale/=(ScaleCorrectFactor*ScaleCorrectFactor*ScaleCorrectFactor);
                    // pixels at the input
                    int yindIn = y + GrainScheduler[g][1];
                    int xindIn = x + GrainScheduler[g][0];
                    bIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn)) ] ;
                    gIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 1)];
                    rIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 2)];
                    //pixels at the output
                    // output indexes
                    int xindOut = x + GrainScheduler[g][3];
                    int yindOut = y + GrainScheduler[g][4];
                    
                    
                    bOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] ;
                    gOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1];
                    rOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2];
                    

                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] = (uchar)(bOut + Scale*bIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1]=(uchar)(gOut +Scale*gIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2]=(uchar)(rOut + Scale*rIn);
                    
                    
                
                
                }
            }
            
            
            
            
            
        }
    }
    
    }
    
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
    AuxDrawImage2.draw(0, AuxDrawImage.height);
   

    
    
    ofc++;
    if(ofc>OutFramesN){ofc=0;}
    fc++;
    if(fc>BlockDims.z-1){fc=0;}
    
}


void UpdateScheduler(){
// calculates the grains output position, creates a table of grain input output pairs.
    GrainScheduler.clear();
    
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
    // output jump
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
    
    
    for (int y = 0; y < (BlockDims.y - Gs); y+=jp) {
        for (int x=0; x< (BlockDims.x - Gs); x+=jp) {
            for (int t=0; t < (BlockDims.z -Gs); t+=jp) {
                vector <int > tempStorage;
                tempStorage.clear();
                tempStorage.push_back(x);
                tempStorage.push_back(y);
                tempStorage.push_back(t);
                // first case is  a 1 to 1 mapping
                // each input grain has an output grain;
                // refinding the index
                tempStorage.push_back(ojp*x/jp);
                tempStorage.push_back(ojp*y/jp);
                tempStorage.push_back(ojp*t/jp);
                GrainScheduler.push_back(tempStorage);
   
            }
        }
    }

}


void GrainRandomizer(){
    vector<vector<int> > Grainshuffled;
    vector<int> tempVec;
    Grainshuffled = GrainScheduler;
    int N = Grainshuffled.size();
    for (int k =0; k <N; k++){
        
        int cual = floor(ofRandom(GrainScheduler.size()));
        tempVec = GrainScheduler[cual];
        Grainshuffled[k][3] = tempVec[3];
        Grainshuffled[k][4] = tempVec[4];
        Grainshuffled[k][5] = tempVec[5];
        GrainScheduler.erase(GrainScheduler.begin()+cual);
    }
    GrainScheduler =Grainshuffled;
    
}

void UpdateSchedulerUpDown(int Upf, int Downf){
    // calculates the grains output position, creates a table of grain input output pairs.
    GrainScheduler.clear();
    
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
    // output jump
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
    
    // First Filling with Upsampling:
    
    for (int y = 0; y < Upf*ojp*(BlockDims.y-Gs)/jp; y+=ojp) {
        for (int x=0; x< Upf*ojp*(BlockDims.x-Gs)/jp; x+=ojp) {
            for (int t=0; t < Upf*ojp*(BlockDims.z-Gs)/jp; t+=ojp) {
                vector <int > tempStorage;
                
                    tempStorage.clear();
                    tempStorage.push_back(jp*(x/ojp)/Upf);
                    tempStorage.push_back(jp*(y/ojp)/Upf);
                if (boolfrozesomegrains) {
                    float chances = ofRandom(100);
                    if(chances>40){
                    tempStorage.push_back(jp*(t/ojp)/Upf);
                    }
                    else{
                        int randplace = floor(ofRandom((BlockDims.z-Gs)/jp))*ojp;
                     tempStorage.push_back(jp*(randplace/ojp)/Upf);
                    }
                }
                else{
                tempStorage.push_back(jp*(t/ojp)/Upf);
                }
                    
                    tempStorage.push_back(x);
                    tempStorage.push_back(y);
                    tempStorage.push_back(t);
                    GrainScheduler.push_back(tempStorage);
                
            }
        }
    }
    
}





float CalculateScalecorrect(void){

    vector<float> Tempovec(10*Gs);
    // initializing
    for (int k=0; k<(10*Gs); k++) {
        Tempovec[k]=0.0;
    }
    
    // output jump
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
    
    // overlaping
    for (int k=0; k<(9*Gs); k+=ojp) {
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

//void testApp::FillGrainList(float GrainSize){
//
//    int Osamp = (int)(ceil(GrainSize*inOp/100.0));
//    int jp = GrainSize - Osamp;
//    vector<Grain> Tgrains;
//    vector<vector<Grain> >XandTgrains;
//    GrainArray.clear();
//    for (int y = 0; y < (BlockDims.y - GrainSize); y+=jp) {
//        XandTgrains.clear();
//        for (int x=0; x< (BlockDims.x - GrainSize); x+=jp) {
//             Tgrains.clear();
//             for (int t=0; t < (BlockDims.z -GrainSize); t+=jp) {
//                 Tgrains.push_back(Grain(x,y,t,GrainSize));
//             }
//             XandTgrains.push_back(Tgrains); 
//        }
//        GrainArray.push_back(XandTgrains);
//    }
//}


//vector<Grain> testApp::ListGrainsThatContribute(float x, float y, float t, float GrainSize){
//
//    vector<Grain> TobeReturned;
//    TobeReturned.clear();
//    int Osamp = (int)(ceil(GrainSize*inOp/100.0));
//    int jp = GrainSize - Osamp;
//    vector<int> Indx;
//    int TemInd = floor(x/(float)jp);
//    Indx.push_back(jp*TemInd);
//    if (((x-jp*TemInd) < Osamp)&&(TemInd!=0)) {
//        Indx.push_back(jp*(TemInd-1));
//    }
//    
//    vector<int> Indy;
//    TemInd = floor(y/(float)jp);
//    Indy.push_back(jp*TemInd);
//    if (((y-jp*TemInd) < Osamp)&&(TemInd!=0)) {
//        Indy.push_back(jp*(TemInd-1));
//    }
//    
//    vector<int> Indt;
//    TemInd = floor(t/(float)jp);
//    Indt.push_back(jp*TemInd);
//    if (((t-jp*TemInd) < Osamp)&&(TemInd!=0)) {
//        Indt.push_back(jp*(TemInd-1));
//    }
//
//    for (int k=0; k<Indx.size(); k++) {
//        for (int q=0; q<Indy.size(); q++) {
//            for (int r=0; r < Indt.size(); r++) {
//                TobeReturned.push_back(Grain(Indx[k],Indy[q],Indt[r],GrainSize));
//            }
//        }
//    }
//    
//    return TobeReturned;
//    
//    
//    
//}









//--------------------------------------------------------------
void testApp::keyPressed(int key){
    
    switch (key) {
        case 'u':
            TheUf++;
            if(TheUf>5){TheUf=1;}
            UpdateBool =true;

            break;
            
        case 'r':
            boolshufflegrains = !boolshufflegrains;
            UpdateBool =true;
            
            break;
        case OF_KEY_RIGHT:
            NewoutOp++;
            UpdateBool = true;
            
            break;
        case OF_KEY_LEFT:
            NewoutOp--;
            UpdateBool = true;
            
            break;
        case 't':
            boolfrozesomegrains = !boolfrozesomegrains;
            UpdateBool = true;
            break;
            
        case OF_KEY_RETURN:
            boolfrozesomegrains = false;
            boolshufflegrains = false;
            TheUf = 1;
            NewoutOp = inOp;
            UpdateBool = true;
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