#include "testApp.h"


const int Nx = 352;
const int Ny = 288;
// A vector with all the frames.
vector <cv::Mat> TheFramesInput;

vector<ofVec3f> FirstGrains;
vector<ofVec3f> SecondGrains;
vector<ofVec3f> OutGrains;
vector<ofVec3f> GrainMask;
vector<float> randAng;


const int Nlast = 400;
ofVec3f BlockDims;
float inOp;//  overlaping porcentage
float outOp;//  overlaping porcentage
float NewoutOp;
int Gs = 51; // for fixed grain size in all dimensions
cv::Mat CurrentOutput;
int fc=0; // frame counter


void InitGrains(void);

float CalculateScalecorrect(void);
float ScaleCorrectFactor=1.0;


void GrainRandomizer();

int circuIndex = Nlast-1;


// Control global variables:
float SLD1 = 0;
float SLD2 = 0;
float SLD3 = 0;
float SLD4 = 0;
float SLD5 = 1.0;
float SLD6 = 0.0;
bool boolFixedDelay = true;
bool boolRotGrains = false;
bool boolDelaygrains = false;
bool UpdateVar = true; // variable to update all the info from control and interface.

bool boolhueshift = false;
bool boolSatChange = false;
bool boolGeomGrains = false;
bool boolShiftGrain = false;




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
    
    
   // vidGrabber.setVerbose(true);
 //   vidGrabber.initGrabber(Nx,Ny);
    
    ofSetFrameRate(15);
    inOp =50;
    outOp =50;
    NewoutOp = outOp;

    InitGrains();
//    UpdateSchedulerUpDown(TheUf,1);
    ofSetLogLevel(OF_LOG_VERBOSE);
    vidGrabber1.setDeviceID(0);
    vidGrabber1.initGrabber(Nx,Ny);
  //  vidGrabber2.setDeviceID(1);
   // vidGrabber2.initGrabber(Nx,Ny);
    
    //vidGrabber.initGrabber(Nx,Ny);

//    std::exit(1);
    colorImg.allocate(Nx,Ny);
	grayImage.allocate(Nx,Ny);
    
    
    // Midi
    
    // open port by number
	midiIn.openPort(0);
	//midiIn.openPort("IAC Pure Data In");	// by name
	//midiIn.openVirtualPort("ofxMidiIn Input");	// open a virtual port
	
	// don't ignore sysex, timing, & active sense messages,
	// these are ignored by default
	midiIn.ignoreTypes(false, false, false);
	
	// add testApp as a listener
	midiIn.addListener(this);
	
	// print received messages to the console
	midiIn.setVerbose(false);
    
    
    
    // scale correction  factor 50% overlap
    ScaleCorrectFactor =CalculateScalecorrect();
    
    
}

//--------------------------------------------------------------

// TODO
// change the grain vector to only x,y info
// but two arrays: current grains and next grains
// the t value can be used as offset when ploting
// update cada jp the "next grain info"









void testApp::update(){
    bool bNewFrame = false;
    
    vidGrabber1.update();
    bNewFrame = vidGrabber1.isFrameNew();
    //
    if (bNewFrame){
        //
        // circular buffer
        colorImg.setFromPixels(vidGrabber1.getPixels(), Nx,Ny);
        cv::Mat AuxMat;
        AuxMat = colorImg.getCvImage();
        AuxMat.copyTo(TheFramesInput[circuIndex]);
        fc = circuIndex;
        // Fc has the current sample index
        circuIndex--;
        if (circuIndex<0){
            circuIndex = Nlast-1;
        }
        
     // In case update is needed:
        if(UpdateVar){
        
            // if changes in rotation controls or interface:
            if (boolDelaygrains){
            
                for (int g =0; g < SecondGrains.size(); g++) {
                    if(GrainMask[g].z ==1.0){
                        if (boolFixedDelay){
                            SecondGrains[g].z = (int)((SLD2/127.0)*(Nlast/2.0));
                        
                        }
                        else{
                            SecondGrains[g].z = (int)(ofRandom((SLD2/127.0)*(Nlast/2.0)));
                        }
                    }
                    else{
                        SecondGrains[g].z = 0.0;
                    }
                }
            
            
            }
        
        
        
        
            UpdateVar = false;
        }
        
        
        
        
        
        
        
    // calculate new output size:
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
 
    // input jump size
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;

    
    int outYsize = Gs + ojp*((BlockDims.y-Gs)/jp);
    int outXsize = Gs + ojp*((BlockDims.x-Gs)/jp);
    
   // creating an empty frame for the output

    CurrentOutput =cv::Mat::zeros(outYsize,outXsize, CV_8UC3);
 
    // pointer to output data:
    
    unsigned char *output = (unsigned char*)(CurrentOutput.data);
    
    // running throught the list of grains
    
    for (int g =0; g < FirstGrains.size(); g++) {
        
        
        // two halves
        for (int m=0; m<2; m++) {
          // copy the grain from the input
            int IntiFrame;
            if(m==1){
                    IntiFrame = (fc + (int)FirstGrains[g].z)%Nlast;
             }
            else{
                    IntiFrame = ((fc + (int)SecondGrains[g].z)%Nlast);
            }
            
            //int IntiFrame = (m==1)?(fc + (int)FirstGrains[g].z)%Nlast:((fc + (int)SecondGrains[g].z)%Nlast);
            
            unsigned char *input = (unsigned char*)(TheFramesInput[IntiFrame].data);
            for (int x =0; x<Gs; x++) {
                for (int y=0; y<Gs; y++) {
                    
                    float rIn,gIn,bIn;
                    float rOut,gOut,bOut;
                    float Scale;
                    Scale = ((0.5*0.5*0.5)*(1.0 -cosf(2*PI*x/(float)(Gs-1)))*
                             (1.0 -cosf(2*PI*y/(float)(Gs-1)))*
                             (1.0 -cosf(2*PI*(m*(Gs-1)/2.0 +(fc%jp))/(float)(Gs-1))));
                    // m =0 is the next grain
                    
                    Scale/=(ScaleCorrectFactor*ScaleCorrectFactor*ScaleCorrectFactor);
                    // pixels at the inut
                    int yindIn;
                    int xindIn;
                    
                    if ((!boolRotGrains)||(GrainMask[g].z==0.0)) {
                    
                        if (m==1) {
                            xindIn = x + FirstGrains[g].x;
                            yindIn = y + FirstGrains[g].y;
                        }
                        else {
                            xindIn = x + SecondGrains[g].x;
                            yindIn = y + SecondGrains[g].y;
                        }
                    }
                    if ((boolRotGrains)&&(GrainMask[g].z==1.0)){
                    
                        if (m==1) {
                            float newX = (x + FirstGrains[g].x) - (FirstGrains[g].x+Gs/2.0);
                            float newY = (y + FirstGrains[g].y) - (FirstGrains[g].y+Gs/2.0);
                            xindIn = (int)(newX*cos(PI*SLD1/127.0) - newY*sin(PI*SLD1/127.0)) + (FirstGrains[g].x+Gs/2.0);
                            yindIn = (int)(newX*sin(PI*SLD1/127.0) + newY*cos(PI*SLD1/127.0)) + (FirstGrains[g].y+Gs/2.0);
                        }
                        else {
                            float newX = (x + SecondGrains[g].x) - (SecondGrains[g].x+Gs/2.0);
                            float newY = (y + SecondGrains[g].y) - (SecondGrains[g].y+Gs/2.0);
                            xindIn = (int)(newX*cos(PI*SLD1/127.0) - newY*sin(PI*SLD1/127.0)) + (SecondGrains[g].x+Gs/2.0);
                            yindIn = (int)(newX*sin(PI*SLD1/127.0) + newY*cos(PI*SLD1/127.0)) + (SecondGrains[g].y+Gs/2.0);

                        }
                    
                    
                    }
                    
                    
                    // Geometrical transformations:
                    // center
                    
                    if ((boolGeomGrains)&&(GrainMask[g].z==1.0)){
                        
                        float newX,newY;
                        if (m==1) {
                            newX = (xindIn - FirstGrains[g].x - Gs/2.0);
                            newY = (yindIn - FirstGrains[g].y - Gs/2.0);
                        }
                        
                        else{
                            newX = (xindIn - SecondGrains[g].x - Gs/2.0);
                            newY = (yindIn - SecondGrains[g].y - Gs/2.0);
                        }

                        // First converting to polar and normalizing magnitude
                    
                        float Rf = sqrtf(2*(newX*newX+newY*newY))/Gs;
                        float newR = powf(Rf, 0.6+2*SLD5/127.0);
                    
                        newY = (Rf!=0)?(newR*newY/Rf):0.0;
                        newX = (Rf!=0)?(newR*newX/Rf):0.0;
                    
                        if (m==1) {
                            xindIn = (int)(newX) + (FirstGrains[g].x+Gs/2.0);
                            yindIn = (int)(newY) + (FirstGrains[g].y+Gs/2.0);
                        }
                    
                        else{
                            xindIn = (int)(newX) + (SecondGrains[g].x+Gs/2.0);
                            yindIn = (int)(newY) + (SecondGrains[g].y+Gs/2.0);
                        }
                    
                    }
                    
  
                    
                    xindIn = (xindIn<0)?0:xindIn;
                    xindIn = (xindIn>outXsize-1)?BlockDims.x:xindIn;
                    
                    yindIn = (yindIn<0)?0:yindIn;
                    yindIn = (yindIn>outYsize-1)?BlockDims.y:yindIn;
                    
                    
                    
                    
                    
                    
                    bIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn)) ] ;
                    gIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 1)];
                    rIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 2)];
                    
                    
                    // Hue and saturation transformations:
                    
                    
                    float Ang;
                    if(boolhueshift&&(GrainMask[g].z==1.0)){
                        Ang = SLD3/127.0*2*PI;
                    }
                    else{
                        Ang =0;
                    }
               
                    float S;
                    if(boolSatChange&&(GrainMask[g].z==1.0)){
                      S = 2.1*SLD4/127.0;
                    }
                    else{
                        S = 1.0;
                    }
                    float SU = S*cos(Ang);
                    float SW = S*sin(Ang);
                    float rMed,gMed,bMed;
                    
                  
                    rMed = (.299+.701*SU+.168*SW)*rIn
                    + (.587-.587*SU+.330*SW)*gIn
                    + (.114-.114*SU-.497*SW)*bIn;
                    gMed = (.299-.299*SU-.328*SW)*rIn
                    + (.587+.413*SU+.035*SW)*gIn
                    + (.114-.114*SU+.292*SW)*bIn;
                    bMed = (.299-.3*SU+1.25*SW)*rIn
                    + (.587-.588*SU-1.05*SW)*gIn
                    + (.114+.886*SU-.203*SW)*bIn;
                   
                    
                    // from:
                    // http://beesbuzz.biz/code/hsv_color_transforms.php
                    
                    
                    
                    
                    
                    
                    
                    //pixels at the output
                    // output indexes
                    
                    int xindOut;
                    int yindOut;
                    xindOut = x + OutGrains[g].x ;
                    yindOut = y + OutGrains[g].y;
                    if((boolShiftGrain)&&(GrainMask[g].z==1.0)){
                        xindOut = x + OutGrains[g].x + SLD6/127.0*Gs*cos(randAng[g]);
                        yindOut = y + OutGrains[g].y + SLD6/127.0*Gs*sin(randAng[g]);
                    }

                  
                    xindOut = (xindOut<0)?0:xindOut;
                    xindOut = (xindOut>outXsize-1)?outXsize-1:xindOut;
                    
                    yindOut = (yindOut<0)?0:yindOut;
                    yindOut = (yindOut>outYsize-1)?outYsize-1:yindOut;
                    
                    
                    
                    bOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] ;
                    gOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1];
                    rOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2];
                    

                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] = (uchar)(bOut + Scale*bMed);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1]=(uchar)(gOut +Scale*gMed);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2]=(uchar)(rOut + Scale*rMed);
                    
                    
                
                
                } //y
            } //x
            
        }//m
        
        if ((fc%jp) == jp-1) {
            FirstGrains[g] = SecondGrains[g];
        }
        
        
        
    }//g
    
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
   

    
    // Ploting the  mask
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    
    for (int g =0; g < FirstGrains.size(); g++) {
        ofSetColor(0, 0, 0,100);
        ofNoFill();
        ofCircle(2*(FirstGrains[g].x +Gs/2.0) + 2*Nx, 2*(FirstGrains[g].y +Gs/2.0), 2*Gs/2.0);
        ofSetColor(60 + 160*(GrainMask[g].z), 0, 60,100);
        ofFill();
        ofCircle(2*(FirstGrains[g].x +Gs/2.0) + 2*Nx, 2*(FirstGrains[g].y +Gs/2.0), 2*Gs/2.0);
    }
    
    glDisable(GL_BLEND);

    
    fc++;
    if(fc>BlockDims.z-1){fc=0;}
    
}


void testApp::exit(){

  
    vidGrabber1.close();
   // vidGrabber2.close();


}

void InitGrains(){
// calculates the grains input and output positions.
    FirstGrains.clear();
    SecondGrains.clear();
    OutGrains.clear();
    
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
    // output jump
    int outOsamp = (int)(ceil(Gs*outOp/100.0));
    int ojp = Gs - outOsamp;
    
    
    for (int y = 0; y < (BlockDims.y - Gs); y+=jp) {
        for (int x=0; x< (BlockDims.x - Gs); x+=jp) {
                ofVec3f tempStorage;
                tempStorage.x = x;
                tempStorage.y=y;
                tempStorage.z=0;
                FirstGrains.push_back(tempStorage);
                SecondGrains.push_back(tempStorage);
                OutGrains.push_back(tempStorage);
                GrainMask.push_back(tempStorage);
                randAng.push_back(ofRandom(2*PI));
        }
    }

}


void GrainRandomizer(){
    for (int g =0; g < SecondGrains.size(); g++) {
        
        if(GrainMask[g].z ==1.0){
            SecondGrains[g].z = (int)ofRandom(Nlast/4.0);
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



//--------------------------------------------------------------
void testApp::keyPressed(int key){
    
    switch (key) {
        case 'u':


            break;
            
        case 'r':
            boolDelaygrains = !boolDelaygrains;

            
            break;
        case OF_KEY_RIGHT:
            NewoutOp++;

            
            break;
        case OF_KEY_LEFT:
            NewoutOp--;
            
            break;
        case 't':
 

            break;
            
        case OF_KEY_RETURN:
   
            boolDelaygrains = false;
            NewoutOp = inOp;

            break;
        case 'q':
            boolRotGrains =!boolRotGrains;
            break;
        default:
            break;
    }
    
    
}


//
void testApp::newMidiMessage(ofxMidiMessage& msg) {
    
	
	switch(msg.control){
            // slidders:
            case 0: // slidder 1 rotation
                SLD1 = msg .value;
            break;
            case 1: // slidder 2 time delay
                SLD2 = msg .value;
                UpdateVar =true;
            break;
            case 2: // slidder 3 hue
            SLD3 = msg .value;
            break;
            case 3: // slidder 4 sat
            SLD4 = msg .value;
            break;
            case 4: // slidder 4 sat
            SLD5 = msg .value;
            break;
            case 5:
            SLD6 = msg .value;
            break;
            // Buttons:
            
            // Select none
            case 45:
            for (int g=0; g < GrainMask.size(); g++) {
                GrainMask[g].z = 0.0;
            }
            UpdateVar =true;
            
            
            break;
            
            // Select all
        case 41:
            for (int g=0; g < GrainMask.size(); g++) {
               GrainMask[g].z = 1.0;
            }
            UpdateVar =true;
            
            
            break;
            
            
            
            case 59:
                if(msg.value >100){
                    boolFixedDelay = !boolFixedDelay;
                    UpdateVar = true;
                }
            break;
            case 64:
            boolRotGrains = true;
            break;
            case 48:
            boolRotGrains = false;
            break;
            case 65:
            if(msg.value >100){
                boolDelaygrains = true;
                UpdateVar =true;
            }
            break;
            
        case 50:
            boolhueshift = false;
            break;
        case 66:
            if(msg.value >100){
                boolhueshift = true;
            }
            break;
        case 51:
            boolSatChange = false;
            break;
        case 67:
            if(msg.value >100){
                boolSatChange = true;
            }
            break;
        case 53:
            boolShiftGrain = false;
            break;
        case 69:
            if(msg.value >100){
                boolShiftGrain = true;
            }
            break;
            
        case 52:
            boolGeomGrains = false;
            break;
        case 68:
            if(msg.value >100){
                boolGeomGrains = true;
            }
            break;
            
            case 49:
            if(msg.value >100){
                if(boolDelaygrains){
                    // reset delays
                    for (int g =0; g < SecondGrains.size(); g++) {
                        SecondGrains[g].z = 0.0;
                    }
                }
                boolDelaygrains = false;
            }
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
    
    for (int g=0; g < GrainMask.size(); g++) {
        
        float distX = 2*(GrainMask[g].x +Gs/2.0) + 2*Nx - x;
        float distY = 2*(GrainMask[g].y +Gs/2.0)-y;
        if((distX*distX + distY*distY) < 30*30){
            GrainMask[g].z = (GrainMask[g].z==0);
        }
        
        
    }
    UpdateVar =true;

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