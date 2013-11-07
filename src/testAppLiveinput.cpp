#include "testApp.h"

const int Nx = 352;
const int Ny = 288;




// A vector with all the frames.
vector <cv::Mat> TheFramesInput;

//vector<Grain> GrainList;
const int Nlast = 300;
ofVec3f BlockDims;

float inOp;//  overlaping porcentage
int Gs = 31; // for fixed grain size in all dimensions
cv::Mat CurrentOutput;
int fc=0; // frame counter
int ofc =0; // output frame counter

vector<ofVec3f > GrainScheduler;
void UpdateScheduler(void);

int OutFramesN;
bool UpdateBool;
float CalculateScalecorrect(void);
float ScaleCorrectFactor=1.0;
ofVec3f GrainCamera;
ofVec3f GrainUpVector;
int outYsize;
int outXsize;
 float minX,maxX,minY,maxY,minZ,maxZ;
ofVec3f GrainViewTransform(ofVec3f InPutPoint, ofVec3f Eye, ofVec3f Target, ofVec3f UpVec);
ofVec3f RotateAux(ofVec3f input, float latitud, float azimuth, ofVec3f pivot);
int viewerDir =0;
int circuIndex = 0;

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
    UpdateScheduler();
    //UpdateSchedulerUpDown(TheUf,1);
//    GrainCamera =  RotateAux(ofVec3f(0.0,0.0,-1.1* sqrt(BlockDims.x*BlockDims.x+BlockDims.y*BlockDims.y+
//                                                        BlockDims.z*BlockDims.z)/2.0), 85.0, 0.0,
//                             ofVec3f(BlockDims.x/2.0,BlockDims.y/2.0,BlockDims.z/2.0));
  
 //   GrainUpVector = ofVec3f(0.0,-1.0,0.0);
    //side
  //  GrainCamera = ofVec3f(1.2*BlockDims.x,BlockDims.y/2.0,BlockDims.z/2.0);
    //vertex
   // GrainCamera = ofVec3f(BlockDims.x,0.0,0.0);
// corner
 //    GrainCamera = ofVec3f(BlockDims.x,BlockDims.y/2.0,0.0);
  // corner side
   //  GrainCamera = ofVec3f(BlockDims.x,0.0,BlockDims.z/2);
    // up
   // GrainCamera = ofVec3f(BlockDims.x/2.0,0.0,BlockDims.z/2);
    //GrainUpVector = ofVec3f(-1.0,0.0,0.0);
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
    

    if (UpdateBool){
        
        switch (viewerDir) {
                case 0: // front
                GrainUpVector = ofVec3f(0.0,-1.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x/2.0,BlockDims.y/2.0,- BlockDims.z/2.0);
                break;
                
                case 1: // side
                GrainUpVector = ofVec3f(0.0,-1.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x,BlockDims.y/2.0,BlockDims.z/2.0);
                break;
                
                case 2: // vertex
                GrainUpVector = ofVec3f(0.0,-1.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x,0.0,0.0);
                break;
                
                case 3: // corner
                GrainUpVector = ofVec3f(0.0,-1.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x,BlockDims.y/2.0,0.0);
                break;
                
                case 4: // corner side
                GrainUpVector = ofVec3f(0.0,-1.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x,0.0,BlockDims.z/2);
                break;
                
                case 5: // up
                GrainUpVector = ofVec3f(-1.0,0.0,0.0);
                GrainCamera = ofVec3f(BlockDims.x/2.0,0.0,BlockDims.z/2);
                break;
                
                
            default:
                break;
        }
        
        
        
        
        // input jump size
        int Osamp = (int)(ceil(Gs*inOp/100.0));
        int jp = Gs - Osamp;

         // calculate the new ouptput size and video length;
    
         ofVec3f TheExtremes;
    
        // checking the coordenates of the 8 vertex of the video.
    
    
        minX =10000000; maxX = -1000000; minY = 1000000; maxY = -1000000; minZ = 1000000; maxZ=-100000;
        for (int k =0; k<8; k++) {
            ofVec3f TheInputVec(floor(k/4)*BlockDims.x,(((int)(floor(k/2))%2)*BlockDims.y),(k%2)*BlockDims.z);
            
            
            
            
            TheExtremes  = GrainViewTransform(TheInputVec,GrainCamera,
                ofVec3f(BlockDims.x/2.0,BlockDims.y/2.0,BlockDims.z/2.0),
                GrainUpVector);
//             cout<<TheInputVec.x<<","<<TheInputVec.y<<","<<TheInputVec.z<<endl;
//            cout<<TheExtremes.x<<","<<TheExtremes.y<<","<<TheExtremes.z<<endl;
   
            if (TheExtremes.x<minX) {
                minX = TheExtremes.x;
            }
            if (TheExtremes.x>maxX) {
                maxX = TheExtremes.x;
            }
            if (TheExtremes.y<minY) {
                minY = TheExtremes.y;
            }
            if (TheExtremes.y>maxY) {
                maxY = TheExtremes.y;
            }
            if (TheExtremes.z<minZ) {
                minZ = TheExtremes.z;
            }
            if (TheExtremes.z>maxZ) {
                maxZ = TheExtremes.z;
            }
 //       cout <<minX<<","<<maxX<<","<<minY<<","<<maxY<<","<<minZ<<","<<maxZ<<endl;
    
        }
        
    //   cout <<minX<<","<<maxX<<","<<minY<<","<<maxY<<","<<minZ<<","<<maxZ<<endl;
        outYsize = ceil(maxY-minY);
        outXsize = ceil(maxX-minX);
        // calculate the number of output frames.
        OutFramesN = ceil(maxZ-minZ);
    
      //  if(ofc>OutFramesN){ofc = floor(minZ);}
    ofc = floor(minZ);
        // creating an empty frame for the output

        
 
        UpdateBool = false;
    
    //cout<<outYsize<<endl;
    }
    CurrentOutput =cv::Mat::zeros(outYsize,outXsize, CV_8UC3);
    // pointer to output data:
//    cout<<outYsize<<endl;
//    cout<<outXsize<<endl;
    unsigned char *output = (unsigned char*)(CurrentOutput.data);
    
    // running throught the list of grains
//    cout<<"OFC:"<<ofc<<endl;
//    cout<<"GRAIN SIZE:"<<GrainScheduler.size()<<endl;
    
    for (int g =0; g < GrainScheduler.size(); g++) {
        
        // check if in the current frame the current grain needs to be ploted
        
       
       ofVec3f TransGrain =
        GrainViewTransform(ofVec3f(GrainScheduler[g].x, GrainScheduler[g].y,GrainScheduler[g].z),
                           GrainCamera, ofVec3f(BlockDims.x/2.0,BlockDims.y/2.0,BlockDims.z/2.0),
                           GrainUpVector);
        
        
        int indT = ofc - TransGrain.z;
        if ((indT >=0)&&(indT<Gs)) {
          // copy the grain from the input
          // getting a pointer to the input frame
//            cout<<"TransGrainX:"<<TransGrain.x<< "minX:"<<minX<<endl;
//            cout<<"TransGrainX- minX:"<<TransGrain.x - minX<<endl;
//            cout<<"floor(TransGrainX-minX)"<<floor(TransGrain.x-minX)<<endl;
//            cout<<"X:"<<GrainScheduler[g].x<<"Y:"<<GrainScheduler[g].y<<"Z:"<<GrainScheduler[g].z<<endl;
//            cout<<"outX:"<<floor(TransGrain.x-minX)<<"outY:"<<floor(TransGrain.y - minY)<<"outZ:"<<ofc - TransGrain.z<<endl;
            
             unsigned char *input = (unsigned char*)(TheFramesInput[GrainScheduler[g].z+indT].data);
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
                    int yindIn = y + GrainScheduler[g].y;
                    int xindIn = x + GrainScheduler[g].x;
                    
                  
                  
                    bIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn)) ] ;
                    gIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 1)];
                    rIn = (float)input[(int)(3*BlockDims.x * (yindIn) + 3*(xindIn) + 2)];
                    //pixels at the output
                    // output indexes
                    int xindOut = x + floor(TransGrain.x-minX);
                    xindOut = (xindOut<0)?0:xindOut;
                    xindOut = (xindOut>outXsize-1)?outXsize-1:xindOut;
                    int yindOut = y + floor(TransGrain.y - minY);
                    yindOut = (yindOut<0)?0:yindOut;
                    yindOut = (yindOut>outYsize-1)?outYsize-1:yindOut;
                    
                    bOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] ;
                    gOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1];
                    rOut = (float)output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2];
                    

                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) ] = (uchar)(bOut + Scale*bIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 1]=(uchar)(gOut +Scale*gIn);
                    output[3*CurrentOutput.cols * (yindOut) + 3*(xindOut) + 2]=(uchar)(rOut + Scale*rIn);
                    
                    
                
                
                }// for y
            } // for x
         } // if indT
    }// for g
    

    
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
    
    if(ofc>OutFramesN){ofc = floor(minZ);}
    fc++;
    if(fc>BlockDims.z-1){fc=0;}
    
}


ofVec3f GrainViewTransform(ofVec3f InPutPoint, ofVec3f Eye, ofVec3f Target, ofVec3f UpVec){

// Transform the cordinates of a point InPutPoint. For a camera placed at Eye and looking at target with Up vector UpVec
    
    ofMatrix4x4 ViewMat;
    ofVec3f zaxis(Target - Eye);
    zaxis.normalize();
    ofVec3f xaxis = UpVec.getCrossed(zaxis);
    xaxis.normalize();
    ofVec3f yaxis = zaxis.getCrossed(xaxis);
    // translate vect
    
    InPutPoint -= Eye;
    
    // project on each axis
    
    ofVec3f Outpos;
    
    Outpos.x = -xaxis.dot(InPutPoint);
    Outpos.y = -yaxis.dot(InPutPoint);
    Outpos.z = zaxis.dot(InPutPoint);
    
    return Outpos;
    
}


ofVec3f RotateAux(ofVec3f input, float latitud, float azimuth, ofVec3f pivot){


    ofVec3f theOut = input.getRotated(azimuth, ofVec3f(1, 0, 0));
    
    theOut = theOut.getRotated(latitud, ofVec3f(0, 1, 0));
    
    theOut += pivot;

    return theOut;


}


void UpdateScheduler(){
// calculates the grains output position, creates a table of grain input output pairs.
    GrainScheduler.clear();
    
    // input jump
    int Osamp = (int)(ceil(Gs*inOp/100.0));
    int jp = Gs - Osamp;
    
       
    
    for (int y = 0; y < (BlockDims.y - Gs); y+=jp) {
        for (int x=0; x< (BlockDims.x - Gs); x+=jp) {
            for (int t=0; t < (BlockDims.z -Gs); t+=jp) {
                ofVec3f tempStorage;               
                tempStorage.x =x;
                tempStorage.y = y;
                tempStorage.z =t;
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
        case 'u':

            break;
            
        case 's':
      
            break;
        case OF_KEY_RIGHT:
            viewerDir++;
            if (viewerDir >5) {
                viewerDir=0;
            }
            UpdateBool =true;
            break;
        case OF_KEY_LEFT:
            viewerDir--;
            if (viewerDir<0) {
                viewerDir =5;
            }
            UpdateBool =true;
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