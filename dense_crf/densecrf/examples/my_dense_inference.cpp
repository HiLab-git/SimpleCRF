/*
 Copyright (c) 2013, Philipp Krähenbühl
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of the Stanford University nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY Philipp Krähenbühl ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL Philipp Krähenbühl BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "densecrf.h"
#include <cstdio>
#include <stdio.h>
#include <cmath>
#include "lodepng.h"
#include "ppm.h"
#include "common.h"
#include <iostream>
#include <fstream>
using namespace std;
// Certainty that the groundtruth is correct
const float GT_PROB = 0.5;

// Simple classifier that is 50% certain that the annotation is correct
MatrixXf computeUnary( const VectorXs & lbl, int M ){
    const float u_energy = -log( 1.0 / M );
    const float n_energy = -log( (1.0 - GT_PROB) / (M-1) );
    const float p_energy = -log( GT_PROB );
    MatrixXf r( M, lbl.rows() );
    r.fill(u_energy);
    //printf("%d %d %d \n",im[0],im[1],im[2]);
    for( int k=0; k<lbl.rows(); k++ ){
        // Set the energy
        if (lbl[k]>=0){
            r.col(k).fill( n_energy );
            r(lbl[k],k) = p_energy;
        }
    }
    return r;
}

MatrixXf computeUnary( unsigned char *prob, int W, int H){
    const float u_energy = -log( 1.0 / 2 );
//    const float n_energy = -log( (1.0 - GT_PROB) / (M-1) );
//    const float p_energy = -log( GT_PROB );
    MatrixXf r( 2, W*H );
    r.fill(u_energy);
    //printf("%d %d %d \n",im[0],im[1],im[2]);
    for( int k=0; k<W*H; k++ ){
        // Set the energy
        float p = prob[k]/255.0;
        r(0,k) = -log(1-p);
        r(1,k) = -log(p);
//        if (lbl[k]>=0){
//            r.col(k).fill( n_energy );
//            r(lbl[k],k) = p_energy;
//        }
    }
    return r;
}

unsigned char * readPNG( char * infile, int& W, int&H, LodePNGColorType format)
{
    std::vector<unsigned char> image; //the raw pixels
    unsigned width, height;
    
    unsigned error = lodepng::decode(image, width, height, infile, format, 8);
    
    if(error)
    {
        std::cout << "load png error " << std::endl;
        return 0;
    }
    W = width;
    H = height;
    std::cout<<"png image loaded, size is "<<W<<" "<<H<<" "<<image.size()<<std::endl;
    std::cout<<"min and max "<<W<<" "<<H<<" "<<image.size()<<std::endl;
    int imsize;
    if(format==LCT_GREY){
        imsize = W*H;
    }
    else{
        imsize = W*H*3;
    }
    unsigned char * r = new unsigned char[imsize];
    memcpy(r, &image[0], sizeof(unsigned char)*imsize);
    return r;
}

int main( int argc, char* argv[]){
    if (argc<5){
        printf("Usage: %s param image probability output\n", argv[0] );
        return 1;
    }
    float w1, alpha, beta, w2, gamma;
    int iter;
    // read parameters from file
    ifstream paramfile(argv[1]);
    if(paramfile.is_open())
    {
        paramfile>>w1>>alpha>>beta>>w2>>gamma>>iter;
        paramfile.close();
    }
    else{
        cout<<"unable to open param file: "<<argv[1]<<endl;
        return 0;
    }
    cout<<"params"<<endl;
    cout<<"w1    "<<w1   <<endl;
    cout<<"alpha "<<alpha<<endl;
    cout<<"beta  "<<beta <<endl;
    cout<<"w2    "<<w2   <<endl;
    cout<<"gamma "<<gamma<<endl;
    cout<<"iter  "<<iter<<endl;
    
    // Number of labels
    const int M = 2;
    // Load the color image and some crude annotations (which are used in a simple classifier)
    int W, H, GW, GH;
    unsigned char * im = readPNG(argv[2], W, H, LCT_RGB);
    if (!im){
        printf("Failed to load image!\n");
        return 1;
    }
    unsigned char * anno = readPNG( argv[3], GW, GH, LCT_GREY );
    if (!anno){
        printf("Failed to load annotations!\n");
        return 1;
    }
    if (W!=GW || H!=GH){
        printf("Annotation size doesn't match image!\n");
        return 1;
    }

    /////////// Put your own unary classifier here! ///////////
//    MatrixXf unary = computeUnary( getLabeling( anno, W*H, M ), M );
    MatrixXf unary = computeUnary(anno,W, H);
    ///////////////////////////////////////////////////////////
    
    // Setup the CRF model
    DenseCRF2D crf(W, H, M);
    // Specify the unary potential as an array of size W*H*(#classes)
    // packing order: x0y0l0 x0y0l1 x0y0l2 .. x1y0l0 x1y0l1 ...
    crf.setUnaryEnergy( unary );

    // add a color dependent term (feature = xyrgb)
    // x_stddev = 60
    // y_stddev = 60
    // r_stddev = g_stddev = b_stddev = 20
    // weight = 10
    crf.addPairwiseBilateral( alpha, alpha, beta, beta, beta, im, new PottsCompatibility( w1 ) );
    
    // add a color independent term (feature = pixel location 0..W-1, 0..H-1)
    // x_stddev = 3
    // y_stddev = 3
    // weight = 3
    crf.addPairwiseGaussian( gamma, gamma, new PottsCompatibility( w2 ) );
    
    // Do map inference
    // 	MatrixXf Q = crf.startInference(), t1, t2;
    // 	printf("kl = %f\n", crf.klDivergence(Q) );
    // 	for( int it=0; it<5; it++ ) {
    // 		crf.stepInference( Q, t1, t2 );
    // 		printf("kl = %f\n", crf.klDivergence(Q) );
    // 	}
    // 	VectorXs map = crf.currentMap(Q);
    VectorXs map = crf.map(iter);
    // Store the result
    unsigned char *res = colorize( map, W, H );
    
//    std::vector<unsigned char> resv(res, res+W*H*3);
//    std::cout<<"result size "<<resv.size()<<std::endl;
//    writePPM( argv[3], W, H, res );
    std::vector<unsigned char> png;
    png.resize(W*H*4);
    for(int i=0; i<W*H; i++)
    {
        
        png[4*i+0]=map(i);
        png[4*i+1]=map(i);
        png[4*i+2]=map(i);
        png[4*i+3]=255;
    }

    unsigned error = lodepng::encode(argv[4], png, W, H);
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
    
    delete[] im;
    delete[] anno;
    delete[] res;
}
