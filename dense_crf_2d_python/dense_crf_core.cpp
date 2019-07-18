
#include "dense_crf_core.h"

MatrixXf computeUnary( const float *prob, int W, int H){
    const float u_energy = -log( 1.0 / 2 );
    MatrixXf r(2, W*H );
    r.fill(u_energy);
    //printf("%d %d %d \n",im[0],im[1],im[2]);
    for( int k=0; k<W*H; k++ ){
        // Set the energy
        float p = prob[k];
        if(p<0.01)p=0.01;
        if(p>0.99)p=0.99;
        r(0,k) = -log(1.0-p);
        r(1,k) = -log(p);
    }
    return r;
}

VectorXs dense_crf_inference(const unsigned char * img, const float * prob, int W, int H, CRFParam param)
{
    // Number of labels
    const int M = 2;
    MatrixXf unary = computeUnary(prob,W, H);

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
    crf.addPairwiseBilateral(param.alpha, param.alpha, param.beta, param.beta, param.beta, img, new PottsCompatibility( param.w1 ) );
    
    // add a color independent term (feature = pixel location 0..W-1, 0..H-1)
    // x_stddev = 3
    // y_stddev = 3
    // weight = 3
    crf.addPairwiseGaussian( param.gamma, param.gamma, new PottsCompatibility( param.w2 ) );
    
    // Do map inference
    // 	MatrixXf Q = crf.startInference(), t1, t2;
    // 	printf("kl = %f\n", crf.klDivergence(Q) );
    // 	for( int it=0; it<5; it++ ) {
    // 		crf.stepInference( Q, t1, t2 );
    // 		printf("kl = %f\n", crf.klDivergence(Q) );
    // 	}
    // 	VectorXs map = crf.currentMap(Q);
    VectorXs map = crf.map(param.iter);

    return map;
    
//    *res = colorize(map, W,H);
//    return 0;
}