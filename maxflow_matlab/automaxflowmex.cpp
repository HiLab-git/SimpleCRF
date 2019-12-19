//------------------------------------------------------------------------
// guotai wang
// guotai.wang.14@ucl.ac.uk
// 11 Jan, 2016
// max flow with possibility input
//------------------------------------------------------------------------

#include "mex.h"
#include "maxflow-v3.0/graph.h"
#include <iostream>
#include <cmath>

using namespace std;

// [flow label]=automaxflowmex(I,fgProb, bgProb,lambda,sigma)
void mexFunction(int			nlhs, 		/* number of expected outputs */
                 mxArray		*plhs[],	/* mxArray output pointer array */
                 int			nrhs, 		/* number of inputs */
                 const mxArray	*prhs[]		/* mxArray input pointer array */)
{
    // input checks
    if (nrhs != 5 )
    {
        mexErrMsgTxt ("USAGE: [flow label]=automaxflowmex(I,fgProb, bgProb,lambda,sigma);");
    }
    const mxArray *I = prhs[0];
    const mxArray *fgProb = prhs[1];
    const mxArray *bgProb = prhs[2];
    double lamda= * mxGetPr(prhs[3]);
    double sigma= * mxGetPr(prhs[4]);
    
    unsigned char * IPr=(unsigned char *)mxGetPr(I);
    double * fgProbPr=mxGetPr(fgProb);
    double * bgProbPr=mxGetPr(bgProb);
    // size of image
    mwSize m = mxGetM(I);//height
    mwSize n = mxGetN(I);//width
    
    //construct graph
    typedef Graph<float,float,float> GraphType;
    GraphType *g = new GraphType(/*estimated # of nodes*/ m*n, /*estimated # of edges*/ 2*m*n);
    g->add_node(m*n);
    
    float maxWeight=-10000;
    for(int x=0;x<n;x++)
    {
        for(int y=0;y<m;y++)
        {
            //n-link
            float pValue=(float)*(IPr+x*m+y);
            int uperPointx=x;
            int uperPointy=y-1;
            int LeftPointx=x-1;
            int LeftPointy=y;
            float n_weight=0;
            if(uperPointy>=0 && uperPointy<m)
            {
                float qValue=(float)*(IPr+uperPointx*m+uperPointy);
                n_weight=lamda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*m+y;
                int qIndex=uperPointx*m+uperPointy;
                
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(n_weight>maxWeight)
            {
                maxWeight=n_weight;
            }
            
            if(LeftPointx>=0 && LeftPointx<n)
            {
                float qValue=(float)*(IPr+LeftPointx*m+LeftPointy);
                n_weight=lamda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*m+y;
                int qIndex=LeftPointx*m+LeftPointy;
                
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(n_weight>maxWeight)
            {
                maxWeight=n_weight;
            }
        }
    }
    
    for(int x=0;x<n;x++)
    {
        for(int y=0;y<m;y++)
        {
            float forePosibility = (float)*(fgProbPr+x*m+y);
            if(forePosibility<0.001) forePosibility=0.001;
            float backPosibility = (float)*(bgProbPr+x*m+y);
            if(backPosibility<0.001) backPosibility=0.001;
            float s_weight=-log(backPosibility);
            float t_weight=-log(forePosibility);
            
            int pIndex=x*m+y;
            g->add_tweights(pIndex,s_weight,t_weight);
        }
    }
    // return the results
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* flow = mxGetPr(plhs[0]);
    *flow = g->maxflow();
    //printf("max flow: %f\n",*flow);
    // figure out segmentation
    plhs[1] = mxCreateNumericMatrix(m, n, mxUINT8_CLASS, mxREAL);
    unsigned char * labels = (unsigned char*)mxGetData(plhs[1]);
    for (int x = 0; x < n; x++)
    {
        for (int y=0;y<m;y++)
        {
            int Index=x*m+y;
            labels[Index] = g->what_segment(Index);
        }
    }
    // cleanup
    delete g;
}

