#include <Python.h>
#include "numpy/arrayobject.h"
//#include "Numeric/arrayobject.h"
#include "maxflow-v3.0/graph.h"
#include <iostream>
using namespace std;

// example to use numpy object: http://blog.debao.me/2013/04/my-first-c-extension-to-numpy/
// write a c extension ot Numpy: http://folk.uio.no/hpl/scripting/doc/python/NumPy/Numeric/numpy-13.html
static PyObject *
max_flow_wrapper(PyObject *self, PyObject *args)
{
    PyObject *I=NULL, *fP=NULL, *bP=NULL, *param=NULL;
    PyArrayObject *arr_I=NULL, *arr_fP=NULL, *arr_bP=NULL;
    
    if (!PyArg_ParseTuple(args, "OOOO", &I, &fP, &bP, &param)) return NULL;
    
    arr_I = (PyArrayObject*)PyArray_FROM_OTF(I, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_I == NULL) return NULL;
    
    arr_fP = (PyArrayObject*)PyArray_FROM_OTF(fP, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_fP == NULL) return NULL;

    arr_bP = (PyArrayObject*)PyArray_FROM_OTF(bP, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_bP == NULL) return NULL;

    
    /*vv* code that makes use of arguments *vv*/
    
    int nd = PyArray_NDIM(arr_I);   //number of dimensions
    npy_intp * shape = PyArray_DIMS(arr_I);  // npy_intp array of length nd showing length in each dim.
    npy_intp * shape_fP = PyArray_DIMS(arr_fP);
    npy_intp * shape_bP = PyArray_DIMS(arr_bP);
//    cout<<"crf input shape: ";
    for(int i=0; i<nd; i++)
    {
//        cout<<shape[i]<<" ";
        if(shape[i] !=shape_fP[i] || shape[i]!=shape_bP[i])
        {
            cout<<"input shape does not match"<<endl;
            return NULL;
        }
    }
//    cout<<endl;
//    for (int i=0; i<nd; ++i)
//    std::cout<<" "<<shape[i];
//    std::cout<<std::endl;
    
//    std::cout<<"stride "<<arr_I->strides[0]<<" "<<arr_I->strides[1]<<endl;
    
//    for (int i=0; i<arr_fP->nd; ++i)
//    std::cout<<" "<<arr_fP->dimensions[i];
//    std::cout<<std::endl;
    
//    Py_ssize_t len=PyTuple_GET_SIZE(param);
//    cout<<"param size "<<len<<std::endl;
//    for(int i=0;i<len;i++)
//    {
//        PyObject * temp_object =PyTuple_GET_ITEM(param, i);
//        cout<<PyFloat_AsDouble(temp_object)<<" ";
//    }
    double lambda = PyFloat_AsDouble(PyTuple_GET_ITEM(param, 0));
    double sigma  = PyFloat_AsDouble(PyTuple_GET_ITEM(param, 1));
//    cout<<"lambda "<<lambda<<", sigma "<<sigma<<endl;
    
    typedef Graph<float, float, float> GraphType;
//    /*estimated # of nodes*/ /*estimated # of edges*/
    GraphType * g = new GraphType(shape[0]*shape[1], 2*shape[0]*shape[1]);
    g->add_node(shape[0]*shape[1]);
    float max_weight = -100000;
    for(int x=0; x<shape[0];x++)
    {
        for(int y=0; y<shape[1];y++)
        {
            float pValue = *(float *)(arr_I->data + x*arr_I->strides[0] + y*arr_I->strides[1]);
            int uperPointx=x;
            int uperPointy=y-1;
            int LeftPointx=x-1;
            int LeftPointy=y;
            float n_weight=0;
            if(uperPointy>=0 && uperPointy<shape[1])
            {
                float qValue=*(float *)(arr_I->data + uperPointx*arr_I->strides[0] + uperPointy*arr_I->strides[1]);
                n_weight=lambda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*shape[1]+y;
                int qIndex=uperPointx*shape[1]+uperPointy;
                
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(LeftPointx>=0 && LeftPointx<shape[0])
            {
                float qValue=*(float *)(arr_I->data + LeftPointx*arr_I->strides[0] + LeftPointy*arr_I->strides[1]);
                n_weight=lambda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*shape[1]+y;
                int qIndex=LeftPointx*shape[1]+LeftPointy;
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(n_weight>max_weight)
            {
                max_weight=n_weight;
            }
        }
    }

    for(int x=0;x<shape[0];x++)
    {
        for(int y=0;y<shape[1];y++)
        {
            float forePosibility=*(float *)(arr_fP->data + x*arr_fP->strides[0] + y*arr_fP->strides[1]);
            float backPosibility=*(float *)(arr_bP->data + x*arr_bP->strides[0] + y*arr_bP->strides[1]);
            float s_weight=-log(backPosibility);
            float t_weight=-log(forePosibility);
            int pIndex=x*shape[1]+y;
            g->add_tweights(pIndex,s_weight,t_weight);
        }
    }
    double flow;
    flow = g->maxflow();
//    printf("max flow: %f\n",flow);

    int outshape[2];
    outshape[0]=shape[0];
    outshape[1]=shape[1];
    PyArrayObject * labels = (PyArrayObject*)  PyArray_FromDims(2, outshape, NPY_INT8);
    for(int x=0; x<shape[0]; x++)
    {
        for(int y=0; y<shape[1]; y++)
        {
            int Index = x*shape[1]+y;
            *(labels->data + x*labels->strides[0] + y*labels->strides[1]) = 1-g->what_segment(Index);
        }
    }
    delete g;

    Py_DECREF(arr_I);
    Py_DECREF(arr_fP);
    Py_DECREF(arr_bP);
    Py_INCREF(labels);
    return PyArray_Return(labels);
}

static PyObject *
interactive_max_flow_wrapper(PyObject *self, PyObject *args)
{
    PyObject *I=NULL, *fP=NULL, *bP=NULL, *Seed=NULL, *param=NULL;
    PyArrayObject *arr_I=NULL,  *arr_fP=NULL, *arr_bP=NULL, *arr_Seed=NULL;
    
    if (!PyArg_ParseTuple(args, "OOOOO", &I, &fP, &bP, &Seed, &param)) return NULL;
    
    arr_I = (PyArrayObject*)PyArray_FROM_OTF(I, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_I == NULL) return NULL;
    
    arr_fP = (PyArrayObject*)PyArray_FROM_OTF(fP, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_fP == NULL) return NULL;
    
    arr_bP = (PyArrayObject*)PyArray_FROM_OTF(bP, NPY_FLOAT32, NPY_IN_ARRAY);
    if (arr_bP == NULL) return NULL;
    
    arr_Seed = (PyArrayObject*)PyArray_FROM_OTF(Seed, NPY_UINT8, NPY_IN_ARRAY);
    if (arr_Seed == NULL) return NULL;
    
    /*vv* code that makes use of arguments *vv*/
    
    int nd = PyArray_NDIM(arr_I);   //number of dimensions
    npy_intp * shape = PyArray_DIMS(arr_I);  // npy_intp array of length nd showing length in each dim.
    npy_intp * shape_fP = PyArray_DIMS(arr_fP);
    npy_intp * shape_bP = PyArray_DIMS(arr_bP);
    npy_intp * shape_seed = PyArray_DIMS(arr_Seed);
    //    cout<<"crf input shape: ";
    for(int i=0; i<nd; i++)
    {
        //        cout<<shape[i]<<" ";
        if(shape[i] !=shape_fP[i] || shape[i]!=shape_bP[i] || shape[i]!=shape_seed[i])
        {
            cout<<"input shape does not match"<<endl;
            return NULL;
        }
    }

    double lambda = PyFloat_AsDouble(PyTuple_GET_ITEM(param, 0));
    double sigma  = PyFloat_AsDouble(PyTuple_GET_ITEM(param, 1));
    //    cout<<"lambda "<<lambda<<", sigma "<<sigma<<endl;
    
    typedef Graph<float, float, float> GraphType;
    //    /*estimated # of nodes*/ /*estimated # of edges*/
    GraphType * g = new GraphType(shape[0]*shape[1], 2*shape[0]*shape[1]);
    g->add_node(shape[0]*shape[1]);
    float max_weight = -100000;
    for(int x=0; x<shape[0];x++)
    {
        for(int y=0; y<shape[1];y++)
        {
            float pValue = *(float *)(arr_I->data + x*arr_I->strides[0] + y*arr_I->strides[1]);
            int uperPointx=x;
            int uperPointy=y-1;
            int LeftPointx=x-1;
            int LeftPointy=y;
            float n_weight=0;
            if(uperPointy>=0 && uperPointy<shape[1])
            {
                float qValue=*(float *)(arr_I->data + uperPointx*arr_I->strides[0] + uperPointy*arr_I->strides[1]);
                n_weight=lambda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*shape[1]+y;
                int qIndex=uperPointx*shape[1]+uperPointy;
                
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(n_weight>max_weight)
            {
                max_weight=n_weight;
            }
            if(LeftPointx>=0 && LeftPointx<shape[0])
            {
                float qValue=*(float *)(arr_I->data + LeftPointx*arr_I->strides[0] + LeftPointy*arr_I->strides[1]);
                n_weight=lambda*exp(-(pValue-qValue)*(pValue-qValue)/(2*sigma*sigma));
                int pIndex=x*shape[1]+y;
                int qIndex=LeftPointx*shape[1]+LeftPointy;
                g->add_edge(qIndex,pIndex,n_weight,n_weight);
            }
            if(n_weight>max_weight)
            {
                max_weight=n_weight;
            }
        }
    }
    
    max_weight =  1000 * max_weight;
    for(int x=0;x<shape[0];x++)
    {
        for(int y=0;y<shape[1];y++)
        {
            unsigned char label = * (arr_Seed->data + x*arr_Seed->strides[0] + y*arr_Seed->strides[1]);
            float s_weight = 1e-3;
            float t_weight = 1e-3;
            if(label == 127)
            {
                s_weight = max_weight;
            }
            else if(label == 255)
            {
                t_weight = max_weight;
            }
            else
            {
                float forePosibility=*(float *)(arr_fP->data + x*arr_fP->strides[0] + y*arr_fP->strides[1]);
                if(forePosibility<0.001)
                {
                    forePosibility = 0.001;
                }
                float backPosibility=*(float *)(arr_bP->data + x*arr_bP->strides[0] + y*arr_bP->strides[1]);
                if(backPosibility<0.001)
                {
                    backPosibility = 0.001;
                }
                s_weight=-log(backPosibility);
                t_weight=-log(forePosibility);
            }
            int pIndex=x*shape[1]+y;
            g->add_tweights(pIndex,s_weight,t_weight);
        }
    }
    double flow;
    flow = g->maxflow();
    //    printf("max flow: %f\n",flow);
    
    int outshape[2];
    outshape[0]=shape[0];
    outshape[1]=shape[1];
    PyArrayObject * labels = (PyArrayObject*)  PyArray_FromDims(2, outshape, NPY_INT8);
    for(int x=0; x<shape[0]; x++)
    {
        for(int y=0; y<shape[1]; y++)
        {
            int Index = x*shape[1]+y;
            *(labels->data + x*labels->strides[0] + y*labels->strides[1]) = 1-g->what_segment(Index);
        }
    }
    delete g;
    
    Py_DECREF(arr_I);
    Py_DECREF(arr_fP);
    Py_DECREF(arr_bP);
    Py_INCREF(labels);
    return PyArray_Return(labels);
}

static PyMethodDef Methods[] = {
    {"max_flow",  max_flow_wrapper, METH_VARARGS, "computing max flow"},
    {"interactive_max_flow",  interactive_max_flow_wrapper, METH_VARARGS, "computing max flow"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initmax_flow(void) {
    (void) Py_InitModule("max_flow", Methods);
    import_array();
}
