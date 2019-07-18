#include "dense_crf3d.cpp"

PyMODINIT_FUNC
initdense_crf3d(void) {
   (void) Py_InitModule("dense_crf3d", Methods);
   import_array();
}
