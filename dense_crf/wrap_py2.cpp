#include "dense_crf.cpp"

PyMODINIT_FUNC
initdense_crf(void) {
    (void) Py_InitModule("dense_crf", Methods);
    import_array();
}
