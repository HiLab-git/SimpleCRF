#include "max_flow3d.cpp"


PyMODINIT_FUNC
initmax_flow(void) {
    (void) Py_InitModule("max_flow3d", Methods);
    import_array();
}
