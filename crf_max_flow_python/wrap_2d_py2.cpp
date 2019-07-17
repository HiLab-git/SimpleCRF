#include "max_flow.cpp"


PyMODINIT_FUNC
initmax_flow(void) {
    (void) Py_InitModule("max_flow", Methods);
    import_array();
}
