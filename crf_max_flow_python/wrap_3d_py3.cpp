#include "max_flow3d.cpp"

static struct PyModuleDef sMaxFlow3d = 
{
    PyModuleDef_HEAD_INIT,
    "max_flow3d", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    Methods
};

PyMODINIT_FUNC PyInit_max_flow3d(void) {
    import_array();
    return PyModule_Create(&sMaxFlow3d);
}
