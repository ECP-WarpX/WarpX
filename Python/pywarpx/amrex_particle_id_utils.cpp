#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_3kcompat.h"
#include <math.h>

static PyMethodDef GetIDMethods[] = {
    {NULL, NULL, 0, NULL}
};

static long uint64_t_to_long_id (uint64_t idcpu) {
    long r = 0;

    uint64_t sign = idcpu >> 63;  // extract leftmost sign bit
    uint64_t val  = ((idcpu >> 24) & 0x7FFFFFFFFF);  // extract next 39 id bits

    long lval = static_cast<long>(val);  // bc we take -
    r = (sign) ? lval : -lval;
    return r;
}

static long uint64_t_to_long_cpu (uint64_t idcpu) {
    return idcpu & 0x00FFFFFF;
}

/* The loop definition must precede the PyMODINIT_FUNC. */

static void get_id (char **args, const npy_intp *dimensions,
                    const npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in = args[0], *out = args[1];
    npy_intp in_step = steps[0], out_step = steps[1];

    for (i = 0; i < n; i++) {
        *((long *)out) = uint64_t_to_long_id( *((uint64_t *)in));
        in += in_step;
        out += out_step;
    }
}

/* The loop definition must precede the PyMODINIT_FUNC. */

static void get_cpu (char **args, const npy_intp *dimensions,
                    const npy_intp *steps, void *data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in = args[0], *out = args[1];
    npy_intp in_step = steps[0], out_step = steps[1];

    for (i = 0; i < n; i++) {
        *((long *)out) = uint64_t_to_long_cpu( *((uint64_t *)in));
        in += in_step;
        out += out_step;
    }
}

/* This a pointer to the above function */
PyUFuncGenericFunction id_funcs[1] = {&get_id};

/* This a pointer to the above function */
PyUFuncGenericFunction cpu_funcs[1] = {&get_cpu};

/* These are the input and return dtypes of the get_id functions.*/
static char id_types[2] = {NPY_UINT64, NPY_LONG};

/* These are the input and return dtypes of the get_cpu function*/
static char cpu_types[2] = {NPY_UINT64, NPY_INT};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "amrex_particle_id_utils",
    NULL,
    -1,
    GetIDMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_amrex_particle_id_utils (void)
{
    PyObject *m, *get_id, *get_cpu, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    get_id = PyUFunc_FromFuncAndData(id_funcs, NULL, id_types, 1, 1, 1,
                                     PyUFunc_None, "get_id",
                                     "get_id_docstring", 0);

    get_cpu = PyUFunc_FromFuncAndData(cpu_funcs, NULL, cpu_types, 1, 1, 1,
                                      PyUFunc_None, "get_cpu",
                                      "get_cpu_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "get_id", get_id);
    PyDict_SetItemString(d, "get_cpu", get_cpu);

    Py_DECREF(get_id);
    Py_DECREF(get_cpu);

    return m;
}
