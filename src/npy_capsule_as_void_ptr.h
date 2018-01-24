/*
 * extracted from:
 *   https://github.com/numpy/numpy/blob/v1.13.3/numpy/core/include/numpy/npy_3kcompat.h
 * licensed under BSD-3-Clause:
 *   https://github.com/numpy/numpy/blob/v1.13.3/LICENSE.txt
 */
#ifndef _NPY_CAPSULE_AS_VOID_PTR_H_
#define _NPY_CAPSULE_AS_VOID_PTR_H_

#include <Python.h>

/*
 * PyCObject functions adapted to PyCapsules.
 *
 * The main job here is to get rid of the improved error handling
 * of PyCapsules. It's a shame...
 */

#if PY_VERSION_HEX >= 0x03000000

static /*NPY_INLINE*/ void *
NpyCapsule_AsVoidPtr(PyObject *obj)
{
    void *ret = PyCapsule_GetPointer(obj, NULL);
    if (ret == NULL) {
        PyErr_Clear();
    }
    return ret;
}

#else

static /*NPY_INLINE*/ void *
NpyCapsule_AsVoidPtr(PyObject *ptr)
{
    return PyCObject_AsVoidPtr(ptr);
}

#endif

#endif
