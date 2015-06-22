#include <Python.h>

static PyObject *_check_parallel(PyObject *self, PyObject *args);

static PyObject *_check_parallel(PyObject *self, PyObject *args)
{
	int can_parallelize;
	#if defined (_OPENMP)
	can_parallelize = 1;
	#else
	can_parallelize = 0;	
	#endif

	return Py_BuildValue("i", can_parallelize); 
}

static char _check_parallel_doc[] = """Checks whether parallelization is possible.""";

static PyMethodDef _check_parallel_methods[] = {
  {"_check_parallel", _check_parallel, METH_VARARGS, _check_parallel_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _check_parallel_module = {
		PyModuleDef_HEAD_INIT,
		"_check_parallel",
		_check_parallel_doc,
		-1, 
		_check_parallel_methods
	};

	PyMODINIT_FUNC
	PyInit__check_parallel(void)
	{
		return PyModule_Create(&_check_parallel_module);
	}
#else

	void init_check_parallel(void)
	{
	  Py_InitModule("_check_parallel", _check_parallel_methods);
	}
#endif

