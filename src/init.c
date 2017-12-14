#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void edgeselect(void *, void *, void *, void *, void *, void *);
extern void fm123(void *, void *, void *, void *, void *, void *, void *);
extern void fmij0(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fmijk0(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fmk0(void *, void *, void *, void *, void *, void *, void *);
extern void int1(void *, void *, void *, void *, void *, void *, void *, void *);
extern void int2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void int3(void *, void *, void *, void *, void *);
extern void pointinashape(void *, void *, void *, void *, void *, void *, void *, void *);
extern void sortm(void *, void *, void *, void *);
extern void triangleNormals(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP sortbycolumn(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"edgeselect",      (DL_FUNC) &edgeselect,       6},
    {"fm123",           (DL_FUNC) &fm123,            7},
    {"fmij0",           (DL_FUNC) &fmij0,           15},
    {"fmijk0",          (DL_FUNC) &fmijk0,          13},
    {"fmk0",            (DL_FUNC) &fmk0,             7},
    {"int1",            (DL_FUNC) &int1,             8},
    {"int2",            (DL_FUNC) &int2,            14},
    {"int3",            (DL_FUNC) &int3,             5},
    {"pointinashape",   (DL_FUNC) &pointinashape,    8},
    {"sortm",           (DL_FUNC) &sortm,            4},
    {"triangleNormals", (DL_FUNC) &triangleNormals,  8},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"sortbycolumn", (DL_FUNC) &sortbycolumn, 3},
    {NULL, NULL, 0}
};

void R_init_alphashape3d(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
