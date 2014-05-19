#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>
#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

/*     C FUNCTIONS IN PACKAGE alphashape3d  */
/*     Thomas Lafarge <thowas.lafarge@gmail.com> */
/*     Beatriz Pateiro-L�pez <beatriz.pateiro@usc.es> */
/* ========================================================== */

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int intcmp(const void *aa, const void *bb)
{
    const int *a = aa, *b = bb;
    return (*a < *b) ? -1 : (*a > *b);
}


void sort(int n, int *a) {
    qsort(a, n, sizeof(int), intcmp);
}

/*     SUBROUTINE SORTM (ORDM,NR,NC,M,A) */
/*     SORTS THE ROWS OF THE INTEGER MATRIX M IN INCREASING ORDER */
void sortm(int *ordm, int *nr, int *nc, int * m) {

	int m_dim1, i1, i2;

	int i, j, k;

	m_dim1 = *nr;

	/* Function Body */
	int *a = malloc(*nc * sizeof(int));
	i1 = *nr;
	for (i = 0; i < i1; ++i) {
		i2 = *nc;
		for (k = 0; k < i2; ++k) {
			a[k] = m[i + k * m_dim1];
		}
		sort(*nc, a);
		for (j = 0; j < i2; ++j) {
			ordm[i + j * m_dim1] = a[j];
		}
	}
	free(a);
}

/*     SUBROUTINE Fm123(x,nrx,t1,t2,t3,nt,m123) */
/*     COMPUTES DETERMINANTS Mijk */
/*     Maybe use ti -1;       */
void fm123(double *x, int *nrx, int *t1, int *t2, int *t3, int *nt, double *m123) {
	/* System generated locals */
	int x_dim1, i1;

	/* Local variables */
	int i, j;
	double a1, a2, a3, aux[9];


	x_dim1 = *nrx;

	/* Function Body */
	i1 = *nt;
	for (i = 0; i < i1; ++i) {
		for (j = 0; j < 3; ++j) {
			aux[0 + j * 3] = x[t1[i]-1 + j * x_dim1];
			aux[1 + j * 3] = x[t2[i]-1 + j * x_dim1];
			aux[2 + j * 3] = x[t3[i]-1 + j * x_dim1];
		}
		a1 = aux[0] * (aux[4] * aux[8] - aux[7] * aux[5]);
		a2 = -aux[3] * (aux[1] * aux[8] - aux[7] * aux[2]);
		a3 = aux[6] * (aux[1] * aux[5] - aux[4] * aux[2]);
		m123[i] = a1 + a2 + a3;
	}
}

/*     SUBROUTINE Fmk0(x,nrx,e1,e2,ned,mk0,num) */
/*     COMPUTES DETERMINANTS Mk0, k=1,2,3 AND NUM.RHO2 */
/*Maybe use ei -1*/
void fmk0(double *x, int *nrx, int *e1, int *e2, int *ned, double *mk0,
		double *num) {

	int x_dim1, mk0_dim1, i1;
	double d1, d2, d3;

	int i;

	/* Parameter adjustments */
	x_dim1 = *nrx;
	mk0_dim1 = *ned;

	/* Function Body */
	i1 = *ned;
	for (i = 0; i < i1; ++i) {
		mk0[i + mk0_dim1 * 0] = x[e1[i]-1 + x_dim1 * 0] - x[e2[i]-1 + x_dim1 * 0];
		mk0[i + mk0_dim1 * 1] = x[e1[i]-1 + x_dim1 * 1] - x[e2[i]-1 + x_dim1 * 1];
		mk0[i + mk0_dim1 * 2] = x[e1[i]-1 + x_dim1 * 2] - x[e2[i]-1 + x_dim1 * 2];

		/* Computing num(i) */
		d1 = mk0[i + mk0_dim1 * 0];
		d2 = mk0[i + mk0_dim1 * 1];
		d3 = mk0[i + mk0_dim1 * 2];
		num[i] = d1 * d1 + d2 * d2 + d3 * d3;
	}
}

/*     SUBROUTINE Fmij0(x,n,t1,t2,t3,nt,ta,ie,ned,m0,m23,m13,m12,nu2,nu3) */
/*     COMPUTES DETERMINANTS Mij0, AND NUM.RHO3 */
void fmij0(double *x, int *n, int *t1, int *t2, int *t3, int *nt, int *ta,
		int *ie, int *ned, double *m0, double *m23, double *m13, double *m12,
		double *nu2, double *nu3) {
	/* System generated locals */
	int x_dim1, ta_dim1, m0_dim1, i1;

	/* Local variables */
	int i;
	double au1, au2, au3;

	/* Parameter adjustments */
	x_dim1 = *n;
	ta_dim1 = *nt;
	m0_dim1 = *ned;

	/* Function Body */
	i1 = *nt;
	for (i = 0; i < i1; ++i) {
		au1 = x[t1[i]-1 + x_dim1 * 1] * m0[ie[ta[i + ta_dim1 * 2]-1]-1 + m0_dim1 * 2];
		au2 = -x[t2[i]-1 + x_dim1 * 1]
		         * m0[ie[ta[i + ta_dim1 * 1]-1]-1 + m0_dim1 * 2];
		au3 = x[t3[i]-1 + x_dim1 * 1] * m0[ie[ta[i + ta_dim1 * 0]-1]-1 + m0_dim1 * 2];
		m23[i] = au1 + au2 + au3;
		au1 = x[t1[i]-1 + x_dim1 * 0] * m0[ie[ta[i + ta_dim1 * 2]-1]-1 + m0_dim1 * 2];
		au2 = -x[t2[i]-1 + x_dim1 * 0]
		         * m0[ie[ta[i + ta_dim1 * 1]-1]-1 + m0_dim1 * 2];
		au3 = x[t3[i]-1 + x_dim1 * 0] * m0[ie[ta[i + ta_dim1 * 0]-1]-1 + m0_dim1 * 2];
		m13[i] = au1 + au2 + au3;
		au1 = x[t1[i]-1 + x_dim1 * 0] * m0[ie[ta[i + ta_dim1 * 2]-1]-1 + m0_dim1 * 1];
		au2 = -x[t2[i]-1 + x_dim1 * 0]
		         * m0[ie[ta[i + ta_dim1 * 1]-1]-1 + m0_dim1 * 1];
		au3 = x[t3[i]-1 + x_dim1 * 0] * m0[ie[ta[i + ta_dim1 * 0]-1]-1 + m0_dim1 * 1];
		m12[i] = au1 + au2 + au3;
		nu3[i] = nu2[ie[ta[i + ta_dim1 * 0]-1]-1] * nu2[ie[ta[i + ta_dim1 * 1]-1]-1]
		                                            * nu2[ie[ta[i + ta_dim1 * 2]-1]-1];
	}
}

/*     SUBROUTINE Fmijk0(x,n,tc,ntc,tca,it,nt,m23,m13,m12,m234,m134,m124) */
/*     COMPUTES DETERMINANTS Mijk0 */
void fmijk0(double *x, int *n, int *tc, int *ntc, int *tca, int *it, int *nt,
		double *m23, double *m13, double *m12, double *m234, double *m134,
		double *m124) {
	/* System generated locals */
	int tc_dim1, tca_dim1, i1;

	/* Local variables */
	int i;
	double au1, au2, au3, au4;

	/* Parameter adjustments */
	tca_dim1 = *ntc;
	tc_dim1 = *ntc;

	/* Function Body */
	i1 = *ntc;
	for (i = 0; i < i1; ++i) {

		au1 =  x[tc[i + tc_dim1 * 0]-1] * m23[it[tca[i + tca_dim1 * 3]-1]-1];
		au2 = -x[tc[i + tc_dim1 * 1]-1] * m23[it[tca[i + tca_dim1 * 2]-1]-1];
		au3 =  x[tc[i + tc_dim1 * 2]-1] * m23[it[tca[i + tca_dim1 * 1]-1]-1];
		au4 = -x[tc[i + tc_dim1 * 3]-1] * m23[it[tca[i + tca_dim1 * 0]-1]-1];
		m234[i] = au1 + au2 + au3 + au4;

		au1 =  x[tc[i + tc_dim1 * 0]-1] * m13[it[tca[i + tca_dim1 * 3]-1]-1];
		au2 = -x[tc[i + tc_dim1 * 1]-1] * m13[it[tca[i + tca_dim1 * 2]-1]-1];
		au3 =  x[tc[i + tc_dim1 * 2]-1] * m13[it[tca[i + tca_dim1 * 1]-1]-1];
		au4 = -x[tc[i + tc_dim1 * 3]-1] * m13[it[tca[i + tca_dim1 * 0]-1]-1];
		m134[i] = au1 + au2 + au3 + au4;

		au1 =  x[tc[i + tc_dim1 * 0]-1] * m12[it[tca[i + tca_dim1 * 3]-1]-1];
		au2 = -x[tc[i + tc_dim1 * 1]-1] * m12[it[tca[i + tca_dim1 * 2]-1]-1];
		au3 =  x[tc[i + tc_dim1 * 2]-1] * m12[it[tca[i + tca_dim1 * 1]-1]-1];
		au4 = -x[tc[i + tc_dim1 * 3]-1] * m12[it[tca[i + tca_dim1 * 0]-1]-1];
		m124[i] = au1 + au2 + au3 + au4;

	}
}

/*     SUBROUTINE int3(ntri2,rf1,rf2,l3,u3) */
/*     COMPUTES INTERVALS OF ALPHA VALUES FOR TRIANGLES */
void int3(int *ntri2, double *rf1, double *rf2, double *l3, double *u3) {
	/* System generated locals */
	int i1;

	/* Local variables */
	int i;

	/* Function Body */
	i1 = *ntri2;
	for (i = 0; i < i1; ++i) {
		if (rf1[i] > rf2[i]) {
			u3[i] = rf1[i];
			l3[i] = rf2[i];
		} else {
			u3[i] = rf2[i];
			l3[i] = rf1[i];
		}
	}
}

/*     SUBROUTINE int2(dup,ned,ntri,itro,l3,u3,l2,u2,tra,aux,ie,nrh,m,ia) */
/*     COMPUTES INTERVALS OF ALPHA VALUES FOR EDGES */
void int2(int *dup, int *ned, int *ntri, int *itro, double *l3, double *u3,
		double *l2, double *u2, int *tra, int *aux, int *ie, double *nrh,
		double *m, int *ia) {
	/* System generated locals */
	int m_dim1, tra_dim1, i1, i2;
	double d1, d2;

	/* Local variables */
	int h, i, j, k, ed[3];
	double at = 0.0;
	double aa1, aa2, aa3;
	int edo[3];

	/* Parameter adjustments */
	m_dim1 = *ned;
	tra_dim1 = *ntri;

	/* Function Body */
	k = 0;
	i1 = *ned;
	for (i = 0; i < i1; ++i) {
		ia[i] = 0;
		l2[i] = 1e16;
		u2[i] = -1e16;
		i2 = dup[i];
		for (j = 0; j < i2; ++j) {
			/* Computing MIN */
			d1 = l2[i];
			d2 = l3[itro[k]-1];
			l2[i] = MIN(d1, d2);
			/* Computing MAX */
			d1 = u2[i];
			d2 = u3[itro[k]-1];
			u2[i] = MAX(d1, d2);
			if (ia[i] == 0) {
				for (h = 0; h < 3; ++h) {
					ed[h] = tra[itro[k]-1 + h * tra_dim1];
					edo[h] = ie[ed[h]-1];
				}
				if (aux[k] == 1) {
					/* Computing 2nd power */
					d1 = m[edo[1]-1 + m_dim1 * 0] + m[edo[2]-1 + m_dim1 * 0];
					aa1 = d1 * d1;
					/* Computing 2nd power */
					d1 = m[edo[1]-1 + m_dim1 * 1] + m[edo[2]-1 + m_dim1 * 1];
					aa2 = d1 * d1;
					/* Computing 2nd power */
					d1 = m[edo[1]-1 + m_dim1 * 2] + m[edo[2]-1 + m_dim1 * 2];
					aa3 = d1 * d1;
					at = nrh[i] - aa1 - aa2 - aa3;
				} else if (aux[k] == 2) {
					/* Computing 2nd power */
					d1 = m[edo[0]-1 + m_dim1 * 0] - m[edo[2]-1 + m_dim1 * 0];
					aa1 = d1 * d1;
					/* Computing 2nd power */
					d1 = m[edo[0]-1 + m_dim1 * 1] - m[edo[2]-1 + m_dim1 * 1];
					aa2 = d1 * d1;
					/* Computing 2nd power */
					d1 = m[edo[0]-1 + m_dim1 * 2] - m[edo[2]-1 + m_dim1 * 2];
					aa3 = d1 * d1;
					at = nrh[i] - aa1 - aa2 - aa3;
				} else if (aux[k] == 3) {
					/* Computing 2nd power */
					d1 = -m[edo[0]-1 + m_dim1 * 0] - m[edo[1]-1 + m_dim1 * 0];
					aa1 = d1 * d1;
					/* Computing 2nd power */
					d1 = -m[edo[0]-1 + m_dim1 * 1] - m[edo[1]-1 + m_dim1 * 1];
					aa2 = d1 * d1;
					/* Computing 2nd power */
					d1 = -m[edo[0]-1 + m_dim1 * 2] - m[edo[1]-1 + m_dim1 * 2];
					aa3 = d1 * d1;
					at = nrh[i] - aa1 - aa2 - aa3;
				}
				if (at > 0.f) {
					ia[i] = 1;
				}
			}
			++k;
		}
	}
}

/*     SUBROUTINE int1(dup,nvt,ned,iedo,l2,u2,l1,u1) */
/*     COMPUTES INTERVALS OF ALPHA VALUES FOR VERTICES */
void int1(int *dup, int *nvt, int *ned, int *iedo, double *l2, double *u2,
		double *l1, double *u1) {
	/* System generated locals */
	int i1, i2;
	double d1, d2;

	/* Local variables */
	int i, j, k;

	/* Function Body */
	k = 0;
	i1 = *nvt;
	for (i = 0; i < i1; ++i) {
		l1[i] = 1e16;
		u1[i] = -1e16;
		i2 = dup[i];
		for (j = 0; j < i2; ++j) {
			/* Computing MIN */
			d1 = l1[i];
			d2 = l2[iedo[k]-1];
			l1[i] = MIN(d1, d2);
			/* Computing MAX */
			d1 = u1[i];
			d2 = u2[iedo[k]-1];
			u1[i] = MAX(d1, d2);
			++k;
		}
	}
}

/*     SUBROUTINE edgeSelect(n, x, ed, nb, alpha, nbfinal) */
/*     SELECTION OF EDGES */
void edgeselect(int *n, double *x, int *ed, int *nb, double *alpha, int *nbfinal) {
	/* System generated locals */
	double d1;
	int x_dim1, ed_dim1;

	/* Local variables */
	int i;
	double dist;

	/* Function Body */
	*nbfinal = 0;
	ed_dim1 = *nb;
	x_dim1 = *n;
	for (i = 0; i < ed_dim1; ++i) {
		/* Computing 2nd power */
		d1 = x[ed[i + ed_dim1 * 0]-1 + x_dim1 * 0]
		       - x[ed[i + ed_dim1 * 1]-1 + x_dim1 * 0];
		dist = d1 * d1;
		/* Computing 2nd power */
		d1 = x[ed[i + ed_dim1 * 0]-1 + x_dim1 * 1]
		       - x[ed[i + ed_dim1 * 1]-1 + x_dim1 * 1];
		dist += d1 * d1;
		/* Computing 2nd power */
		d1 = x[ed[i + ed_dim1 * 0]-1 + x_dim1 * 2]
		       - x[ed[i + ed_dim1 * 1]-1 + x_dim1 * 2];
		dist = sqrt(dist + d1 * d1);
		if (dist < *alpha * 2) {
			ed[*nbfinal] = i+1;
			++(*nbfinal);
		}
	}
}



struct myData {
    int data[3];
    int orig_pos; // will hold the original position of data on your array
};

int myData_compare (const void* a, const void* b) {
	int value = 1;
	if (((struct myData*)a)->data[0] < ((struct myData*)b)->data[0])
		value= -1;
	else if (((struct myData*)a)->data[0] == ((struct myData*)b)->data[0] && ((struct myData*)a)->data[1] < ((struct myData*)b)->data[1])
		value= -1;
	else if (((struct myData*)a)->data[0] == ((struct myData*)b)->data[0] && ((struct myData*)a)->data[1] == ((struct myData*)b)->data[1] &&  ((struct myData*)a)->data[2] < ((struct myData*)b)->data[2])
		value= -1;

	if (((struct myData*)a)->data[0] == ((struct myData*)b)->data[0] && ((struct myData*)a)->data[1] == ((struct myData*)b)->data[1] &&  ((struct myData*)a)->data[2] == ((struct myData*)b)->data[2])
		value= 0;

	return value;
}







SEXP sortbycolumn(SEXP t1,SEXP t2,SEXP t3) {
	int i, n;
	int *v1,* v2,* v3;
	v1 = INTEGER(t1);
	v2 = INTEGER(t2);
	v3 = INTEGER(t3);
	n=length(t1);
    struct myData* array = (struct myData*) malloc(n * sizeof(struct myData));
	for(i=0;i<n;i++)
	{
		array[i].data[0]=v1[i] ;
		array[i].data[1]=v2[i] ;
		array[i].data[2]=v3[i] ;
		array[i].orig_pos=i;
	}
    qsort(array, n, sizeof(struct myData), &myData_compare);

    SEXP index;
    PROTECT (index = allocVector(INTSXP,n));
    for (size_t i = 0; i < n; i++) {
        INTEGER(index)[i] = array[i].orig_pos +1;
    }
    UNPROTECT(1);

    return(index);
}



