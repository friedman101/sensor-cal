#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define SPLIT_VARS \
double a1, a2, a3, a4, a5, a6; \
double b1, b2, b3; \
a1 = A[0][0]; a2 = A[0][1]; a3 = A[0][2]; \
a4 = A[1][1]; a5 = A[1][2]; a6 = A[2][2]; \
b1 = b[0]; b2 = b[1]; b3 = b[2];

#define SPLIT_INPUTS \
double x1, x2, x3; \
x1 = x[i][0]; \
x2 = x[i][1]; \
x3 = x[i][2];

#define PACK_B \
gsl_matrix_set(B, 0, 0, a1); \
gsl_matrix_set(B, 1, 0, a2); \
gsl_matrix_set(B, 2, 0, a3); \
gsl_matrix_set(B, 3, 0, a4); \
gsl_matrix_set(B, 4, 0, a5); \
gsl_matrix_set(B, 5, 0, a6); \
gsl_matrix_set(B, 6, 0, b1); \
gsl_matrix_set(B, 7, 0, b2); \
gsl_matrix_set(B, 8, 0, b3); \

#define UNPACK_B \
A[0][0] = gsl_matrix_get(B, 0, 0); \
A[0][1] = gsl_matrix_get(B, 1, 0); \
A[0][2] = gsl_matrix_get(B, 2, 0); \
A[1][1] = gsl_matrix_get(B, 3, 0); \
A[1][2] = gsl_matrix_get(B, 4, 0); \
A[2][2] = gsl_matrix_get(B, 5, 0); \
b[0] = gsl_matrix_get(B, 6, 0); \
b[1] = gsl_matrix_get(B, 7, 0); \
b[2] = gsl_matrix_get(B, 8, 0);


double calc_residual(double (*x)[3], unsigned int x_len,
        double x_mag, double A[3][3], double b[3], gsl_matrix *residual) {

    SPLIT_VARS;
    double e = 0;
    for (unsigned int i = 0; i < x_len; i++) {
        SPLIT_INPUTS;
        double p1 = a6*x3+b3;
        double p2 = a5*x3+a4*x2+b2;
        double p3 = a3*x3+a2*x2+a1*x1+b1;
        double my_mag = p1*p1 + p2*p2 + p3*p3;
        double my_residual = x_mag*x_mag - my_mag;
        gsl_matrix_set(residual, i, 0, my_residual);
        e += fabs(my_residual);
    }

    return e;
}

unsigned int cal(double (*x)[3], unsigned int x_len,
        double x_mag, double A[3][3], double b[3],
        unsigned int max_iter, double tol) {

    A[1][0] = 0;
    A[2][0] = 0;
    A[2][1] = 0;

    gsl_permutation *p = gsl_permutation_alloc(9);
    gsl_matrix *J = gsl_matrix_alloc(x_len, 9);
    gsl_matrix *Jt = gsl_matrix_alloc(9, x_len);
    gsl_matrix *JtJ = gsl_matrix_alloc(9, 9);
    gsl_matrix *JtJinv = gsl_matrix_alloc(9, 9);
    gsl_matrix *JtJinvJt = gsl_matrix_alloc(9, x_len);
    gsl_matrix *residual = gsl_matrix_alloc(x_len, 1);
    gsl_matrix *B = gsl_matrix_alloc(9, 1);
    gsl_matrix *JtJinvJtr = gsl_matrix_alloc(9, 1);

    unsigned int j;
    double e_old = 0;
    for (j = 0; j < max_iter; j++) {
        double e = calc_residual(x, x_len, x_mag, A, b, residual);
        double pct_improve = (e - e_old)/e_old;
        e_old = e;
        if (pct_improve < tol && pct_improve >= 0)
            break;

        SPLIT_VARS;
        for (unsigned int i = 0; i < x_len; i++) {
            SPLIT_INPUTS;
            double J_mat[9] = {2*x1*(a3*x3+a2*x2+a1*x1+b1),
                2*x2*(a3*x3+a2*x2+a1*x1+b1),
                2*x3*(a3*x3+a2*x2+a1*x1+b1),
                2*x2*(a5*x3+a4*x2+b2),
                2*x3*(a5*x3+a4*x2+b2),
                2*x3*(a6*x3+b3),
                2*(a3*x3+a2*x2+a1*x1+b1),
                2*(a5*x3+a4*x2+b2),
                2*(a6*x3+b3)};
            for (unsigned int k = 0; k < 9; k++)
                gsl_matrix_set(J, i, k, J_mat[k]);
        }

        gsl_matrix_transpose_memcpy(Jt, J);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Jt, J, 0.0, JtJ); 
        int s;
        gsl_linalg_LU_decomp(JtJ, p, &s);
        gsl_linalg_LU_invert(JtJ, p, JtJinv);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, JtJinv, Jt, 0.0, JtJinvJt); 
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, JtJinvJt, residual, 0.0, JtJinvJtr); 
        PACK_B;
        gsl_matrix_add(B, JtJinvJtr);
        UNPACK_B;

    }

    gsl_matrix_free(JtJinvJtr);
    gsl_matrix_free(B);
    gsl_matrix_free(residual);
    gsl_matrix_free(JtJinvJt);
    gsl_matrix_free(JtJinv);
    gsl_matrix_free(JtJ);
    gsl_matrix_free(Jt);
    gsl_matrix_free(J);
    gsl_permutation_free(p);

    return j;
}





