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
x1 = in_data[i][0]; \
x2 = in_data[i][1]; \
x3 = in_data[i][2];

double calc_residual(double (*in_data)[3], unsigned int in_data_len,
        double out_mag, double A[3][3], double b[3], gsl_matrix *residual) {

    SPLIT_VARS
    double e = 0;
    for (unsigned int i = 0; i < in_data_len; i++) {
        SPLIT_INPUTS
        double my_mag = pow(a6*x3+b3, 2)+pow(a5*x3+a4*x2+b2, 2)+pow(a3*x3+a2*x2+a1*x1+b1, 2);
        double my_residual = out_mag*out_mag - my_mag;
        gsl_matrix_set(residual, i, 0, my_residual);
        e += my_residual*my_residual;
    }

    return e;
}

unsigned int cal(double (*in_data)[3], unsigned int in_data_len,
        double out_mag, double A[3][3], double b[3],
        unsigned int max_iter) {

    gsl_permutation *p = gsl_permutation_alloc(9);
    gsl_matrix *J = gsl_matrix_alloc(in_data_len, 9);
    gsl_matrix *Jt = gsl_matrix_alloc(9, in_data_len);
    gsl_matrix *JtJ = gsl_matrix_alloc(9, 9);
    gsl_matrix *JtJinv = gsl_matrix_alloc(9, 9);
    gsl_matrix *JtJinvJt = gsl_matrix_alloc(9, in_data_len);
    gsl_matrix *residual = gsl_matrix_alloc(in_data_len, 1);
    gsl_matrix *B = gsl_matrix_alloc(9, 1);
    gsl_matrix *JtJinvJtr = gsl_matrix_alloc(9, 1);

    for (unsigned int j = 0; j < max_iter; j++) {
        calc_residual(in_data, in_data_len, out_mag, A, b, residual); 

        SPLIT_VARS
        gsl_matrix_set(B, 0, 0, a1);
        gsl_matrix_set(B, 1, 0, a2);
        gsl_matrix_set(B, 2, 0, a3);
        gsl_matrix_set(B, 3, 0, a4);
        gsl_matrix_set(B, 4, 0, a5);
        gsl_matrix_set(B, 5, 0, a6);
        gsl_matrix_set(B, 6, 0, b1);
        gsl_matrix_set(B, 7, 0, b2);
        gsl_matrix_set(B, 8, 0, b3);
        for (unsigned int i = 0; i < in_data_len; i++) {
            SPLIT_INPUTS
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
        gsl_matrix_add(B, JtJinvJtr);

        A[0][0] = gsl_matrix_get(B, 0, 0);
        A[0][1] = gsl_matrix_get(B, 1, 0);
        A[0][2] = gsl_matrix_get(B, 2, 0);
        A[1][1] = gsl_matrix_get(B, 3, 0);
        A[1][2] = gsl_matrix_get(B, 4, 0);
        A[2][2] = gsl_matrix_get(B, 5, 0);
        b[0] = gsl_matrix_get(B, 6, 0);
        b[1] = gsl_matrix_get(B, 7, 0);
        b[2] = gsl_matrix_get(B, 8, 0);

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

    return 1;
}





