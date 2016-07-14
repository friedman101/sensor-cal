#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "sensor_cal.h"

#define TEST_DATA_FILE "../test/test.bin"
#define TEST_LEN 100000
#define NUM_ITER 100
#define MAG 3
#define TOL 1e-4


#ifdef __MACH__
#include <sys/time.h>
#define CLOCK_MONOTONIC 0
#define UNUSED(expr) do { (void)(expr); } while (0)
//clock_gettime is not implemented on OSX
int clock_gettime(int trsh, struct timespec* t) {
    UNUSED(trsh);
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

double timespec2double(struct timespec in) {
    double out;
    out = (double) in.tv_sec;
    out += ( (double) in.tv_nsec)*1.0e-9;

    return out;
}

int main() {
    double test_data[TEST_LEN][3];
    FILE *f = fopen(TEST_DATA_FILE, "r");
    assert(f != NULL);
    double buffer[3];
    unsigned int count = 0;
    while (fread(buffer, sizeof(buffer), 1, f)) {
        memcpy(test_data[count], buffer, sizeof(buffer));
        count++;
    }
    fclose(f);

    double A[3][3] = { {1.5, 0, 0}, {2.5, 3.5, 0}, {4.5, 5.5, 6.5} };
    double b[3] = {7.5, 8.5, 9.5};

    printf("before cal:\n");
    printf("A = [%f %f %f;\n %f %f %f;\n %f %f %f]\n",
            A[0][0], A[0][1], A[0][2], A[1][0], A[1][1],
            A[1][2], A[2][0], A[2][1], A[2][2]);
    printf("b = [%f; %f; %f]\n", b[0], b[1], b[2]);

    struct timespec tstart, tend;
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    unsigned int c = cal(test_data, TEST_LEN, MAG, A, b, NUM_ITER, TOL);
    clock_gettime(CLOCK_MONOTONIC, &tend);

    double t1 = timespec2double(tstart);
    double t2 = timespec2double(tend);

    printf("took %0.2f ms to optimize %i samples\n", (t2-t1)*1e3,
            TEST_LEN);
    printf("done in %i iterations\n", c);
    printf("after cal:\n");
    printf("A = [%f %f %f;\n %f %f %f;\n %f %f %f]\n",
            A[0][0], A[0][1], A[0][2], A[1][0], A[1][1],
            A[1][2], A[2][0], A[2][1], A[2][2]);
    printf("b = [%f; %f; %f]\n", b[0], b[1], b[2]);

    return 0;

}

