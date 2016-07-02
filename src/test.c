#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "sensor_cal.h"

#define TEST_DATA_FILE "../test/test.bin"
#define TEST_LEN 1000
#define NUM_ITER 100
#define MAG 3

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

    double A[3][3] = { {5, 0, 0}, {0, 6, 0}, {0, 0, 7} };
    double b[3] = {3,3,3};
    
    printf("before cal:\n");
    printf("A = [%f %f %f;\n %f %f %f;\n %f %f %f]\n",
    A[0][0], A[0][1], A[0][2], A[1][0], A[1][1],
    A[1][2], A[2][0], A[2][1], A[2][2]);



    cal(test_data, TEST_LEN, MAG, A, b, NUM_ITER);


    printf("after cal:\n");
    printf("A = [%f %f %f;\n %f %f %f;\n %f %f %f]\n",
    A[0][0], A[0][1], A[0][2], A[1][0], A[1][1],
    A[1][2], A[2][0], A[2][1], A[2][2]);

    printf("b = [%f; %f; %f]\n", b[0], b[1], b[2]);

}

