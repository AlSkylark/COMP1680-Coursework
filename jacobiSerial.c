#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 2000

int main(int argc, char *argv[]) {

    if (argc < 4 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        fprintf(stderr, "Usage: %s <nx> <ny> <tolerance> <OPTIONAL -wop (without printing result)>\n", argv[0]);
        return 1;
    }

    bool withPrinting = true;

    const int nx = atoi(argv[1]);
    const int ny = atoi(argv[2]);
    const double tol = atof(argv[3]);

    double (*t)[ny + 2] = malloc((nx + 2) * sizeof *t);
    double (*tnew)[ny + 2] = malloc((nx + 2) * sizeof *tnew);
    if (t == NULL || tnew == NULL) {
        fprintf(stderr, "Could not allocate memory for a %d x %d grid\n", nx, ny);
        return 1;
    }

    printf("--- 2D Jacobi Solver ---\n");
    printf("Grid size (interior): %d x %d\n", nx, ny);
    printf("Tolerance: %g\n\n", tol);

    // Initial guess for all interior points
    for (int i = 0; i < nx + 2; i++) {
        for (int j = 0; j < ny + 2; j++) {
            t[i][j] = 30.0;
        }
    }

    // Set boundary conditions
    for (int i = 1; i <= nx; i++) {
        t[i][0]      = 30.0; //left
        t[i][ny + 1] = 28.0; //right
    }
    for (int j = 1; j <= ny; j++) {
        t[0][j]      = 15.0; //top
        t[nx + 1][j] = 5.0; //bottom
    }

    int iter = 0;
    double difmax = 1e6; //1,000,000

    //we'll need to parallelise this most likely
    while (difmax > tol && iter < MAX_ITER) {
        iter++;

        for (int i = 1; i <= nx; i++) {
            for (int j = 1; j <= ny; j++) {
                tnew[i][j] = (t[i-1][j] + t[i+1][j] + t[i][j-1] + t[i][j+1]) / 4.0;
            }
        }

        difmax = 0.0;
        for (int i = 1; i <= nx; i++) {
            for (int j = 1; j <= ny; j++) {
                double diff = fabs(tnew[i][j] - t[i][j]);
                if (diff > difmax) {
                    difmax = diff;
                }
            }
        }

        for (int i = 1; i <= nx; i++) {
            for (int j = 1; j <= ny; j++) {
                t[i][j] = tnew[i][j];
            }
        }
    }

    if (difmax > tol) {
        printf("Reached the iteration limit (%d) with a final max difference of %g\n", iter, difmax);
    } else {
        printf("Solver converged in %d iterations with a final max difference of %g\n", iter, difmax);
    }

    if (argc > 4 && strcmp(argv[4], "-wop") == 0) {
        withPrinting = false;
    }

    if(withPrinting){
        printf("\n--- Final Temperature ---\n");
        for (int i = 0; i < nx + 2; i++) {
            for (int j = 0; j < ny + 2; j++) {
                printf("%8.3f ", t[i][j]);
            }
            printf("\n");
        }
    }

    free(t);
    free(tnew);
    return 0;
}