//MFC #include "stdafx.h"
#include <complex>
#include "matrix.h"
#include "PI.h"
#include "Ident3D.h"

const int MAX_LANDMARKS = 32;

int main() {
    double a, b, c, d;
    double x1, x2, x3;
    std::complex<double> z1;
    std::complex<double> z2;
    std::complex<double> z3;
    int numRoots;
    R3Matrix m;
    double coeff[4];
    int numReals, numComplex;
    double realEigenvalues[3];
    int numEigenvectors;
    double eigenvalues[3];
    std::complex<double> complexEigenvalues[3];
    R3Vector realEigenvectors[3];

    int numLandmarks;
    double landmarkWeight[MAX_LANDMARKS];
    R3Vector landmarkCenter[MAX_LANDMARKS];
    R3Vector centroid;
    R3Vector secondMoment;

    double landmarkWeight1[MAX_LANDMARKS];
    R3Vector landmarkCenter1[MAX_LANDMARKS];
    R3Vector centroid1;
    double mass;
    R3Vector secondMoment1;

    R3Vector rotationAxis;
    double rotationAngleDegrees, rotationAngle;
    R3Matrix rotation;
    R3Vector shift;

    R3Matrix rotationDetected;
    R3Vector shiftDetected;

    int command = 1;

    while (command > 0) {
        printf(
            "Test of 3D identification package.\n"
            "Commands:\n"
            "   1 -- test of cubic equation;\n"
            "   2 -- test of eigenvalues & eigenvectors computing;\n"
            "   3 -- test of position identification;\n"
            "   0 -- quit.\n"
            "Input a command:\n"
        );
        scanf("%d", &command);
        if (command > 3 || command < 0)
            continue;
        if (command == 0)
            break;
        if (command == 1)
            goto LCubic;
        else if (command == 2)
            goto LEigen;
        else if (command == 3)
            goto LPosition;
        continue;

        LCubic: ;
        printf("Test of cubic equation.\n");
        numRoots = 1;
        while (numRoots != 0) {
            printf("Input a, b, c, d (a == 0 for the end):\n");
            scanf("%lf%lf%lf%lf", &a, &b, &c, &d);
            numRoots = SolveCubicEquation(
                a, b, c, d, x1, x2, x3, z1, z2, z3
            );
            printf("Number of real roots = %d\n", numRoots);
            if (numRoots == 3) {
                printf("Roots: %.3lf, %.3lf, %.3lf\n", x1, x2, x3);
            } else if (numRoots == 1) {
                printf(
                    "Real root: %.3lf,\n"
                    "complex roots: (%.3lf, %.3lf), (%.3lf, %.3lf)\n",
                    x1, z2.real(), z2.imag(), z3.real(), z3.imag()
                );
            }
        } // end while
        continue;

        LEigen: ;
        printf("Test of eigenvalues computing.\n");
        while (1) {
            printf("Input a 3*3 matrix (zero matrix for the end):\n");
            int zeroMatrix = 1;
            int j;
            for (j = 0; j < 9; ++j) {
                scanf("%lf", &(m[j/3][j%3]));
                if (m[j/3][j%3] != 0.)
                    zeroMatrix = 0;
            }
            if (zeroMatrix)
                break;
            printf("Matrix M:\n");
            int i;
            for (i = 0; i < 3; ++i) {
                printf(
                    "%.3lf %.3lf %.3lf\n",
                    m[i][0], m[i][1], m[i][2]
                );
            }
            printf("----");
            defineCharacteristicEquation(m, coeff);
            printf("Coeff. of characteristic equation:\n");
            printf(
                "%.3lf %.3lf %.3lf %.3lf\n",
                coeff[0], coeff[1], coeff[2], coeff[3]
            );
            computeEigenvalues(
                m, numReals, realEigenvalues,
                numComplex, complexEigenvalues
            );
            printf("Number of real eigenvalues: %d\nEigenvalues:", numReals);
            for (j = 0; j < numReals; ++j) {
                printf(" %.3lf", realEigenvalues[j]);
            }
            printf("\n");

            computeRealEigenvectors(
                m, numReals, realEigenvalues,
                numEigenvectors, realEigenvectors, eigenvalues
            );
            printf("Eigenvectors:\n");
            for (j = 0; j < numEigenvectors; ++j) {
                printf(
                    "(%.3lf, %.3lf, %.3lf), eigenvalue %.3lf\n",
                    realEigenvectors[j].x,
                    realEigenvectors[j].y,
                    realEigenvectors[j].z,
                    eigenvalues[j]
                );
            }

            printf("Result checking:\n");
            for (j = 0; j < numEigenvectors; ++j) {
                printf("Eigenvalue u: %.3lf\n", eigenvalues[j]);
                printf(
                    "Eigenvector v: (%.3lf, %.3lf, %.3lf)\n",
                    realEigenvectors[j].x,
                    realEigenvectors[j].y,
                    realEigenvectors[j].z
                );
                R3Vector v = m * realEigenvectors[j];
                printf(
                    "        M * v: (%.3lf, %.3lf, %.3lf)\n",
                    v.x, v.y, v.z
                );
                v = realEigenvectors[j] * eigenvalues[j];
                printf(
                    "        u * v: (%.3lf, %.3lf, %.3lf)\n",
                    v.x, v.y, v.z
                );
            }
        } // end while
        continue;

        LPosition: ;
        printf("Test of models position identification.\n");
        while (1) {
            printf("Input a number of landmarks (0 for the end):\n");
            scanf("%d", &numLandmarks);
            if (numLandmarks <= 0)
                break;
            if (numLandmarks > MAX_LANDMARKS)
                numLandmarks = MAX_LANDMARKS;
            int i;
            for (i = 0; i < numLandmarks; ++i) {
                printf("Input a 3D-center of %d-th landmark:\n", i+1);
                scanf(
                    "%lf%lf%lf",
                    &(landmarkCenter[i].x),
                    &(landmarkCenter[i].y),
                    &(landmarkCenter[i].z)
                );
                printf("Input a weight of %d-th landmark:\n", i+1);
                scanf("%lf", &(landmarkWeight[i]));
            }

            printf("Input a shift vector:\n");
            scanf(
                "%lf%lf%lf",
                &(shift.x), &(shift.y), &(shift.z)
            );

            while (1) {
                R3Vector rotationAxis;
                printf("Input a rotation axis (zero for end):\n");
                scanf(
                    "%lf%lf%lf",
                    &(rotationAxis.x),
                    &(rotationAxis.y),
                    &(rotationAxis.z)
                );
                if (rotationAxis == R3Vector(0., 0., 0.))
                    break;
                printf("Input an angle of rotation in degrees:\n");
                scanf("%lf", &rotationAngleDegrees);
                rotationAngle = rotationAngleDegrees * PI / 180.;

                rotation = R3Matrix(rotationAxis, rotationAngle);
                printf("Rotation matrix:\n");
                for (i = 0; i < 3; ++i) {
                    printf(
                        "%12.3lf %12.3lf %12.3lf\n",
                        rotation[i][0], rotation[i][1], rotation[i][2]
                    );
                }

                // Calculate a centroid
                mass = 0.;
                centroid = R3Vector(0., 0., 0.);
                for (i = 0; i < numLandmarks; ++i) {
                    mass += landmarkWeight[i];
                    centroid += landmarkCenter[i] * landmarkWeight[i];
                }
                if (mass > 0.) {
                    centroid *= 1./mass;
                }

                // Define transformed model
                for (i = 0; i < numLandmarks; ++i) {
                    landmarkWeight1[i] = landmarkWeight[i];
                    landmarkCenter1[i] =
                        rotation * (landmarkCenter[i] - centroid)
                        + centroid
                        + shift;
                }

                // Ok, now we are ready to perform the main test.
                define3DTransform(
                    numLandmarks,
                    landmarkCenter,
                    landmarkWeight,

                    numLandmarks,
                    landmarkCenter1,
                    landmarkWeight1,

                    rotationDetected,
                    shiftDetected
                );

                printf("Shift vector detected:\n");
                printf(
                    "(%.3lf, %.3lf, %.3lf)\n",
                    shiftDetected.x, shiftDetected.y, shiftDetected.z
                );

                printf("Rotation matrix detected:\n");
                for (i = 0; i < 3; ++i) {
                    printf(
                        "%12.3lf %12.3lf %12.3lf\n",
                        rotationDetected[i][0],
                        rotationDetected[i][1],
                        rotationDetected[i][2]
                    );
                }
                printf("----\n");
            } // end while (rotation input...

        } // end while

        continue;
    } // end while (command...

    return 0;
}
