//MFC #include "stdafx.h"
#include <math.h>
#include "matrix.h"
#include "PI.h"
#include "Ident3D.h"

void define3DTransform(
    int numLandmarks0,
    const R3Vector* landmarkPos0,
    const double* landmarkWeight0,
    int numLandmarks1,
    const R3Vector* landmarkPos1,
    const double* landmarkWeight1,
    R3Matrix& rotation,
    R3Vector& shift
) {
    R3Vector center0, center1;
    double mass0, mass1;
    R3Matrix inertionMatrix0, inertionMatrix1;
    R3Vector secondMoment0, secondMoment1;
    double eigenvalues0[3], eigenvalues1[3];
    R3Vector axes0[3], axes1[3];

    define3DCentroid(
        numLandmarks0, landmarkPos0, landmarkWeight0,
        center0, mass0
    );

#   ifdef TEST_PRINT
    int i;
    printf("Centroid of first model:\n");
    printf("(%.3lf, %.3lf, %.3lf)\n", center0.x, center0.y, center0.z);
#   endif

    define3DCentroid(
        numLandmarks1, landmarkPos1, landmarkWeight1,
        center1, mass1
    );

#   ifdef TEST_PRINT
    printf("Centroid of second model:\n");
    printf("(%.3lf, %.3lf, %.3lf)\n", center1.x, center1.y, center1.z);
#   endif

    shift = center1 - center0;

#   ifdef TEST_PRINT
    printf("Shift vector calculated:\n");
    printf("(%.3lf, %.3lf, %.3lf)\n", shift.x, shift.y, shift.z);
#   endif

    define3DInertionMatrix(
        numLandmarks0, landmarkPos0, landmarkWeight0,
        center0,
        inertionMatrix0,
        secondMoment0
    );

#   ifdef TEST_PRINT
    printf("Inertia matrix of first model:\n");
    for (i = 0; i < 3; ++i) {
        printf(
            "%.3lf %.3lf %.3lf\n",
            inertionMatrix0[i][0],
            inertionMatrix0[i][1],
            inertionMatrix0[i][2]
        );
    }
    printf("Second moment vector of first model:\n");
    printf(
        "(%.3lf, %.3lf, %.3lf)\n",
        secondMoment0.x, secondMoment0.y, secondMoment0.z
    );
#   endif

    define3DInertionMatrix(
        numLandmarks1, landmarkPos1, landmarkWeight1,
        center1,
        inertionMatrix1,
        secondMoment1
    );

#   ifdef TEST_PRINT
    printf("Inertia matrix of second model:\n");
    for (i = 0; i < 3; ++i) {
        printf(
            "%.3lf %.3lf %.3lf\n",
            inertionMatrix1[i][0],
            inertionMatrix1[i][1],
            inertionMatrix1[i][2]
        );
    }
    printf("Second moment vector of second model:\n");
    printf(
        "(%.3lf, %.3lf, %.3lf)\n",
        secondMoment1.x, secondMoment1.y, secondMoment1.z
    );
#   endif

    define3DInertionAxes(
        inertionMatrix0,
        eigenvalues0,
        axes0
    );

#   ifdef TEST_PRINT
    printf("Eigenvalues of first model:\n");
    printf(
        "%.3lf, %.3lf, %.3lf\n",
        eigenvalues0[0], eigenvalues0[1], eigenvalues0[2]
    );
    printf("Inertia axes of first model:\n");
    for (i = 0; i < 3; ++i) {
        printf(
            "(%.3lf, %.3lf, %.3lf)\n",
            axes0[i].x, axes0[i].y, axes0[i].z
        );
    }
    printf("----");
#   endif

    define3DInertionAxes(
        inertionMatrix1,
        eigenvalues1,
        axes1
    );

#   ifdef TEST_PRINT
    printf("Eigenvalues of second model:\n");
    printf(
        "%.3lf, %.3lf, %.3lf\n",
        eigenvalues1[0], eigenvalues1[1], eigenvalues1[2]
    );
    printf("Inertia axes of second model:\n");
    for (i = 0; i < 3; ++i) {
        printf(
            "(%.3lf, %.3lf, %.3lf)\n",
            axes1[i].x, axes1[i].y, axes1[i].z
        );
    }
    printf("----");
#   endif

    define3DRotation(
        eigenvalues0, axes0, secondMoment0,
        eigenvalues1, axes1, secondMoment1,
        rotation
    );
}

void define3DCentroid(
    int numLandmarks,
    const R3Vector* landmarkPos,
    const double* landmarkWeight,
    R3Vector& center,
    double& mass
) {
    int i;
    center = R3Vector(0., 0., 0.);
    mass = 0.;
    for (i = 0; i < numLandmarks; ++i) {
        center += landmarkPos[i] * landmarkWeight[i];
        mass += landmarkWeight[i];
    }
    if (mass > 0.) {
        center *= 1. / mass;
    }
}

void define3DInertionMatrix(
    int numLandmarks,
    const R3Vector* landmarkPos,
    const double* landmarkWeight,
    const R3Vector& center,
    R3Matrix& inertionMatrix,
    R3Vector& secondMoment
) {
    int i;
    double Jxx = 0., Jxy = 0., Jxz = 0.,
        Jyy = 0., Jyz = 0.,
        Jzz = 0.;
    double mass2 = 0;
    secondMoment = R3Vector(0., 0., 0.);
    for (i = 0; i < numLandmarks; ++i) {
        R3Vector dv = landmarkPos[i] - center;
        double w = landmarkWeight[i];
        Jxx += (dv.y*dv.y + dv.z*dv.z)*w;
        Jyy += (dv.x*dv.x + dv.z*dv.z)*w;
        Jzz += (dv.x*dv.x + dv.y*dv.y)*w;
        Jxy += dv.x*dv.y*w;
        Jxz += dv.x*dv.z*w;
        Jyz += dv.y*dv.z*w;
        secondMoment += dv * (dv.Length() * (w*w));
        mass2 += (w*w);
    }
    if (mass2 > 0.) {
        secondMoment *= 1. / mass2;
    }
    /*???
    if (numLandmarks > 0) {
        Jxx /= (double) numLandmarks;
        Jyy /= (double) numLandmarks;
        Jzz /= (double) numLandmarks;
        Jxy /= (double) numLandmarks;
        Jxz /= (double) numLandmarks;
        Jyz /= (double) numLandmarks;
    }
    ???*/
    inertionMatrix[0][0] = Jxx;
    inertionMatrix[1][1] = Jyy;
    inertionMatrix[2][2] = Jzz;
    inertionMatrix[0][1] = (-Jxy);
    inertionMatrix[1][0] = (-Jxy);
    inertionMatrix[0][2] = (-Jxz);
    inertionMatrix[2][0] = (-Jxz);
    inertionMatrix[1][2] = (-Jyz);
    inertionMatrix[2][1] = (-Jyz);
}

void define3DInertionAxes(
    const R3Matrix& inertionMatrix,
    double eigenvalues[3],
    R3Vector axes[3]
) {
    int numReals, numComplex;
    double realEigenvalues[3];
    std::complex<double> complexEigenvalues[3];
    int numEigenvectors;
    R3Vector realEigenvectors[3];
    double eigenvals[3];

    computeEigenvalues(
        inertionMatrix, numReals, realEigenvalues,
        numComplex, complexEigenvalues
    );

    computeRealEigenvectors(
        inertionMatrix, numReals, realEigenvalues,
        numEigenvectors, realEigenvectors, eigenvals
    );

    //... ASSERT(numEigenvectors == 3);
    int i;
    for (i = 0; i < numEigenvectors; ++i) {
        eigenvalues[i] = eigenvals[i];
        axes[i] = realEigenvectors[i];
        axes[i].Normalize();
    }
    // Make a positive orientation
    if (axes[2] * axes[0].VectorProduct(axes[1]) < 0.) {
        axes[2] *= (-1.);   // Invert a direction of axis
    }
}

void define3DRotation(
    const double eigenvalues0[3],
    const R3Vector axes0[3],
    const R3Vector& secondMoment0,
    const double eigenvalues1[3],
    const R3Vector axes1[3],
    const R3Vector& secondMoment1,
    R3Matrix& rotation
) {
    R3Vector moment0(secondMoment0);
    moment0.Normalize();
    R3Vector moment1(secondMoment1);
    moment1.Normalize();

    computeTransitionMatrix(
        axes0, axes1, rotation
    );

    // We should consider 4 cases
    double cosines[4], maxCosine;
    int maxCosineIdx;
    //... R3Matrix rotationX(
    //...     1., 0.,  0.,
    //...     0., -1., 0.,
    //...     0., 0., -1.
    //... );
    R3Matrix rotationX(axes1[0], PI);

    //... R3Matrix rotationZ(
    //...     -1., 0.,  0.,
    //...     0.,  -1., 0.,
    //...     0.,  0.,  1.
    //... );
    R3Matrix rotationZ(axes1[2], PI);

    R3Vector moment = rotation * moment0;
    cosines[0] = moment1 * moment;
    cosines[1] = moment1 * (rotationX * moment);
    cosines[2] = moment1 * (rotationZ * moment);
    cosines[3] = moment1 * (rotationZ * (rotationX * moment));

    maxCosine = cosines[0]; maxCosineIdx = 0;

    int i;
    for (i = 1; i < 4; ++i) {
        if (cosines[i] > maxCosine) {
            maxCosine = cosines[i]; maxCosineIdx = i;
        }
    }
    if (maxCosineIdx == 1) {
        rotation = rotationX * rotation;
    } else if (maxCosineIdx == 2) {
        rotation = rotationZ * rotation;
    } else if (maxCosineIdx == 3) {
        rotation = (rotationZ * rotationX) * rotation;
    }
}

// Solve a cubic equation.
// Return: number of real roots found
//         x1, x2, x3 are real roots (if exist),
//         z1, z2, z3 are complex roots,
//             where z1 is always real, z1 == x1.
int SolveCubicEquation(
    double a, double b, double c, double d,
    double& x1, double& x2, double& x3,
    std::complex<double>& z1,
    std::complex<double>& z2,
    std::complex<double>& z3
) {
    if (a == 0.)
        return 0;

    // New variable y = x + b/(3a)
    double y1, y2, y3;
    double xyShift = b/(3.*a);
    double q = (
        2.*b*b*b/(27.*a*a*a) -
        b*c/(3.*a*a) +
        d/a
    ) / 2.;
    double p = (
        (3.*a*c - b*b)/(3.*a*a)
    ) / 3.;

    // Equivalent equation:
    // y^2 + 3*p*y + 2*q = 0

    // Discriminant:
    double D = q*q + p*p*p;
    if (D > 0) {
        // 1 real root and 2 imaginary
        // Using the Cardano formula
        D = sqrt(D);
        double u = (-q) + D;
        if (u > 0.) {
            u = pow(u, 1./3.);
        } else if (u < 0.) {
            u = (-pow(-u, 1./3.));
        }

        double v = (-q) - D;
        if (v > 0.) {
            v = pow(v, 1./3.);
        } else if (v < 0.) {
            v = (-pow(-v, 1./3.));
        }

        // Real root
        y1 = u + v;

        // Complex roots
        double sqrt3_2 = sqrt(3.)/2.;
        std::complex<double> i(0., 1.);
        std::complex<double> eps1 = (-0.5) + i*sqrt3_2;
        std::complex<double> eps2 = (-0.5) - i*sqrt3_2;
        std::complex<double> zy2 = (eps1*u) + (eps2*v);
        std::complex<double> zy3 = (eps2*u) + (eps1*v);

        x1 = y1 - xyShift;
        z1 = x1;
        z2 = zy2 - xyShift;
        z3 = zy3 - xyShift;

        // For equal real roots: D may be very small,
        // so this case will be considered as complex
        if (
            D <= R3_EPSILON &&
            fabs(z2.imag()) <= R3_EPSILON &&
            fabs(z3.imag()) <= R3_EPSILON
        ) {
            x2 = z2.real();
            x3 = z3.real();
            goto LSort;
        }

        return 1;
    } else if (D == 0. && p == 0.) {
        if (p == 0.) {
            // 3 zero roots
            y1 = 0.; y2 = 0.; y3 = 0.;
        }
    } else if (D <= 0.) {
        ASSERT(p < 0.);
        double r = sqrt(fabs(p));
        if (q < 0.) {
            r = (-r);
        }
        double cosPhi = q/(r*r*r);
        ASSERT(fabs(cosPhi) <= 1.);
        double phi = acos(cosPhi);
        y1 = (-2.*r*cos(phi/3.));
        y2 = (2.*r*cos(PI/3. - phi/3.));
        y3 = (2.*r*cos(PI/3. + phi/3.));
    }

    x1 = y1 - xyShift;
    x2 = y2 - xyShift;
    x3 = y3 - xyShift;
    z1 = x1; z2 = x2; z3 = x3;

    // Sort real roots in ascending order
    LSort: ;
    double tmp;
    std::complex<double> ztmp;
    if (x2 < x1) {
        tmp = x1; x1 = x2; x2 = tmp;
        ztmp = z1; z1 = z2; z2 = ztmp;
    }
    if (x3 < x1) {
        tmp = x1; x1 = x3; x3 = tmp;
        ztmp = z1; z1 = z3; z3 = ztmp;
    }
    if (x3 < x2) {
        tmp = x2; x2 = x3; x3 = tmp;
        ztmp = z2; z2 = z3; z3 = ztmp;
    }
    return 3;
}

void defineCharacteristicEquation(
    const R3Matrix& m,
    double coeff[4]     // Descending degree order
) {
    coeff[0] = 1.;
    coeff[1] = (-(m[0][0] + m[1][1] + m[2][2]));    // -trace
    coeff[2] = (-(      // Decompose determinant via a first row
        (-m[1][1]*m[2][2] + m[2][1]*m[1][2]) +
        m[0][0]*(-m[1][1] - m[2][2]) -
        m[0][1]*(-m[1][0]) +
        m[0][2]*m[2][0]
    ));
    coeff[3] = (-m.Determinant());
}

void computeEigenvalues(
    const R3Matrix& m,
    int& numReals,
    double realEigenvalues[3],      // Must be all different
    int& numComplex,
    std::complex<double> eigenValues[3]
) {
    double coeff[4];
    std::complex<double> z[3];
    defineCharacteristicEquation(m, coeff);
    numReals = SolveCubicEquation(
        coeff[0], coeff[1], coeff[2], coeff[3],
        realEigenvalues[0], realEigenvalues[1], realEigenvalues[2],
        z[0], z[1], z[2]
    );
    numComplex = 3;

    // Throw out repeated values
    int i = 0;
    while (i < numReals-1) {
        int j = i+1;
        while (j < numReals) {
            if (
                fabs(realEigenvalues[i] - realEigenvalues[j]) <=
                R3_EPSILON
            ) {
                // Remove j-th value
                int k;
                for (k = j; k < numReals-1; ++k) {
                    realEigenvalues[k] = realEigenvalues[k+1];
                    z[k] = z[k+1];
                }
                --numReals;
                --numComplex;
            } else {
                ++j;
            }
        } // end while (j...
        ++i;
    } // end while (i...
}

void GaussMethod(
    R3Matrix& a,
    int& rank
) {
    int n = 3;  // Number of columns
    int m = 3;  // Number of rows
    double r;
    int i, j, k, l;
    i = 0; j = 0;
    while (i < m && j < n) {
        // Looking up a maximal element in j-th column,
        // starting from i-th row
        r = fabs(a[i][j]); l = i;
        for (k = i + 1; k < m; ++k) {
            if (fabs(a[k][j]) > r) {
                l = k;      
                r = fabs(a[k][j]); 
            }
        }
        if (r <= R3_EPSILON) {
            // Make j-th column zero, starting from line i
            for (k = i; k < m; ++k) {
                a[k][j] = 0.0;
            }
            ++j;      // Go to the next column
            continue; // Next iteration
        }

        if (l != i) {
            // Swap rows i and l
            for (k = j; k < n; ++k) {
                r = a[i][k];
                a[i][k] = a[l][k];
                a[l][k] = (-r); // Change sign to preserve determinant
            }
        }

        ASSERT(fabs(a[i][j]) > R3_EPSILON);

        // Make j-th column zero, starting from row i+1
        for (k = i+1; k < m; ++k) {
            r = (-a[k][j] / a[i][j]);

            // Add i-th row, multiplied by r, to k-th row
            a[k][j] = 0.0;
            for (l = j+1; l < n; ++l) {
                a[k][l] += r * a[i][l];
            }
        }

        ++i; ++j;   // To the next minor of matrix
    }
    rank = i;
}

void computeRealEigenvectors(
    const R3Matrix& m,
    int numReals,
    const double realEigenvalues[3],
    int& numEigenvectors,
    R3Vector realEigenvectors[3],
    double eigenvalues[3]
) {
    R3Vector eigenvectors[3];
    int vectorsFound = 0;

    int i, j;
    for (i = 0; i < numReals && vectorsFound < 3; ++i) {
        double lambda = realEigenvalues[i];
        // Skip eigenvalue, if it was considered at previous steps
        int wasProcessed = 0;
        for (j = 0; j < i; ++j) {
            if (fabs(lambda - realEigenvalues[j]) <= R3_EPSILON) {
                wasProcessed = 1; break;
            }
        }
        if (wasProcessed)
            continue;

        R3Matrix a(m);
        a[0][0] -= lambda;
        a[1][1] -= lambda;
        a[2][2] -= lambda;
        int dimension;
        solveHomogeniousEquation(
            a, dimension, eigenvectors
        );
        for (j = 0; j < dimension && vectorsFound < 3; ++j) {
            realEigenvectors[vectorsFound] =
                eigenvectors[j];
            eigenvalues[vectorsFound] = lambda;
            ++vectorsFound;
        }
    }
    numEigenvectors = vectorsFound;
}

void solveHomogeniousEquation(
    R3Matrix& a,
    int& dimension,             // Dimension of solution subspace
    R3Vector solutionBasis[3]   // Basis of solution subspace
) {
    int rank;
    int i, j;
    int n = 3;  // Number of columns
    int m = 3;  // Number of rows
    int found;

    GaussMethod(a, rank);
    if (rank == 3) {
        dimension = 0;
        return;
    } else if (rank == 2) {
        dimension = 1;
        // Look up first non-zero element in row 1
        i = 1;
        found = 0;
        for (j = 1; j < n; ++j) {
            if (fabs(a[i][j]) >= R3_EPSILON) {
                found = 1; break;
            }
        }
        ASSERT(found);
        if (!found) {   // Should not be so
            dimension = 0; return;
        }
        int nonZero1 = j;
        if (nonZero1 == 1) {
            // Free variable z
            solutionBasis[0].z = 1.0;
            solutionBasis[0].y = (
                -a[i][2]*solutionBasis[0].z / a[i][1]
            );
            ASSERT(fabs(a[0][0]) >= R3_EPSILON);
            solutionBasis[0].x = (
                -(
                    a[0][1]*solutionBasis[0].y +
                    a[0][2]*solutionBasis[0].z
                ) / a[0][0]
            );
            return;
        } else {
            ASSERT(nonZero1 == 2);
            solutionBasis[0].z = 0.0;

            // Look up first non-zero element in row 0
            i = 0;
            found = 0;
            for (j = 0; j < n-1; ++j) {
                if (fabs(a[i][j]) >= R3_EPSILON) {
                    found = 1; break;
                }
            }
            ASSERT(found);
            if (!found) {   // Should not be so
                dimension = 0; return;
            }
            int nonZero0 = j;
            if (nonZero0 == 0) {
                // Free variable y
                solutionBasis[0].y = 1.0;
                solutionBasis[0].x = (
                    -(
                        a[0][1]*solutionBasis[0].y +
                        a[0][2]*solutionBasis[0].z      // Actually, zero
                    ) / a[0][0]
                );
                return;
            } else {
                // Free variable x
                solutionBasis[0].y = 0.0;
                solutionBasis[0].x = 1.0;
                return;
            }
        }
    } else if (rank == 1) {
        dimension = 2;
        // Look up first non-zero element in row 0
        i = 0;
        found = 0;
        for (j = 0; j < n; ++j) {
            if (fabs(a[i][j]) >= R3_EPSILON) {
                found = 1; break;
            }
        }
        ASSERT(found);
        if (!found) {   // Should not be so
            dimension = 0; return;
        }
        int nonZero0 = j;
        if (nonZero0 == 0) {
            // Free variables y, z
            solutionBasis[0].y = 1.0;
            solutionBasis[0].z = 0.0;
            solutionBasis[0].x = (
                -(
                    a[0][1]*solutionBasis[0].y +
                    a[0][2]*solutionBasis[0].z
                ) / a[0][0]
            );

            solutionBasis[1].y = 0.0;
            solutionBasis[1].z = 1.0;
            solutionBasis[1].x = (
                -(
                    a[0][1]*solutionBasis[1].y +
                    a[0][2]*solutionBasis[1].z
                ) / a[0][0]
            );
            return;
        } else if (nonZero0 == 1) {
            // Free variables x, z
            solutionBasis[0].x = 1.0;
            solutionBasis[0].z = 0.0;
            solutionBasis[0].y = 0.0;

            solutionBasis[1].x = 0.0;
            solutionBasis[1].z = 1.0;
            solutionBasis[1].y = (
                -(
                    a[0][2]*solutionBasis[1].z
                ) / a[0][1]
            );
            return;
        } else {
            ASSERT(nonZero0 == 2);

            // Free variables x, y
            solutionBasis[0].x = 1.0;
            solutionBasis[0].y = 0.0;
            solutionBasis[0].z = 0.0;

            solutionBasis[1].x = 0.0;
            solutionBasis[1].y = 1.0;
            solutionBasis[1].z = 0.0;
            return;
        }
    }
}

void computeTransitionMatrix(
    const R3Vector u[3],    // Orthonormal bases
    const R3Vector v[3],    // m: u_i -> v_i
    R3Matrix& m
) {
    // Standard basis
    R3Vector e[3];
    e[0] = R3Vector(1., 0., 0.);
    e[1] = R3Vector(0., 1., 0.);
    e[2] = R3Vector(0., 0., 1.);

    R3Matrix alpha, beta;
    // Express basis vectors via u_i
    // e_i = \sum_i alpha_{ij} u_j
    alpha[0][0] = e[0] * u[0];
    alpha[0][1] = e[0] * u[1];
    alpha[0][2] = e[0] * u[2];

    alpha[1][0] = e[1] * u[0];
    alpha[1][1] = e[1] * u[1];
    alpha[1][2] = e[1] * u[2];

    alpha[2][0] = e[2] * u[0];
    alpha[2][1] = e[2] * u[1];
    alpha[2][2] = e[2] * u[2];

    // Express v_i via basis vectors
    // v_j = \sum_k beta_{jk} e_k
    beta[0][0] = v[0] * e[0];
    beta[0][1] = v[0] * e[1];
    beta[0][2] = v[0] * e[2];

    beta[1][0] = v[1] * e[0];
    beta[1][1] = v[1] * e[1];
    beta[1][2] = v[1] * e[2];

    beta[2][0] = v[2] * e[0];
    beta[2][1] = v[2] * e[1];
    beta[2][2] = v[2] * e[2];

    // m * e_i = \sum_j alpha_{ij} v_j  =
    //         = \sum_j alpha_{ij} \sum_k beta_{jk} e_k
    // Let gamma = alpha * beta
    // Rows of gamma are images of basis vectors
    // Copy them to columns of resulting matrix
    R3Matrix gamma = alpha * beta;
    m[0][0] = gamma[0][0];
    m[1][0] = gamma[0][1];
    m[2][0] = gamma[0][2];

    m[0][1] = gamma[1][0];
    m[1][1] = gamma[1][1];
    m[2][1] = gamma[1][2];

    m[0][2] = gamma[2][0];
    m[1][2] = gamma[2][1];
    m[2][2] = gamma[2][2];
}

void start3DInertionMatrixDefinition(
    R3Matrix& inertionMatrix
) {
    inertionMatrix.erase();
}

void finish3DInertionMatrixDefinition(
    const R3Matrix& inertionMatrix
) {
    // Nothing to do.
}

void addPointTo3DInertionMatrix(
    const R3Vector& point,
    double pointMass,
    const R3Vector& center,
    R3Matrix& inertionMatrix
) {
    R3Vector dv = point - center;

    double Jxx = inertionMatrix[0][0];
    double Jyy = inertionMatrix[1][1];
    double Jzz = inertionMatrix[2][2];
    double Jxy = (-inertionMatrix[0][1]);
    double Jxz = (-inertionMatrix[0][2]);
    double Jyz = (-inertionMatrix[1][2]);
    double w = pointMass;

    Jxx += (dv.y*dv.y + dv.z*dv.z)*w;
    Jyy += (dv.x*dv.x + dv.z*dv.z)*w;
    Jzz += (dv.x*dv.x + dv.y*dv.y)*w;
    Jxy += dv.x*dv.y*w;
    Jxz += dv.x*dv.z*w;
    Jyz += dv.y*dv.z*w;

    inertionMatrix[0][0] = Jxx;
    inertionMatrix[1][1] = Jyy;
    inertionMatrix[2][2] = Jzz;
    inertionMatrix[0][1] = (-Jxy);
    inertionMatrix[1][0] = (-Jxy);
    inertionMatrix[0][2] = (-Jxz);
    inertionMatrix[2][0] = (-Jxz);
    inertionMatrix[1][2] = (-Jyz);
    inertionMatrix[2][1] = (-Jyz);
}

static const double MATRIX_NORMALIZATION = 1.;
// Multiply the values of matrix elements so
// that all of them be in the range 0. - 1.
void normalize3DInertionMatrix(
    R3Matrix& inertionMatrix
) {
    double maxElem = 0.;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (fabs(inertionMatrix[i][j]) > maxElem)
                maxElem = fabs(inertionMatrix[i][j]);
        }
    }

    if (maxElem > R3_EPSILON) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                inertionMatrix[i][j] /= maxElem;
            }
        }
    }
}
