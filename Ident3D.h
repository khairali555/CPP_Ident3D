#ifndef IDENT3D_H
#define IDENT3D_H

#include "MFC.h"        // Fight againts Microsoft

#include <complex>
#include "matrix.h"

void define3DTransform(
    int numLandmarks0,
    const R3Vector* landmarkPos0,
    const double* landmarkWeight0,
    int numLandmarks1,
    const R3Vector* landmarkPos1,
    const double* landmarkWeight1,
    R3Matrix& rotation,
    R3Vector& shift
);

void define3DCentroid(
    int numLandmarks,
    const R3Vector* landmarkPos,
    const double* landmarkWeight,
    R3Vector& center,
    double& mass
);

void define3DInertionMatrix(
    int numLandmarks,
    const R3Vector* landmarkPos,
    const double* landmarkWeight,
    const R3Vector& center,
    R3Matrix& inertionMatrix,
    R3Vector& secondMoment
);

void start3DInertionMatrixDefinition(
    R3Matrix& inertionMatrix
);

void finish3DInertionMatrixDefinition(
    const R3Matrix& inertionMatrix
);

void normalize3DInertionMatrix(
    R3Matrix& inertionMatrix
);

void addPointTo3DInertionMatrix(
    const R3Vector& point,
    double pointMass,
    const R3Vector& center,
    R3Matrix& inertionMatrix
);

void defineCharacteristicEquation(
    const R3Matrix& m,
    double coeff[4]     // Descending degree order
);

void computeEigenvalues(
    const R3Matrix& m,
    int& numReals,
    double realEigenvalues[3],      // Must be all different
    int& numComplex,
    std::complex<double> complexEigenvalues[3]
);

void computeRealEigenvectors(
    const R3Matrix& m,
    int numEigenvalues,             // Must be all different
    const double realEigenvalues[3],
    int& numEigenvectors,           // Number of eigenvectors, that are
    R3Vector realEigenvectors[3],   // linearly independent
    double eigenvalues[3]           // Their values
);

void define3DInertionAxes(
    const R3Matrix& inertionMatrix,
    double eigenvalues[3],
    R3Vector axes[3]
);

void define3DRotation(
    const double eigenvalues0[3],
    const R3Vector axes0[3],
    const R3Vector& secondMoment0,
    const double eigenvalues1[3],
    const R3Vector axes1[3],
    const R3Vector& secondMoment1,
    R3Matrix& rotation
);

// Solve a quibic equation.
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
);

void computeRealEigenvectors(
    const R3Matrix& m,
    int numReals,
    const double realEigenvalues[3],
    R3Vector realEigenvectors[3]
);

void GaussMethod(
    R3Matrix& m,
    int& rank
);

void solveHomogeniousEquation(
    R3Matrix& a,
    int& dimension,             // Dimension of solution subspace
    R3Vector solutionBasis[3]   // Basis of solution subspace
);

void computeTransitionMatrix(
    const R3Vector basis0[3],  // Orthonormal bases
    const R3Vector basis1[3],
    R3Matrix& m
);

#endif
