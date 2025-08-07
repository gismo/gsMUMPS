#pragma once

#include <gsCore/gsLinearAlgebra.h>

// Include your MUMPS support classes
#include <gsMUMPS/src/MUMPSSupport.h>
#include <gsMatrix/gsSparseSolver.h>

namespace gismo
{
    // You can add convenience typedefs here if needed
    template<typename MatrixType = gsSparseMatrix<>>
    using gsMUMPSLU = gsEigen::MUMPSLU<MatrixType>;

    template<typename MatrixType = gsSparseMatrix<>, int UpLo = gsEigen::Lower>
    using gsMUMPSLDLT = gsEigen::MUMPSLDLT<MatrixType, UpLo>;

    // template<typename T> class gsEigenMUMPSLU;
    // template<typename T> class gsEigenMUMPSLDLT;

    // GISMO_EIGEN_SPARSE_SOLVER(gsEigenMUMPSLU, MUMPSLU)
    // GISMO_EIGEN_SPARSE_SOLVER(gsEigenMUMPSLDLT, MUMPSLDLT)
}