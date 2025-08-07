// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2018 <darcy.beurle@ibnm.uni-hannover.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_MUMPSSUPPORT_H
#define EIGEN_MUMPSSUPPORT_H
#define Eigen gsEigen
#define eigen_assert( cond ) GISMO_ASSERT( cond, "" )

#include "Eigen/Sparse"

#include <dmumps_c.h>
#include <smumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>

#include <complex>
#include <stdexcept>
#include <iostream>

namespace Eigen
{
#if defined(DCOMPLEX)
#define MUMPS_COMPLEX COMPLEX
#define MUMPS_DCOMPLEX DCOMPLEX
#else
#define MUMPS_COMPLEX std::complex<float>
#define MUMPS_DCOMPLEX std::complex<double>
#endif

/** \ingroup MUMPSSupport_Module
 * \brief Interface to the MUMPS solver
 *
 * This class is used to solve the linear systems A.X = B via the MUMPS library
 * The matrix can be either real or complex, symmetric or unsymmetric.
 *
 * \sa TutorialSparseDirectSolvers
 */
template <typename MatrixT>
class MUMPSLU;

template <typename MatrixT, int Options>
class MUMPSLDLT;

namespace internal
{
template <class Mumps>
struct mumps_traits;

template <typename MatrixT>
struct mumps_traits<MUMPSLU<MatrixT>>
{
    typedef MatrixT MatrixType;
    typedef typename MatrixT::Scalar Scalar;
    typedef typename MatrixT::RealScalar RealScalar;
    typedef typename MatrixT::StorageIndex StorageIndex;
};

template <typename MatrixT, int Options>
struct mumps_traits<MUMPSLDLT<MatrixT, Options>>
{
    typedef MatrixT MatrixType;
    typedef typename MatrixT::Scalar Scalar;
    typedef typename MatrixT::RealScalar RealScalar;
    typedef typename MatrixT::StorageIndex StorageIndex;
};

template <typename T>
struct MUMPSAPIWrapper;

template <>
struct MUMPSAPIWrapper<float>
{
    using MUMPS_STRUC_C = SMUMPS_STRUC_C;

    static void mumps_c(MUMPS_STRUC_C& info) { smumps_c(&info); }
};

template <>
struct MUMPSAPIWrapper<double>
{
    using MUMPS_STRUC_C = DMUMPS_STRUC_C;

    static void mumps_c(MUMPS_STRUC_C& info) { dmumps_c(&info); }
};

template <>
struct MUMPSAPIWrapper<std::complex<float>>
{
    using MUMPS_STRUC_C = CMUMPS_STRUC_C;

    static void mumps_c(MUMPS_STRUC_C& info) { cmumps_c(&info); }
};

template <>
struct MUMPSAPIWrapper<std::complex<double>>
{
    using MUMPS_STRUC_C = ZMUMPS_STRUC_C;

    static void mumps_c(MUMPS_STRUC_C& info) { zmumps_c(&info); }
};

namespace mumps
{
enum ordering { AMD, AMF = 2, Scotch, Pord, Metis, QAMD, automatic };

// Jobs in MUMPS use the following:
// 4   Job = 1 && Job = 2
// 5   Job = 2 && Job = 3
// 6   Job = 1 && Job = 2 && Job = 3
enum job {
    terminate = -2,
    initialisation = -1,
    analysis = 1,
    factorisation = 2,
    back_substitution = 3
};

enum residual { none, expensive, cheap };

enum matrix_property { unsymmetric, SPD, general_symmetric };
}
}

// Base class to interface with MUMPS. Users should not instantiate this class directly.
template <class Derived>
class MumpsBase : public SparseSolverBase<Derived>
{
protected:
    typedef SparseSolverBase<Derived> Base;
    using Base::derived;
    using Base::m_isInitialized;

public:
    using Base::_solve_impl;

    typedef typename internal::mumps_traits<Derived>::MatrixType MatrixT;
    typedef MatrixT MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename MatrixType::StorageIndex StorageIndex;
    typedef Matrix<Scalar, Dynamic, 1> Vector;
    typedef SparseMatrix<Scalar, ColMajor> ColSpMatrix;
    enum {
        ColsAtCompileTime = MatrixType::ColsAtCompileTime,
        MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

    using solver_api_type = internal::MUMPSAPIWrapper<Scalar>;

public:
    MumpsBase() : m_initisOk(false), m_analysisIsOk(false), m_factorizationIsOk(false), m_size(0) {}

    ~MumpsBase()
    {
        m_solver.job = internal::mumps::job::terminate;
        solver_api_type::mumps_c(m_solver);
    }

    /* Solve the system */
    template <typename Rhs, typename Dest>
    bool _solve_impl(const MatrixBase<Rhs>& b, MatrixBase<Dest>& x) const;

    Index cols() const { return m_size; }
    Index rows() const { return m_size; }

    /**
     * \brief Reports whether previous computation was successful.
     *
     * \returns \c Success if computation was succesful,
     *          \c NumericalIssue if the MUMPS reports a problem
     *          \c InvalidInput if the input matrix is invalid
     *
     * \sa iparm()
     */
    ComputationInfo info() const
    {
        eigen_assert(m_isInitialized && "Decomposition is not initialized.");
        return m_info;
    }

    /**
     * \brief Access to the MUMPS integer parameters.
     *
     * \note Follows the 1-based numbering in the manual
     *
     * \returns a reference to the MUMPS integer parameter at index \c i.
     */
    int & ICNTL(int i)
    {
        eigen_assert(m_isInitialized && "Decomposition is not initialized.");
        eigen_assert(i > 0 && i <= sizeof(m_solver.icntl)/sizeof(m_solver.icntl[0]) && "icntl index out of bounds (1 to "<<sizeof(m_solver.icntl)/sizeof(m_solver.icntl[0])<<")");
        return m_solver.icntl[i-1];
    }
    const int & ICNTL(int i) const
    {
        eigen_assert(m_isInitialized && "Decomposition is not initialized.");
        eigen_assert(i > 0 && i <= sizeof(m_solver.icntl)/sizeof(m_solver.icntl[0]) && "icntl index out of bounds (1 to "<<sizeof(m_solver.icntl)/sizeof(m_solver.icntl[0])<<")");
        return m_solver.icntl[i-1];
    }

protected:
    // Initialize the MUMPS data structure, check the matrix
    void init(bool const is_matrix_symmetric);

    // Compute the ordering and the symbolic factorization
    void analyzePattern();

    // Compute the numerical factorization
    void factorize();

    void compute();

protected:
    int m_initisOk;
    int m_analysisIsOk;
    int m_factorizationIsOk;
    mutable ComputationInfo m_info;

    typename internal::MUMPSAPIWrapper<Scalar>::MUMPS_STRUC_C mutable m_solver;

    mutable MUMPS_INT m_size; // Size of the matrix

    mutable Array<MUMPS_INT, Dynamic, 1> m_rows, m_cols; // Coordinate format
    mutable Vector m_coeffs;                             // Coefficients in sparse matrix
};

/**
 * Initialize the MUMPS data structure.
 * A first call to this function fills iparm and dparm with the default MUMPS parameters
 * \sa iparm() dparm()
 */
template <class Derived>
void MumpsBase<Derived>::init(bool const is_matrix_symmetric)
{
    m_solver.job = internal::mumps::job::initialisation;

    // Par determines parallelisation of the computation
    // 0  : Host is not involved with computation
    // 1  : Host is involved
    m_solver.par = 1;

    // Sym is 1 if the matrix is symmetric (triangular) and 0 if the matrix
    // is full.  This determines if LDLT or LU factorisation is used
    // respectively
    m_solver.sym = is_matrix_symmetric;

    // No distributed memory MPI Communicator for horrible FORTRAN compatibility
    m_solver.comm_fortran = -987654;

    solver_api_type::mumps_c(m_solver);

    // ICNTL VALUES, BASED ON userguide_5.8.1.pdf
    // ICNTL1 - the output stream for error messages
    m_solver.icntl[0] = 6; // standard output stream
    // ICNTL2 - the output stream for diagnostic printing and statistics local to each MPI process
    m_solver.icntl[1] = 0; // suppress messages
    // ICNTL3 - the output stream for global information, collected on the host
    m_solver.icntl[2] = 6; // standard output stream
    // ICNTL4 - the level of printing for error, warning, and diagnostic messages
    m_solver.icntl[3] = 2; // errors, warnings and main statistics printed
    // ICNTL5 - controls the matrix input format (see Subsection 5.4.2)
    m_solver.icntl[4] = 0; // assembled format
    // ICNTL6 - permutes the matrix to a zero-free diagonal and/or scale the matrix (see Subsection 3.2 and Subsection 5.5.2).
    m_solver.icntl[5] = 7; // automatic choice done by the package
    // ICNTL7 - computes a symmetric permutation (ordering) to determine the pivot order to be used for the factorization in case of sequential analysis (ICNTL(28)=1). See Subsection 3.2 and Subsection 5.6.
    m_solver.icntl[6] = internal::mumps::ordering::automatic; // automatic ordering
    // ICNTL8 - describes the scaling strategy (see Subsection 5.5)
    m_solver.icntl[7] = 77; // automatic scaling
    // ICNTL9 - computes the solution using A or A^T
    m_solver.icntl[8] = 1; // AX=B is solved
    // ICNTL10 - applies the iterative refinement to the computed solution (see Subsection 5.8).
    m_solver.icntl[9] = 0;
    // ICNTL11 -  computes statistics related to an error analysis of the linear system solved (Ax = b or AT x = b (see ICNTL(9))). See Subsection 5.9
    m_solver.icntl[10] = internal::mumps::residual::none;
    // ICNTL12 - defines an ordering strategy for symmetric matrices (SYM = 2) (see [25] for more details) and is used, in conjunction with ICNTL(6), to add constraints to the ordering algorithm (ICNTL(7) option)
    m_solver.icntl[11] = 0; // automatic
    // ICNTL13 - controls the parallelism of the root node (valid only when multiple processors are being used)
    m_solver.icntl[12] = 0; // parallel factorization on the root node
    // ICNTL14 - controls the percentage increase in the estimated working space
    m_solver.icntl[13] = 20; // 20% increase
    // ICNTL15 - exploits compression of the input matrix resulting from a block format, see Subsection 5.7.
    m_solver.icntl[14] = 0; // no compression
    // ICNTL16 - controls the setting of the number of OpenMP threads, see Subsection 5.21 by MUMPS when the setting of multithreading is not possible outside MUMPS (see Subsection 3.13).
    m_solver.icntl[15] = 0; // nothing is done
    // ICNTL17 - reserved for future use
    // m_solver.icntl[16] = 0;
    // ICNTL18 - defines the strategy for the distributed input matrix (only for assembled matrix, see Subsection 5.4.2).
    m_solver.icntl[17] = 0; // input matrix centralized on the host
    // ICNTL19 - computes the Schur complement matrix (see Subsection 5.18).
    m_solver.icntl[18] = 0; // complete factorization
    // ICNTL20 - determines the format (dense, sparse, or distributed) of the right-hand sides
    m_solver.icntl[19] = 0; // dense right-hand side
    // ICNTL21 - determines the distribution (centralized or distributed) of the solution vectors.
    m_solver.icntl[20] = 0; // assembled centralized format
    // ICNTL22 - controls the in-core/out-of-core (OOC) factorization and solve.
    m_solver.icntl[21] = 0; // in-core factorization
    // ICNTL23 - corresponds to the maximum size of the working memory in MegaBytes that MUMPS can allocate per working process, see Subsection 5.11 for more details.
    m_solver.icntl[22] = 0; // each processor will allocate workspace based on the estimates computed during the analysis
    // ICNTL24 - controls the detection of null pivot rows
    m_solver.icntl[23] = 0; // null pivot row detection disabled
    // ICNTL25 - controls the computation of a null space basis
    m_solver.icntl[24] = 0; // null space basis not computed
    // ICNTL26 -  drives the solution phase if a Schur complement matrix has been computed (ICNTL(19)̸ =0), see Subsection 3.18 for details
    m_solver.icntl[25] = 0; // normal solution phase
    // ICNTL27 - controls the blocking size for multiple right-hand sides
    m_solver.icntl[26] = -32; // automatic choice
    // ICNTL28 - controls the ordering strategy
    m_solver.icntl[27] = 0; // automatic choice between sequential and parallel ordering
    // ICNTL29 - controls the parallel ordering tool
    m_solver.icntl[28] = 0; // automatic choice of parallel ordering tool
    // ICNTL30 - controls the computation of some entries of the inverse
    m_solver.icntl[29] = 0; // inverse entries not computed
    // ICNTL31 - controls which factors may be discarded during factorization
    m_solver.icntl[30] = 0; // factors kept during factorization
    // ICNTL32 - performs forward elimination during factorization
    m_solver.icntl[31] = 0; // standard factorization
    // ICNTL33 - computes the determinant of the input matrix
    m_solver.icntl[32] = 0; // determinant not computed
    // ICNTL34 - controls conservation of OOC files during save/restore
    m_solver.icntl[33] = 0; // out-of-core files marked for deletion
    // ICNTL35 - controls activation of BLR feature
    m_solver.icntl[34] = 0; // standard multifrontal factorization
    // ICNTL36 - controls choice of BLR factorization variant
    m_solver.icntl[35] = 0; // UFSC variant
    // ICNTL37 - controls BLR compression of contribution blocks
    m_solver.icntl[36] = 0; // contribution blocks not compressed
    // ICNTL38 - estimates compression rate of LU factors
    m_solver.icntl[37] = 600; // 60% compression rate estimate
    // ICNTL39 - estimates compression rate of contribution blocks
    m_solver.icntl[38] = 500; // 50% compression rate estimate
    // ICNTL40-47 - reserved in current version
    // m_solver.icntl[39] = 0;
    // m_solver.icntl[40] = 0;
    // m_solver.icntl[41] = 0;
    // m_solver.icntl[42] = 0;
    // m_solver.icntl[43] = 0;
    // m_solver.icntl[44] = 0;
    // m_solver.icntl[45] = 0;
    // m_solver.icntl[46] = 0;
    // ICNTL48 - multithreading with tree parallelism
    m_solver.icntl[47] = 1; // multithreaded tree parallelism activated
    // ICNTL49 - compact workarray at end of factorization
    m_solver.icntl[48] = 0; // nothing done
    // ICNTL50-55 - reserved in current version
    // m_solver.icntl[49] = 0;
    // m_solver.icntl[50] = 0;
    // m_solver.icntl[51] = 0;
    // m_solver.icntl[52] = 0;
    // m_solver.icntl[53] = 0;
    // m_solver.icntl[54] = 0;
    // ICNTL56 - detects pseudo-singularities and uses rank-revealing factorization
    m_solver.icntl[55] = 0; // standard factorization
    // ICNTL57 - reserved in current version
    // m_solver.icntl[56] = 0;
    // ICNTL58 - defines options for symbolic factorization
    m_solver.icntl[57] = 2; // column count based symbolic factorization
    // ICNTL59-60 - not used in current version
    // m_solver.icntl[58] = 0;
    // m_solver.icntl[59] = 0;

    m_size = 0;

    // Check the returned error
    if (m_solver.info[0] < 0)
    {
        m_info = InvalidInput;
        m_initisOk = false;
        throw std::domain_error("Error code " + std::to_string(m_solver.info[0])
                                + " in MUMPS initialization");
    }
    else
    {
        m_info = Success;
        m_initisOk = true;
    }
}

template <class Derived>
void MumpsBase<Derived>::compute()
{
    eigen_assert(m_rows.size() == m_cols.size() && "The input matrix should be square");
    analyzePattern();
    factorize();
}

template <class Derived>
void MumpsBase<Derived>::analyzePattern()
{
    eigen_assert(m_initisOk && "The initialization of MUMPS failed");

    m_solver.job = internal::mumps::job::analysis;

    m_solver.n = m_size;
    m_solver.nz = internal::convert_index<MUMPS_INT8>(m_cols.size());

    m_solver.a = 0;
    m_solver.rhs = 0;

    m_solver.irn = m_rows.data();
    m_solver.jcn = m_cols.data();

    solver_api_type::mumps_c(m_solver);

    // Check the returned error
    if (m_solver.info[0] < 0)
    {
        m_info = NumericalIssue;
        m_analysisIsOk = false;
        throw std::domain_error("Error code " + std::to_string(m_solver.info[0])
                                + " in MUMPS analyzePattern()");
    }
    else
    {
        m_info = Success;
        m_analysisIsOk = true;
    }
}

template <class Derived>
void MumpsBase<Derived>::factorize()
{
    eigen_assert(m_analysisIsOk && "analysePattern() should be called before factorize()");
    eigen_assert(m_cols.size() == m_rows.size() && "Row and column index sizes must match");

    m_solver.n = m_size;
    m_solver.nz = internal::convert_index<MUMPS_INT8>(m_cols.size());

    m_solver.a = m_coeffs.data();
    m_solver.irn = m_rows.data();
    m_solver.jcn = m_cols.data();

    m_solver.job = internal::mumps::job::factorisation;

    solver_api_type::mumps_c(m_solver);

    // Check the returned error
    if (m_solver.info[0] < 0)
    {
        m_info = NumericalIssue;
        m_factorizationIsOk = false;
        m_isInitialized = false;
        throw std::domain_error("Error code " + std::to_string(m_solver.info[0])
                                + " in MUMPS factorize()");
    }
    else
    {
        m_info = Success;
        m_factorizationIsOk = true;
        m_isInitialized = true;
    }
}

template <typename Base>
template <typename Rhs, typename Dest>
bool MumpsBase<Base>::_solve_impl(const MatrixBase<Rhs>& b, MatrixBase<Dest>& x) const
{
    eigen_assert(m_isInitialized && "Call factorize() first");

    EIGEN_STATIC_ASSERT((Dest::Flags & RowMajorBit) == 0,
                        THIS_METHOD_IS_ONLY_FOR_COLUMN_MAJOR_MATRICES);

    // on return, x is overwritten by the computed solution
    x = b;

    m_solver.n = m_size;
    m_solver.nz = internal::convert_index<MUMPS_INT8>(m_cols.size());

    m_solver.a = m_coeffs.data();
    m_solver.irn = m_rows.data();
    m_solver.jcn = m_cols.data();

    m_solver.job = internal::mumps::job::back_substitution;

    m_solver.nrhs = 1;
    m_solver.lrhs = m_size;

    for (Index i = 0; i < b.cols(); i++)
    {
        m_solver.rhs = &x(0, i);
        solver_api_type::mumps_c(m_solver);
    }

    // Check the returned error
    m_info = m_solver.info[0] < 0 ? NumericalIssue : Success;

    if (m_info == NumericalIssue)
    {
        throw std::domain_error("Error code " + std::to_string(m_solver.info[0]) + " in MUMPS solve");
    }
    return m_solver.info[0] >= 0;
}

/**
 * \ingroup MUMPSSupport_Module
 * \class MUMPSLU
 * \brief Sparse direct LU solver based on MUMPS library
 *
 * This class is used to solve the linear systems A.X = B with a multifrontal LU
 * factorization in the MUMPS library. The matrix A should be square and nonsingular
 * MUMPS requires that the matrix A has a symmetric structural pattern.
 * This interface can symmetrize the input matrix otherwise.
 * The vectors or matrices X and B can be either dense or sparse.
 *
 * \tparam MatrixT the type of the sparse matrix A, it must be a SparseMatrix<>
 *
 * \implsparsesolverconcept
 *
 * \sa \ref TutorialSparseSolverConcept, class SparseLU
 */
template <typename MatrixT>
class MUMPSLU : public MumpsBase<MUMPSLU<MatrixT>>
{
public:
    typedef MatrixT MatrixType;
    typedef MumpsBase<MUMPSLU<MatrixType>> Base;
    typedef typename Base::ColSpMatrix ColSpMatrix;
    typedef typename MatrixType::StorageIndex StorageIndex;

public:

    MUMPSLU() : Base() { init(true); }

    explicit MUMPSLU(const MatrixType& matrix) : Base()
    {
        init(true);
        compute(matrix);
    }

    /**
     * Compute the LU multifrontal factorization of \p matrix.
     * \sa analyzePattern() factorize()
     */
    void compute(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::compute();
    }

    /**
     * Compute the LU symbolic factorization of \p matrix using its sparsity pattern.
     * Several ordering methods can be used at this step. See the MUMPS user's manual.
     * The result of this operation can be used with successive matrices having the same pattern as
     * \p matrix \sa factorize()
     */
    void analyzePattern(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::analyzePattern();
    }

    /**
     * Compute the LU multifrontal factorization of \p matrix
     * WARNING The matrix \p matrix should have the same structural pattern
     * as the same used in the analysis phase.
     * \sa analyzePattern()
     */
    void factorize(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::factorize();
    }

protected:
    void allocate_coordinate_format(const MatrixType& matrix)
    {
        m_size = matrix.rows();

        m_rows.resize(matrix.nonZeros());
        m_cols.resize(matrix.nonZeros());

        m_coeffs = matrix.coeffs();

        // Decompress the upper part of the sparse matrix and convert to one
        // based indexing for MUMPS solver
        for (Index k = 0, l = 0; k < matrix.outerSize(); ++k)
        {
            for (typename MatrixType::InnerIterator it(matrix, k); it; ++it, ++l)
            {
                m_rows(l) = it.row();
                m_cols(l) = it.col();
            }
        }
        // Convert to one based indexing
        m_rows += 1;
        m_cols += 1;
    }

protected:
    using Base::init;

    using Base::m_coeffs;
    using Base::m_cols;
    using Base::m_rows;

    using Base::m_size;
};

/**
 * \ingroup MUMPSSupport_Module
 * \class MUMPSLDLT
 * \brief A sparse direct multifrontal Cholesky (LDLT) factorization and solver
 * based on the MUMPS library
 *
 * This class is used to solve the linear systems A.X = B via a LDL^T multifrontal
 * factorization available in the MUMPS library. The matrix A should be symmetric
 * and positive definite
 * WARNING Selfadjoint complex matrices are not supported in the current version
 * of MUMPS.  The vectors or matrices X and B can be either dense or sparse
 *
 * \tparam MatrixT the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo The part of the matrix to use : Lower or Upper. The default is
 * Lower as required by MUMPS
 *
 * \implsparsesolverconcept
 *
 * \sa \ref TutorialSparseSolverConcept, class SimplicialLDLT
 */
template <typename MatrixT, int UpLoT>
class MUMPSLDLT : public MumpsBase<MUMPSLDLT<MatrixT, UpLoT>>
{
public:
    typedef MatrixT MatrixType;
    typedef MumpsBase<MUMPSLDLT<MatrixType, UpLoT>> Base;
    typedef typename Base::ColSpMatrix ColSpMatrix;

public:
    enum { UpLo = UpLoT };
    MUMPSLDLT() : Base() { init(true); }

    explicit MUMPSLDLT(const MatrixType& matrix) : Base()
    {
        init(true);
        compute(matrix);
    }

    /**
     * Compute the L and D factors of the LDL^T factorization of \p matrix
     * \sa analyzePattern() factorize()
     */
    void compute(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::compute();
    }

    /**
     * Compute the LDL^T symbolic factorization of \p matrix using its sparsity
     * pattern. The result of this operation can be used with successive matrices
     * with the same pattern as \p matrix \sa factorize()
     */
    void analyzePattern(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::analyzePattern();
    }

    /** Compute the LDL^T multifrontal numerical factorization of \p matrix */
    void factorize(const MatrixType& matrix)
    {
        allocate_coordinate_format(matrix);
        Base::factorize();
    }

protected:
    void allocate_coordinate_format(const MatrixType& matrix)
    {
        eigen_assert(matrix.rows() == matrix.cols() && "Input matrix must be square");

        m_size = matrix.rows();

        Index sym_nonzeros = 0;
        for (Index k = 0; k < matrix.outerSize(); ++k)
        {
            for (typename MatrixType::InnerIterator it(matrix, k); it; ++it)
            {
                if (it.col() >= it.row()) ++sym_nonzeros;
            }
        }

        m_rows.resize(sym_nonzeros);
        m_cols.resize(sym_nonzeros);
        m_coeffs.resize(sym_nonzeros);

        // Decompress the upper part of the sparse matrix and convert to one
        // based indexing for MUMPS solver
        for (Index k = 0, l = 0; k < matrix.outerSize(); ++k)
        {
            for (typename MatrixType::InnerIterator it(matrix, k); it; ++it)
            {
                if (it.col() >= it.row())
                {
                    m_rows(l) = it.row();
                    m_cols(l) = it.col();
                    m_coeffs(l) = it.value();
                    ++l;
                }
            }
        }
        m_rows += 1;
        m_cols += 1;
    }

protected:
    using Base::init;

    using Base::m_coeffs;
    using Base::m_cols;
    using Base::m_rows;

    using Base::m_size;
};
}

#endif

#undef eigen_assert
#undef Eigen