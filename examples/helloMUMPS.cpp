/** @file helloMUMPS.cpp

    @brief First example of submodule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

//! [Include namespace]
#include <gismo.h>
#include <gsMUMPS/MUMPSSupport.h>

using namespace gismo;
//! [Include namespace]


int main(int argc, char *argv[])
{
    // Size of global sparse matrix
    index_t mat_size = 10;
    std::string spm(""); // sparse matrix from a file

    gsCmdLine cmd("Testing the use of sparse linear solvers.");
    cmd.addInt("n", "size", "Size of the matrices", mat_size);
    cmd.addString("m", "matrix", "Filename to read sparse matrix and right hand side", spm);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

// // Initialize the MPI environment
//     const gsMpi & mpi = gsMpi::init(argc, argv);

//     // Get the world communicator
//     gsMpiComm comm = mpi.worldComm();

//     //Get size and rank of the processor
//     int _rank = comm.rank();
//     int _size = comm.size();

//     if (0==_rank)
//     {
        gsInfo << "Hello MUMPS!\n";
        #ifdef _OPENMP
        gsInfo<<"Running with "<<omp_get_max_threads()<<" threads.\n";
        #else
        gsInfo<<"Running with 1 thread.\n";
        #endif
//     }

    // Initialize MUMPS solver
    gsEigen::MUMPSLDLT<gsSparseMatrix<>,gsEigen::Lower> solver;
    solver.ICNTL(4) = 0;

    gsSparseMatrix<real_t, RowMajor>  Q;
    gsMatrix<>        b, x;

    if (!spm.empty())
    {
        gsFileData<> fd(spm);
        gsSparseMatrix<> Qcm;
        fd.getFirst(Qcm);
        Q = Qcm;
        fd.getFirst(b);
        mat_size = Q.rows();

    }
    else
    {
        Q.resize(mat_size,mat_size);
        b.resize(mat_size,1);
        x.resize(mat_size,1);

        Q.reserve( gsVector<int>::Constant(mat_size,1) ); // Reserve memory for 1 non-zero entry per column
        for (index_t i = 0; i!=mat_size; ++i)
            Q(i,i) = b(i,0) = i+1;

        Q.makeCompressed(); // always call makeCompressed after sparse matrix has been filled
    }

    solver.compute(Q);
    x = solver.solve(b);

    if (mat_size < 200)
        gsInfo <<"Solution: "<< x.transpose() <<"\n";

    gsInfo<<"Check: "<< ( (b-Q*x).squaredNorm()<1e-8 ) <<"\n";

    return EXIT_SUCCESS;
}
