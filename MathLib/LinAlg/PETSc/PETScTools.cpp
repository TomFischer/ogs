/*!
   \file
   \brief Definition of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScTools.h"
#include "PETScVector.h"

namespace MathLib
{
void applyKnownSolution(PETScMatrix& A, PETScVector& b, PETScVector& x,
                        const std::vector<PetscInt>& vec_knownX_id,
                        const std::vector<PetscScalar>& vec_knownX_x)
{
    A.finalizeAssembly();
    x.finalizeAssembly();
    if (vec_knownX_id.size() > 0)
    {
        x.set(vec_knownX_id, vec_knownX_x);
    }
    x.finalizeAssembly();

    const PetscScalar one = 1.0;
    const PetscInt nrows = static_cast<PetscInt>(vec_knownX_id.size());

    PETSc_Mat & mat(A.getRawMatrix());

    // Each process will only zero its own rows.
    // This avoids all reductions in the zero row routines
    // and thus improves performance for very large process counts.
    // See PETSc doc about MAT_NO_OFF_PROC_ZERO_ROWS.
    // MatSetOption(mat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);

    // Keep the non-zero pattern for the assignment operator.
    MatSetOption(mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    Vec& solution(x.getRawVector());
    Vec& rhs(b.getRawVector());
    if (nrows > 0)
    {
        MatZeroRowsColumns(mat, nrows, vec_knownX_id.data(), one, solution, rhs);
    }
    else
    {
        MatZeroRowsColumns(mat, 0, PETSC_NULL, one, solution, rhs);
    }

    A.finalizeAssembly();
    x.finalizeAssembly();
    b.finalizeAssembly();
}

}  // end of namespace MathLib
