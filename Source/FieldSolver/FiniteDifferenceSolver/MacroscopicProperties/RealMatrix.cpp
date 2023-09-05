#include <RealMatrix.H>

RealMatrix
RealMatrix::Inverse ()
{
    RealMatrix A = *this;
    RealMatrix Ainv = RealMatrix();
    amrex::Real det;

#if AMREX_SPACEDIM == 2
    det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (det != 0._rt) {
        Ainv[0][0] = A[1][1] / det;
        Ainv[0][1] = -A[0][1] / det;
        Ainv[1][0] = -A[1][0] / det;
        Ainv[1][1] = A[0][0] / det;
    }
#endif
#if AMREX_SPACEDIM == 3
    det  = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
    det += A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
    det += A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

    if (det != amrex::Real(0.0)) {
        Ainv[0][0] = (A[1][1]*A[2][2]-A[1][2]*A[2][1]) / det;
        Ainv[1][0] = (A[1][2]*A[2][0]-A[1][0]*A[2][2]) / det;
        Ainv[2][0] = (A[1][0]*A[2][1]-A[1][1]*A[2][0]) / det;
        Ainv[0][1] = (A[2][1]*A[0][2]-A[2][2]*A[0][1]) / det;
        Ainv[1][1] = (A[2][2]*A[0][0]-A[2][0]*A[0][2]) / det;
        Ainv[2][1] = (A[2][0]*A[0][1]-A[2][1]*A[0][0]) / det;
        Ainv[0][2] = (A[0][1]*A[1][2]-A[0][2]*A[1][1]) / det;
        Ainv[1][2] = (A[0][2]*A[1][0]-A[0][0]*A[1][2]) / det;
        Ainv[2][2] = (A[0][0]*A[1][1]-A[0][1]*A[1][0]) / det;
    }
#endif
    return Ainv;
}

amrex::RealVect
RealMatrix::Dot (amrex::RealVect v)
{
    amrex::RealVect result;
    int i, k;
    for (i = 0; i < AMREX_SPACEDIM; ++i) {
        for (k = 0; k < AMREX_SPACEDIM; ++k) {
            result[i] += (*this)[i][k] * v[k];
        }
    }
    return result;
}

RealMatrix
RealMatrix::Dot (RealMatrix B)
{
    RealMatrix result;
    int i, j, k;
    for (i = 0; i < AMREX_SPACEDIM; ++i) {
        for (j = 0; j < AMREX_SPACEDIM; ++j) {
            for (k = 0; k < AMREX_SPACEDIM; ++k) {
                result[i][j] += (*this)[i][k] * B[k][j];
            }
        }
    }
    return result;
}