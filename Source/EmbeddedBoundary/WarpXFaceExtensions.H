//
// Created by Lorenzo Giacomel on 26/08/2021.
//

#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEEXTENSIONS_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEEXTENSIONS_H
/**
* \brief For the face of cell pointing in direction idim, return the number of faces
* we need to intrude with the one-way extension. Returns only one or zero: one if the
* face can be extended with the the one-way extension, zeros if it can't.
*
* \param[in] cell \c Dim3 storing the indices of the face to extended
* \param[in] S_ext amount of area needed for the extension
* \param[in] S_red \c Array4 storing the amount of  area each face can still give away
* \param[in] flag_info_face \c Array4 storing face information
* \param[in] flag_ext_face \c Array4 storing face information
* \param[in] idim normal direction to the face in consideration (0 for x, 1 for y, 2 for z)
*/
AMREX_GPU_DEVICE
int ComputeNBorrowOneFaceExtension(const amrex::Dim3 cell, const amrex::Real S_ext,
                                          const amrex::Array4<amrex::Real> &S_red,
                                          const amrex::Array4<int> &flag_info_face,
                                          const amrex::Array4<int> &flag_ext_face, int idim);
/**
* \brief For the face of cell pointing in direction idim, return the number of faces
* we need to intrude with the eight-ways extension.
*
* \param[in] cell \c Dim3 storing the indices of the face to extended
* \param[in] S_ext amount of area needed for the extension
* \param[in] S_red \c Array4 storing the amount of  area each face can still give away
* \param[in] S_red \c Array4 storing the area of face
* \param[in] flag_info_face \c Array4 storing face information
* \param[in] idim normal direction to the face in consideration (0 for x, 1 for y, 2 for z)
*/
AMREX_GPU_DEVICE
int ComputeNBorrowEightFacesExtension(const amrex::Dim3 cell, const amrex::Real S_ext,
                                             const amrex::Array4<amrex::Real> &S_red,
                                             const amrex::Array4<amrex::Real> &S,
                                             const amrex::Array4<int> &flag_info_face,
                                             int idim);
#endif //WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEEXTENSIONS_H
