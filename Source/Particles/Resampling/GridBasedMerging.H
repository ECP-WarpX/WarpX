/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GRID_BASED_MERGING_H_
#define WARPX_GRID_BASED_MERGING_H_

#include "Resampling.H"

#include "Particles/WarpXParticleContainer_fwd.H"

/**
 * \brief This class implements a particle merging scheme wherein particles
 * are clustered in phase space and particles in the same cluster is merged
 * into two remaining particles. The scheme conserves linear momentum and
 * kinetic energy within each cluster.
 */
class GridBasedMerging: public ResamplingAlgorithm {
public:

    /**
     * \brief Default constructor of the GridBasedMerging class.
     */
    GridBasedMerging () = default;

    /**
     * \brief Constructor of the GridBasedMerging class
     *
     * @param[in] species_name the name of the resampled species
     */
    GridBasedMerging (const std::string& species_name);

    /**
     * \brief A method that performs merging for the considered species.
     *
     * @param[in] pti WarpX particle iterator of the particles to resample.
     * @param[in] lev the index of the refinement level.
     * @param[in] pc a pointer to the particle container.
     */
    void operator() (WarpXParIter& pti, int lev, WarpXParticleContainer* pc) const final;

    /**
     * \brief This merging routine requires functionality to sort a GPU vector
     * based on another GPU vector's values. The heap-sort functions below were
     * obtained from https://www.geeksforgeeks.org/iterative-heap-sort/ and
     * modified for the current purpose. It achieves the same as
     * ```
     * std::sort(
     *    sorted_indices_data + cell_start, sorted_indices_data + cell_stop,
     *    [&momentum_bin_number_data](size_t i1, size_t i2) {
     *        return momentum_bin_number_data[i1] < momentum_bin_number_data[i2];
     *    }
     * );
     * ```
     * but with support for device execution.
    */
    struct HeapSort {

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void swap(int &a, int &b) const
        {
            const int temp = b;
            b = a;
            a = temp;
        }

        AMREX_GPU_HOST_DEVICE
        void operator() (int index_array[], const int bin_array[], const int start, const int n) const
        {
            // sort index_array into a max heap structure
            for (int i = 1; i < n; i++)
            {
                // move child through heap if it is bigger than its parent
                if (bin_array[index_array[i+start]] > bin_array[index_array[(i - 1)/2 + start]])
                {
                    int j = i;
                    // swap child and parent until branch is properly ordered
                    do {
                        swap(index_array[j+start], index_array[(j - 1)/2 + start]);
                        j = (j - 1) / 2;
                    } while (bin_array[index_array[j+start]] > bin_array[index_array[(j - 1)/2 + start]]);
                }
            }

            for (int i = n - 1; i > 0; i--)
            {
                // swap value of first (now the largest value) to the new end point
                swap(index_array[start], index_array[i+start]);

                // remake the max heap
                int j = 0, index;
                do {
                    index = (2 * j + 1);

                    // if left child is smaller than right child, point index variable to right child
                    if (index < (i - 1) && bin_array[index_array[index+start]] < bin_array[index_array[index+1+start]]) {
                        index++;
                    }
                    // if parent is smaller than child, swap parent with child having higher value
                    if (index < i && bin_array[index_array[j+start]] < bin_array[index_array[index+start]]) {
                        swap(index_array[j+start], index_array[index+start]);
                    }
                    j = index;
                } while (index < i);
            }
        }
    };

private:
    int m_min_ppc = 100;
    int m_ntheta = 120;
    int m_nphi = 60;
    amrex::ParticleReal m_delta_ur = 1e5;
};
#endif // WARPX_GRID_BASED_MERGING_H_
