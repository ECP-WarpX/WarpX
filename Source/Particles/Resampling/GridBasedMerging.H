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
     * based on another GPU vector. The heap-sort functions below where obtained
     * from https://www.geeksforgeeks.org/heap-sort/ and modified for the
     * current purpose.
    */
    struct HeapSort {

        AMREX_GPU_HOST_DEVICE
        void swap(int *a, int *b) const
        {
            const int temp = *b;
            *b = *a;
            *a = temp;
        }

        AMREX_GPU_HOST_DEVICE
        void heapify(int index_array[], int bin_array[], int start, int N, int i) const
        {

            // Initialize largest as current node
            int largest = i;

            // get child node indices for the current node
            const int l = 2 * i + 1;
            const int r = 2 * i + 2;

            // If left child is larger than the current node, set that child to the
            // largest, and the same
            if (l < N && bin_array[index_array[l+start]] > bin_array[index_array[largest+start]]) {
                largest = l;
            }
            if (r < N && bin_array[index_array[r+start]] > bin_array[index_array[largest+start]]) {
                largest = r;
            }
            // If largest is not root
            if (largest != i) {
                swap(&index_array[largest+start], &index_array[i+start]);

                // Recursively heapify the affected sub-tree
                heapify(index_array, bin_array, start, N, largest);
            }
        }

        AMREX_GPU_HOST_DEVICE
        void operator() (int index_array[], int bin_array[], int start, int stop) const
        {
            const int N = stop - start;

            // Build heap (rearrange array)
            for (int i = N / 2 - 1; i >= 0; i--)
            {
                heapify(index_array, bin_array, start, N, i);
            }

            // One by one extract an element
            // from heap
            for (int i = N - 1; i > 0; i--) {

                // Move current root to end
                swap(&index_array[start], &index_array[i+start]);

                // call max heapify on the reduced heap
                heapify(index_array, bin_array, start, i, 0);
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
