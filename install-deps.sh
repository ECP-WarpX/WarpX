brew update
brew tap openpmd/openpmd
brew install adios2      # for openPMD
brew install ccache
brew install cmake
brew install fftw        # for PSATD
brew install git
brew install hdf5-mpi    # for openPMD
brew install libomp
brew unlink gcc
brew link --force libomp
brew install pkg-config  # for fftw
brew install open-mpi
brew install openblas    # for PSATD in RZ
brew install openpmd-api # for openPMD
