if(WarpX_PSATD)
    # cuFFT  (CUDA)
    #   TODO: check if `find_package` search works

    # rocFFT (HIP)
    if(WarpX_COMPUTE STREQUAL HIP)
        find_package(rocfft REQUIRED)

    # FFTW   (NOACC, OMP, SYCL)
    elseif(NOT WarpX_COMPUTE STREQUAL CUDA)
        # On Windows, try searching for FFTW3(f)Config.cmake files first
        #   Installed .pc files wrongly and unconditionally add -lm
        #   https://github.com/FFTW/fftw3/issues/236

        # On Linux & macOS, note Autotools install bug:
        #   https://github.com/FFTW/fftw3/issues/235
        # Thus, rely on .pc files

        set(WarpX_FFTW_SEARCH_VALUES PKGCONFIG CMAKE)
        set(WarpX_FFTW_SEARCH_DEFAULT PKGCONFIG)
        if(WIN32)
            set(WarpX_FFTW_SEARCH_DEFAULT CMAKE)
        endif()
        set(WarpX_FFTW_SEARCH ${WarpX_FFTW_SEARCH_DEFAULT}
            CACHE STRING "FFTW search method (PKGCONFIG/CMAKE)")
        set_property(CACHE WarpX_FFTW_SEARCH PROPERTY STRINGS ${WarpX_FFTW_SEARCH_VALUES})
        if(NOT WarpX_FFTW_SEARCH IN_LIST WarpX_FFTW_SEARCH_VALUES)
            message(FATAL_ERROR "WarpX_FFTW_SEARCH (${WarpX_FFTW_SEARCH}) must be one of ${WarpX_FFTW_SEARCH_VALUES}")
        endif()
        mark_as_advanced(WarpX_FFTW_SEARCH)

        if(WarpX_FFTW_SEARCH STREQUAL CMAKE)
            if(WarpX_PRECISION STREQUAL "DOUBLE")
                find_package(FFTW3 CONFIG REQUIRED)
            else()
                find_package(FFTW3f CONFIG REQUIRED)
            endif()
        else()
            if(WarpX_PRECISION STREQUAL "DOUBLE")
                find_package(PkgConfig REQUIRED QUIET)
                pkg_check_modules(fftw3 REQUIRED IMPORTED_TARGET fftw3)
            else()
                find_package(PkgConfig REQUIRED QUIET)
                pkg_check_modules(fftw3f REQUIRED IMPORTED_TARGET fftw3f)
            endif()
        endif()
    endif()

    # create an IMPORTED target: WarpX::thirdparty::FFT
    if(WarpX_COMPUTE STREQUAL CUDA)
        # CUDA_ADD_CUFFT_TO_TARGET(WarpX::thirdparty::FFT)
        make_third_party_includes_system(cufft FFT)
    elseif(WarpX_COMPUTE STREQUAL HIP)
        make_third_party_includes_system(roc::rocfft FFT)
    else()
        if(WarpX_PRECISION STREQUAL "DOUBLE")
            if(FFTW3_FOUND)
                # subtargets: fftw3, fftw3_threads, fftw3_omp
                if(WarpX_COMPUTE STREQUAL OMP AND TARGET FFTW3::fftw3_omp)
                    make_third_party_includes_system(FFTW3::fftw3_omp FFT)
                else()
                    make_third_party_includes_system(FFTW3::fftw3 FFT)
                endif()
            else()
                make_third_party_includes_system(PkgConfig::fftw3 FFT)
            endif()
        else()
            if(FFTW3f_FOUND)
                # subtargets: fftw3f, fftw3f_threads, fftw3f_omp
                if(WarpX_COMPUTE STREQUAL OMP AND TARGET FFTW3::fftw3f_omp)
                    make_third_party_includes_system(FFTW3::fftw3f_omp FFT)
                else()
                    make_third_party_includes_system(FFTW3::fftw3f FFT)
                endif()
            else()
                make_third_party_includes_system(PkgConfig::fftw3f FFT)
            endif()
        endif()
    endif()
endif()
