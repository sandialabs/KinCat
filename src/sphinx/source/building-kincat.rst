Building KinCat
===============

For convenience, we explain how to build KinCat using the following environment variables that users can modify according to their working environments. Note that cmake requires use of a separate build directory protecting the source code. 

.. code-block:: bash

   # repositories
   export KINCAT_REPOSITORY_PATH=/where/you/clone/kincat/git/repo   
   export KOKKOS_REPOSITORY_PATH=/where/you/clone/kokkos/git/repo
   export GTEST_REPOSITORY_PATH=/where/you/clone/gtest/git/repo

   # build directories
   export KINCAT_BUILD_PATH=/where/you/build/kincat   
   export KOKKOS_BUILD_PATH=/where/you/build/kokkos
   export GTEST_BUILD_PATH=/where/you/build/gtest

   # install directories
   export KINCAT_INSTALL_PATH=/where/you/install/kincat   
   export KOKKOS_INSTALL_PATH=/where/you/install/kokkos
   export GTEST_INSTALL_PATH=/where/you/install/gtest

System Requirements 
--------------------

- CMake version 3.16 and higher.
- Compiler supporting C++14 and higher standards.
- Boost library 1.75 and higher.
- Jupyter notebook (or lab) and Python3 for visualization and post-processing.
- sphinx for documentation.
- OpenMPI for running KinCat with multiple data.

Mac OSX
^^^^^^^

We use macports for the standard software distributions of required tools and libraries. Install or check the following tools are installed using macports. 

.. code-block:: bash

   sudo port install cmake clang-13 boost python310 py38-sphinx py38-sphinx_rtd_theme openmpi-clang13

   which mpirun-openmpi-clang13 
   mpirun-openmpi-clang13 is /opt/local/bin/mpirun-openmpi-clang13
   
   which clang++-mp-13 
   clang++-mp-13 is /opt/local/bin/clang++-mp-13

   which cmake
   cmake is /opt/local/bin/cmake

   export MYCXX=clang++mp-13

   pip install jupyterlab

Weaver
^^^^^^

Weaver was a Sandia National Laboratories compute resource equipped with IBM Power9 processors accelerated by NVIDIA V100 GPU. Although it was recently taken offline, we include notes relating to Weaver, hoping that they may be a useful guide for building and running on other resources. We used the following modules.

.. code-block:: bash

   module purge
   module load gcc/11.3.0
   module load cuda/11.8.0
   module load openmpi
   module load cmake
   module load openblas
   module load ninja
   module load git
   module load boost
   module load hdf5

   export MYCXX=${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper
   export INCLUDE=${OPENMPI_ROOT}/include:${INCLUDE}

   which gcc
   /projects/ppc64le-pwr9-rhel8/compilers/gcc/11.3.0/gcc/8.3.1/base/tchbki3/bin/gcc
   
   echo ${BOOST_ROOT}
   /projects/ppc64le-pwr9-rhel8/tpls/boost/1.80.0/gcc/11.3.0/base/wrmyc3o

   which cmake
   /projects/ppc64le-pwr9-rhel8/utilities/cmake/3.29.6/gcc/8.5.0/base/7ne4ua7/bin/cmake

   
How to Build
------------

First, clone Kokkos, GTEST and KinCat code repositories. 

.. code-block:: bash

   git clone https://github.com/kokkos/kokkos.git ${KOKKOS_REPOSITORY_PATH}
   git clone https://github.com/google/googletest.git ${GTEST_REPOSITORY_PATH}
   git clone https://github.com/sandialabs/KinCat.git ${KINCAT_REPOSITORY_PATH}

Kokkos
^^^^^^

This builds Kokkos on Intel Haswell architectures and installs Kokkos to ``${KOKKOS_INSTALL_PATH}``. For more details, see [Kokkos github pages](https://github.com/kokkos/kokkos). We can use this script for OSX.

.. code-block:: bash
		
   cd ${KOKKOS_BUILD_PATH}
   cmake \
     -D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_PATH}" \
     -D CMAKE_CXX_COMPILER="${MYCXX}"  \
     -D Kokkos_ENABLE_SERIAL=ON \
     -D Kokkos_ENABLE_OPENMP=ON \
     -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
     -D Kokkos_ARCH_HSW=ON \
     ${KOKKOS_REPOSITORY_PATH}
   make -j install

On Weaver, we compile Kokkos for NVIDIA GPUs. Note that we use Kokkos nvcc_wrapper as its compiler instead of directly using the nvcc compiler. The architecture flag indicates that the host architecture is IBM Power9 and the GPU architecture is Volta70 generation.

.. code-block:: bash
		
   cd ${KOKKOS_BUILD_PATH}
   cmake \
     -D CMAKE_INSTALL_PREFIX="${KOKKOS_INSTALL_PATH}" \
     -D CMAKE_CXX_COMPILER="${KOKKOS_REPOSITORY_PATH}/bin/nvcc_wrapper"  \
     -D Kokkos_ENABLE_SERIAL=ON \
     -D Kokkos_ENABLE_OPENMP=ON \
     -D Kokkos_ENABLE_CUDA:BOOL=ON \
     -D Kokkos_ENABLE_CUDA_UVM:BOOL=OFF \
     -D Kokkos_ENABLE_CUDA_LAMBDA:BOOL=ON \
     -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
     -D Kokkos_ARCH_VOLTA70=ON \
     -D Kokkos_ARCH_POWER9=ON \
     ${KOKKOS_REPOSITORY_PATH}
   make -j install

GTEST
^^^^^

We use GTEST as our testing infrastructure. With the following cmake script, the GTEST can be compiled and installed.

.. code-block:: bash
		
   cd ${GTEST_BUILD_PATH}
   cmake \
     -D CMAKE_INSTALL_PREFIX="${GTEST_INSTALL_PATH}" \
     -D CMAKE_CXX_COMPILER="${MYCXX}"  \
     ${GTEST_REPOSITORY_PATH}
   make -j install

Boost
^^^^^

The Boost library may be installed by the following script.

.. code-block:: bash

    export BOOST_INSTALL_PATH=/where/you/install/boost
    mkdir -p ${BOOST_INSTALL_PATH}
    cd ${BOOST_INSTALL_PATH}
    curl -L https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.bz2 -o boost_1_75_0.tar.bz2
    tar -xvf boost_1_75_0.tar.bz2
    export BOOST_ROOT=${BOOST_INSTALL_PATH}/boost_1_75_0

KinCat
^^^^^^

Build KinCat and link with Kokkos and Gtest. The following script shows how to compile KinCat on Weaver while linking with TPLs explained above. 

.. code-block:: bash
		
   cd ${KINCAT_BUILD_PATH}
   cmake \
        -D CMAKE_INSTALL_PREFIX=${KINCAT_INSTALL_PATH}\
        -D CMAKE_CXX_COMPILER="${MYCXX}" \
        -D CMAKE_CXX_FLAGS="-g -I${OPENMPI_ROOT}/include " \
        -D CMAKE_EXE_LINKER_FLAGS=""\
        -D CMAKE_BUILD_TYPE=RELEASE \
        -D KINCAT_SITE_TYPE="char" \
        -D KINCAT_ENABLE_DEBUG=OFF \
        -D KINCAT_ENABLE_VERBOSE=ON \
        -D KINCAT_ENABLE_TEST=ON \
        -D KINCAT_ENABLE_EXAMPLE=ON \
        -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
        -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
        -D HDF5_INCLUDE_DIRS="${HDF5_INC}" \
        -D HDF5_LIBRARY_DIRS="${HDF5_LIB}" \
        -D HDF5_LIBRARIES="hdf5" \
        ${KINCAT_REPOSITORY_PATH}/src

To install KinCat on OSX, we use the following script. Note that in both scripts the HDF5 related keys are only needed if HDF5 functionality is desired. 

.. code-block:: bash
    
   cmake \
    -D CMAKE_INSTALL_PREFIX=${KINCAT_INSTALL_PATH} \
    -D CMAKE_C_COMPILER="clang-mp-12" \
    -D CMAKE_CXX_COMPILER="clang++-mp-12" \
    -D CMAKE_CXX_FLAGS="-g" \
    -D CMAKE_EXE_LINKER_FLAGS="" \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D KINCAT_ENABLE_DEBUG=OFF \
    -D KINCAT_ENABLE_VERBOSE=ON \
    -D KINCAT_ENABLE_TEST=ON \
    -D KINCAT_ENABLE_EXAMPLE=ON \
    -D KOKKOS_INSTALL_PATH="${HOME}/kokkos/kokkos_install/release" \
    -D GTEST_INSTALL_PATH="/opt/local" \
    -D HDF5_INCLUDE_DIRS="/usr/local/hdf5/include" \
    -D HDF5_LIBRARY_DIRS="/usr/local/hdf5/lib" \
    -D HDF5_LIBRARIES="hdf5"   \
    ${KINCAT_SRC_PATH}


A successful installation creates the following directory structure in ``${KINCAT_INSTALL_PATH}``. Note that the ``site_type`` is determined at compile time. In the above cmake configuration, the type is set ``char`` and the max number of species in KinCat is 256. For a bigger simulation, users can set this ``short`` or ``int``, though this is unlikely to be needed. Also note that currently the ``char`` setting is overwritten to use ``short`` instead. This is due to requiring a signed type since KinCatPy uses -1 to indicate an unspecified site in a configuration. 


.. code-block:: bash

   - bin
     - kincat.x: an executable for solving single problem
     - kincat-batch.x: an executable for solving multiple problem with batch parallelism 
     - plot-dump.ipynb: visualization for jupyter notebook
   - examples
     - RuO2-dictionary.json : auto-generated from KinCatPy
     - RuO2-rates.json : used to set process rates for simulation
     - kincatpy
       - readme.txt
       - ruo2_input.txt
     - non-batch
       - readme.txt
       - ex1-input.json
       - ex2-input.json
       - plot-ex1.py
       - plot-ex2.py
     - batch
       - readme.txt
       - ex3-input.json
       - plot-ex3.py
       - input-rates-override-RuO2.json
   - include
     - kincat
       - header files
   - lib (or lib64)	 
     - cmake: cmake environment when other software interface KinCat via cmake
     - libkincat.a
   - unit-test
     - kincat-test.x: unit test executable
     - test-files: sample files that will be used in test 

Spack
^^^^^
As an alternative to the 'manual' build instructions above, we have also created Spack build instructions. Spack is a package manager intended for scientific software, allowing for software and all of it's dependencies to be built with minimal user interaction. Note that this functionality has not been tested extensively. If Spack is not already available, it needs to be set up:

.. code-block:: bash

    git clone --depth=100 --branch=releases/v0.21 https://github.com/spack/spack.git ~/spack
    cd ~/spack/
    . share/spack/setup-env.sh

While Spack often has recognized software packages, for now KinCat can only be built with a local spack package file.

.. code-block:: bash

    spack repo create kincat-local-repo
    spack repo add {PATH_TO_SPACK_LOCAL}/kincat-local-repo
    mkdir kincat-local-repo/packages/kincat
    cp kincat/spack/kincat-package.py kincat-local-repo/packages/kincat/package.py

This build file needs to be modified depending on the version of KinCat desired. Both the url reference and the version number and sha256 key should be updated. The version of Kokkos to be used may also be modified. Once Spack and the local repo are set up, KinCat should be able to be installed with just a single command:

.. code-block:: bash
    spack install kincat

Because all the dependencies are also being installed, this may take several minutes. Note that the executable may be saved in an obscure directory.

.. autosummary::
   :toctree: generated