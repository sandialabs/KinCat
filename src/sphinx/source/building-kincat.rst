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

Weaver is a Sandia National Laboratories compute resource equipped with IBM Power9 processors accelerated by NVIDIA V100 GPU. We use the following modules. These modules may be changed as the system is upgraded.

.. code-block:: bash

   module purge
   module load devpack/20210226/openmpi/4.0.5/gcc/7.2.0/cuda/10.2.2
   module swap cmake cmake/3.19.3
   module swap boost boost/1.75.0

   which gcc
   /home/projects/ppc64le/gcc/7.2.0/bin/gcc
   
   echo ${BOOST_ROOT}
   /home/projects/ppc64le-pwr9/spack/opt/spack/linux-rhel7-power9le/gcc-7.2.0/boost-1.75.0-qf3d47g2bt3dhlbruldwpqfu3rqkrdtk

   which cmake
   /home/projects/ppc64le/cmake/3.19.3/bin/cmake

   
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

Build KinCat and link with Kokkos and Gtest. The following script shows how to compile KinCat on OSX while linking with TPLs explained above. 

.. code-block:: bash
		
   cd ${KINCAT_BUILD_PATH}
   cmake \
     -D CMAKE_INSTALL_PREFIX=${KINCAT_INSTALL_PATH} \
     -D CMAKE_CXX_COMPILER="${MYCXX}" \
     -D CMAKE_CXX_FLAGS="-g" \
     -D CMAKE_EXE_LINKER_FLAGS="" \
     -D CMAKE_BUILD_TYPE=RELEASE \
     -D KINCAT_SITE_TYPE="char" \
     -D KINCAT_ENABLE_DEBUG=OFF \
     -D KINCAT_ENABLE_VERBOSE=ON \
     -D KINCAT_ENABLE_TEST=ON \
     -D KINCAT_ENABLE_EXAMPLE=ON \
     -D KOKKOS_INSTALL_PATH="${KOKKOS_INSTALL_PATH}" \
     -D GTEST_INSTALL_PATH="${GTEST_INSTALL_PATH}" \
     ${KINCAT_REPOSITORY_PATH}/code/src
   make -j install
   export KINCAT_INSTALL_PATH=${KINCAT_INSTALL_PATH}

To install KinCat on Weaver (GPU platform), replace the C++ compiler with nvcc wrapper, providing ``-D CMAKE_COMPILER="${KOKKOS_INSTALL_PATH}/bin/nvcc_wrapper`` instead.    
A successful installation creates the following directory structure in ``${KINCAT_INSTALL_PATH}``. Note that the ``site_type`` is determined at compile time. In the above cmake configuration, the type is set ``char`` and the max number of species in KinCat is 256. For a bigger simulation, users can set this ``short`` or ``int``.  

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

.. autosummary::
   :toctree: generated


