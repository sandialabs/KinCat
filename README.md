/*====================================================================================== 
kincat version 1.0 Copyright (2023) NTESS https://github.com/sandialabs/kincat Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software. This file is part of KinCat. KinCat is open-source software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. You should have received a copy of the GNU Lesser General Public License along with KinCat. If not, see http://www.gnu.org/licenses/.

Questions? Contact Craig Daniels at cjdanie@sandia.gov

Sandia National Laboratories, Albuquerque, NM, USA ======================================================================================*/

KinCat is an open-source 2D lattice KMC simulator. It is written in Python and C++ with the Kokkos library.

KinCat build instructions, file explanations, and other documentation can be found in the kincat.pdf file. It can also be generated in other formats using sphinx from the source files found in the src/sphinx directory. Example input and post-processing files can be found in the src/examples directory. The core c++ code is stored in src/core. However, the preprocessing code, KinCatPy is found in the kincatpy directory. Note that src/unit-test holds code for unit tests to confirm functionality, but these are deprecated. 
