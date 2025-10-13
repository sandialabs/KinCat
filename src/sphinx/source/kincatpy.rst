KinCatPy
============

**KinCatPy** is an open source Python script to write dictionary input files for use with KinCat. It receives inputs describing the lattice, possible species and process mechanisms, and the range of possible interactions. It outputs a KinCat dictionary file which contains information about all possible symmetrically unique process instances and all configurations needed to define those (given the input parameters). A given process mechanism, such as desorption from the surface, may have different rates due to lateral interactions with species on nearby sites. KinCatPy defines a unique process instance for each process that may occur between two unique initial and final configurations. Thus, rates for each instance can be adjusted to account for the lateral interactions, giving the user as much specificity in setting rates as desired. KinCatPy is adapted from KineCluE, another open source code intended to calculate transport coefficients of mobile clusters. Note that while KineCluE was written to handle bulk crystals, KinCatPy is only intended for 2D lattice descriptions. 

KinCatPy is run in the command line window or terminal as ``python ./path_to_kincatpy/kincatpy.py input.txt``.  

Main Input
----------

Some of the input language and structure remains the same as KinCluE (https://github.com/lukamessina/kineclue). While we outline the inputs here, the user may find KineCluE documentation to be helpful background. It is also helpful to understand how the Bravais lattice construction allows any lattice to be defined by a repeating unit cell and basis vectors within that cell. The following input is extracted from the RuO2 example.  

.. code-block:: javascript

   & CRYSTAL test # creation of the crystal with 2 vectors of 2 components
   6.43 0.0
   0.0 3.12
   & BASIS s 2
   0 0 0 
   1 0.5 0 
   & UNIQUEPOS 2 # describes the sites in the crystal that the atoms and defect will occupy.
   s 0 0 
   s 0.5 0 
   & RANGE 3.0  # the float is the interaction radius
   & SPECIES 2
   2 1 1 O 0.1
   2 1 1 CO 0.1
   & PROCMECH
   #adsorption processes
   %% 1 CO_ads_cus
   s 0.5 0 0 > 2 
   %% 1 CO_ads_br
   s 0 0 0 > 2
   ...
   #desorption processes
   %% 1 CO_des_cus
   s 0.5 0 2 > 0
   %% 1 CO_des_br
   s 0 0 2 > 0
   ...
   ##diffusion processes
   %% 2 CO_cus_cus
   s 0.5 0 2 > 0 
   s 0.5 1 0 > 2 
   %% 2 CO_br_br
   s 0 0 2 > 0
   s 0 1 0 > 2
   ...
   ##Recombination/desorption of CO2
   %% 2 CO_cus_O_cus
   s 0.5 1 1 > 0
   s 0.5 0 2 > 0
   %% 2 CO_br_O_br
   s 0 1 1 > 0
   s 0 0 2 > 0
   ...

It is useful to recognize two possible coordinate systems that KinCatPy will accept. The first is the standard x-y or Cartesian coordinates. The second is in terms of lattice edge vectors which define the unit cell that tiles the plane. These are referred to as orthogonal and supercell coordinate systems and 'o' and 's' respectively are used to specify which is used in the input file. Also recognize that ``&`` is used to indicate a new input command and ``#`` may be used to comment the remainder of the line. 

* ``& CRYSTAL`` is used to define the unit cell of the lattice. The string after the command can be used to name the lattice if the user wishes. The next lines comprise two edge vectors that define the lattice unit cell. These values are always given in the orthogonal coordinate system, and these values are used to map between the orthogonal and supercell coordinate systems moving forward.

* ``& BASIS`` command is used to define the basis vectors of the lattice. If it is not present, then a single basis vector of (0,0) is assumed. The coordinate system of the basis vectors ('o' or 's') must be specified, and the number of basis vectors given. For each basis vector, a line of the form ``site_type a b`` is required. The ``site_type`` is an integer, and can be used to indicate if the site is fundamentally the same or different than another which is important for symmetry considerations. The ``a b`` terms are the coordinates of the basis vector in the specified coordinate system. 

* ``& UNIQUEPOS`` command is similar to the ``& BASIS`` command, and must always be included. The ``BASIS`` vectors are used to define the symmetry of the system, while ``UNIQUEPOS`` vectors define sites the species can occupy. For ``UNIQUEPOS``, only one vector for each symmetrically equivalent site is required. Thus, there should be the same number of vectors here as there are unique site types given in the ``BASIS`` section. In ``UNIQUEPOS``, the number of vectors is again required, but the coordinate system is defined for each vector, in the form ``s(o) a b``. 

* ``& RANGE`` command is used to define the interaction range, and the distance is in the units of the orthogonal coordinate system. Each configuration will include all sites within this range of a process, so increasing it may dramatically increase the complexity of the KMC model generated. 

* ``& SPECIES`` command is used to give information about the possible species in the system. Currently, each species is considered as a single point, so it is impossible to include information about orientation or multi-site binding in the model. The number of species needs to be given, and then each species is specified by a line of the form ``n pos_bools... name size``. The number of each species is intended to be provided by ``n`` (as used in KineCluE). However, the user just needs to ensure that this integer is larger than the number of that species specified in any single ``JUMPMECH``. Next, ``pos_bools`` refers to a set of 0/1 flags indicating if that species can occupy the ``UNIQUEPOS`` given. There should be as many flags are there are ``UNIQUEPOS`` specified. The ``name`` is a unique text string. The last number is ``size`` indicating the implied radius of the species (orthogonal units). Species with sufficient size may 'block' neighboring sites and prevent other species from occupying them. The size may be set arbitrarily small if this functionality is not desired.

* ``& PROCMECH`` command is used to define possible process mechanisms through sets of constraints. After the ``& PROCMECH`` command is given, a new process mechanism definition is denoted by a line of form ``%% n_constraints name``. Each process name should be unique. ``n_constraints`` is the number of constraints that define the process. Each contstraint on the process is given by a line of the form ``s(o) a b ini_species_type > fin_species_type``. The letter ``s(o)`` is used to denote which coordinate system will be used for ``a b``, which are the site coordinates. ``ini_species_type`` and ``fin_species_type`` are the integer indices of the initial and final species respecively. The species inputed using the ``& SPECIES`` command are assigned indicies in their list order, beginning with one. The species index 0 is reserved to indicate a vacant site. Note that each constraint definition consists of a site and the change of species on that site. Thus, defining a diffusion jump involves two constraints, one for the species disapearing from the initial site and another for the species appearing in the final site. So called 'bystander' or 'spectator' species that do not change may also be included in the constraints. For example, a reaction between two 'A' species may be catalyzed by the presence of a 'B' species. The B species would need to be included in the constraints with an approriate coordinate relative to the A species, but the initial and final ``species_type`` of that constraint would be the same. Note that a only single symmetry of a process needs to be included. Symmetric processes are found by accounting for the symmetries of the lattice. 

``& DIRECTORY`` command can be used to specify an output directory where the output files will be generated. If this command is not included, an output directory named 'CALC/' is auto-generated in the working directory.

``& FULLSYM`` command can be used to create a dicitonary where the configurations and process instance dictionaries are not reduced by symmetry operations.

``& UNICONFIG`` command can be used to create a dictionary where the configurations and process instance dictionaries are reduced by symmetry to a single configuration template. Symmetric processes are recovered by symmetry operations.

If neither the ``FULLSYM`` or ``UNICONFIG`` flags are included, KinCatPy will use the default 'Sets' dictionary style. In this style, processes are grouped by involed sites to create multiple configuration templates. Symmetric processes are again recovered by symmetry, and this dictionary style will produce the smallest catalogues.

Note again that the complexity of the model given will greatly depend on 1) the number of species included, 2) the range of interactions accounted for (set by ``RANGE``), and 3) the choice of process definitions. The computational demands of KinCatPy are largely determined by the number of unique configurations that need to be specified. Since it will identify all possible permutations of the species arranged on a set of lattice sites, the number of species and number of sites in the configuration definition will both have an exponential impact on the number of configurations overall. The configuration includes all lattice sites within the interaction range of a site that is changed by a process (the species on that site changes, not just that it is included in a bystander constraint). Increasing the interaction range will increase the number of sites needed for the configuration definition. However, for the reduced symmetry cases, careless process definition may also increase the number of sites defined. The process definitions should use the symmetries that overlap as much as possible. In this example defining CO adsorption to 'cus' site process using site ``0 0.5`` and the CO desorption from 'cus' site process using site ``0 -0.5`` is valid, since the sites are symmetrically equivalent to each other. However, when compared to defining both processes with the same site, this would lead to an expansion of the configuration definition and increase the complexity of the model defined by KinCatPy. Since all symmetries of each process are included in the full symmetry case, it does not matter. We reccomend using the sets case for most simulations, but the full-symmetry case may be easier for the user to interpret.

We also note that for large model dictionaries, a multi-thread script (kincat_event_calc.py) is called. However, this only significantly improves efficiency with large numbers of configurations, and so is not called if the model has fewer than 100,000 configurations. Even models with over a 100,000 configurations and over 1,000,000 process instances only takes a few minutes to generate. We also note that KinCatPy assumes that the kincat_event_calc.py file will be stored in the directory above where KinCatPy is running.


.. autosummary::
   :toctree: generated