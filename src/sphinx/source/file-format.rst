File Format
============

KinCat requires three input files specifying 1) solver, 2) dictionary, and 3) rates. The solver input can be specified by a user to control the KMC simulations while the dictionary input is auto-generated from a pre-processor, KinCatPy. The rates input is separated as we may want to use/test different set of rates that can be computed by multiple sources e.g., RMG.  

Main Input
----------

The following json input is the main input file extracted from the RuO2 example 1.

.. code-block:: javascript

   {
       "kincat": {
           "dictionary" : "${KINCAT_INSTALL_PATH}/src/examples/RuO2-dictionary.json",
           "rates" : "${KINCAT_INSTALL_PATH}/src/examples/RuO2-rates.json",
           "sites" : {
               "lattice" : [ 10, 10, 2 ],
               "type" : "random",
               "random seed" : 12345,
               "random fill ratio" : [0.0, 0.0],
               "filename" : "dump-sites.json"
           },
           "dump" : {
                "dump filename" : "dump-ex1.json",
                "dump interval" : 2.5E-5
             },
           "statistics" : {
                "stats filename" : "stats-ex1.json",
                "types" : ["species_coverage","process_counts"],
                "stats interval" : 2E-5
            },
           "solver" : {
               "type" : "serial-rta",
               "random seed" : 13542,
               "random pool size" : 10000,
               "max number of kmc kernel launches" : 1000,
               "max number of kmc steps per kernel launch" : 50000,
               "time range" : [ 0, 1E-3, 5E-6 ]
           }
       }
   }

* ``kincat`` object includes the following items:

    * ``dictionary`` specifies the file path to the dictionary json file.
    * ``rates`` specifies the file path to the rate json file.
    * ``sites`` object includes lattice information and its initialization method.
    
        * ``lattice``: 2D lattice object. The array input implies ``[ x-length, y-length, n-basis ]``. The lengths are specified in unit cells, and must be integers.
        * ``type``: lattice initialization method. 

            * ``random`` value fills the lattice randomly with percentages specified by ``random fill ratio`` keyword; otherwise, ``random fill ratio`` is ignored. ``random fill ratio`` should be an array with as many entries as the number of species (excluding vacant sites), and in the same order as they were defined in KinCatPy. The sum of these fill ratios should not exceed one. A ``random seed`` is given and used just for the initialization routine.
            * ``file`` value loads site specification from a file given by the restart file ``filename``; otherwise, ``filename`` is ignored. Note that while the initial configuration may be read-in from a file, the initial simulation time is always set by ``time begin`` as set below. 
	
    * ``dump`` object relates to an output file with the full simulation state.

        * ``dump filename`` object specifies the name of the dump output file.  
        * ``dump interval`` object indicates the minimum time between outputs of the full simulation state.  
    * ``statistics`` object relates to an output file of various desired statistics for the KMC simulation.

        * ``filename`` object specifies the name of the statistics output file. 
        * ``types`` object is a list of the various statistics desired by the user. Currently implemented types include ``"process_counts"``, ``"species_coverage"``, and ``"site_species_coverage"``. The ``"process_counts"`` option outputs the number of times each process has occurred up to that point in the KMC simulation which can be used to calculate turn-over frequencies. The ``"species_coverage"`` option outputs the fraction of the lattice sites filled by each species, and ``"site_species_coverage"`` option further breaks down this fraction by each unique type of lattice site.  Further types may be added later, or added by the user.
        * ``stats interval`` object indicates the minimum time between outputs of the statsitics listed in ``types``. 
    * ``solver`` object include the following items to specify the KMC simulation.

        * ``type`` specifies the algorithm selected: ``serial-rta``, ``sublattice``, ``batch-rta``.

            * ``serial-rta`` solver implements the BKL, n-fold, or 1st Order KMC algorithm.
            * ``batch-rta`` solver implements the serial-rta algorithm over multiple samples.
            * ``sublattice`` solver is based on the synchronous sublattice algorithm of Shim and Amar. 

        * ``domain`` specifies the lattice size ``[x-length, y-length]`` processed by a single process (or a group of threads) in a parallel algorithm (``sublattice``). If a parallel algorithm is not used, ``domain`` is ignored.      
        * ``random seed`` specifies the random seed used in the KMC simulations.
        * KinCat uses a random number generator and pre-generates an array of random numbers and ``random pool size`` indicates the array size. The minimum random number array is ``2*(# of subdomains)*(max number of KMC iterations per kernel launch)``.
        * The simulation will complete either after it meets ``max number of kmc kernel launches`` or the simulation reaches ``time end``.
        * Each KMC kernel launch will complete either after the ``max number of kmc steps per kernel launch`` or the ``dt`` is reached. Note that if a kernel completes before reaching ``dt``, the limiting ``dt`` will not be updated and the next kernel will only proceed until the prior ``dt`` is reached. This helps maintain regular output intervals and explains why 'extra' kernels may be desirable.
        * ``time range`` specifies the ``[ time begin, time end, dt (time-increment)]``.

For both the dump and statistics objects, giving a filename of "none" will result in no output file of that type. Filenames must either use ".json" or ".hdf5" extensions. If KinCat is built with HDF5 functionality, it will detect which file extension is used and output that file format. 

Note that a single ``solver.advance()`` function runs until it reaches ``dt`` time step or the maximum number of KMC iterations per kernel launch. When ``dt`` is set zero, the code runs for the specified number of KMC steps. If a user wants to ensure the ``dt`` time step is reached for each kernel launch, then the number of KMC steps needs to be sufficiently large. If the number of steps is not large enough to reach the desired ``dt`` timestep, then another kernel will launch with the same limiting endtime as the prior kernel. For the ``serial-rta``, running the code setting without the ``dt`` constraint does not cause any simulation errors. 

The sublattice algorithm as implemented requires that the lattice size be an even integer multiple of sublattice domain size in both x and y dimensions. Also, the sublattice domain size should be larger than 2X the interaction range of the lattice model (included in the dictionary file) to avoid potential corruptions. The ``sublattice`` algorithm requires synchronization among subdomains. A subset of subdomains are evolved simultaneously while the others are frozen. The algorithm rotates through subsets until all subdomains are synchronized at ``dt``, rejecting the final steps that would extend past ``dt``. The timestep ``dt`` should be sufficiently small so that the kinetics within the subdomain will not be significantly different at the end than at the beginning. Otherwise it will lead to significant errors. On the other hand, shorter timesteps lead to more numerous rejected events and reduced efficiency. KinCat will output a warning if a significant number of events occur before timestep ``dt``. However, this warning level is somewhat arbitrary. The user is strongly encouraged to carefully consider what parameters will balance errors and efficiency for their system. 


Ensemble Input
--------------

To exploit ``kincat-batch.x``, it is required to append an 'ensemble' section to the above main input. We use the input of Example 4 to illustrate it.  

.. code-block:: javascript
	
   "ensemble" : {
        "number of samples" : 4,
        "solver random number variations" :{
            "apply" : "enabled",
        },
        "sites random variations" : {
            "apply" : "enabled",
            "random fill ratio" : {
                "0" : [0.3, 0.5],
                "1" : [0.4, 0.4],
                "2" : [0.5, 0.3],
                "3" : [0.6, 0.2]
            }
        },
        "rates variations" : {
            "apply" : "enabled",
            "type" : "file",
            "processes" :  ["CO_ads_cus", "CO_ads_br"],
            "process rates" : {
                "0" : [1.85e+06 , 2.15e+06 ],
                "1" : [1.90e+06 , 2.10e+06 ],
                "2" : [1.95e+06 , 2.05e+06 ],
                "3" : [2.00e+06 , 2.00e+06 ]
                },
            "process instances" : [0 , 1],
            "instance rates" : {
                "0" : [1.85e+06 , 2.15e+06 ],
                "1" : [1.90e+06 , 2.10e+06 ],
                "2" : [1.95e+06 , 2.05e+06 ],
                "3" : [2.00e+06 , 2.00e+06 ]
                },
            "override filename" : "../input-rates-override-RuO2.json"   
            }
        }   


* ``ensemble`` object includes the following items:  

    * ``number of samples`` specifies the number of concurrent simulations to be run.
    * When ``solver random number variations`` is enabled, each sample uses a different sequence of random numbers to select KMC events. 
    * ``site random variations`` varies the initial configurations of samples according to the random fill ratios.

        * ``random fill ratio`` contains dictionary objects used to define the species fill ratios as in the sites object, but with each sample uniquely defined. If ``site random variations`` is enabled but the ``random fill ratio`` object is not present, each sample will initialize with a different initial configuration, but with the fill ratio set included in the 'sites' object.  
        * When this option is ``disabled``, samples use the same initial configuration, which is specified in the ``sites`` object.

    * ``process rates variations`` allows for samples to use different process rates when it is enabled.

        * ``type`` can be either ``inlined`` or ``file``.

	       * ``inlined`` looks for ``processes`` and ``process instances``, either or both or which may be included.
           * ``file`` type value take ``override filename`` keyword to load the user specified rates for samples. The format of the file should be a json file with the ``processes`` and/or ``process instances`` arrays and their corresponding rates dictionaries.

        * ``processes`` keyword indicates an array of the processes which rates are to be modified. If the ``processes`` keyword is present, then the ``process rates`` must also be present.
        * The ``process rates`` object is a dictonary with the sample index and an array of rates. The rates correspond to the processes listed in the ``processes`` array. There must be the same number of rates provided for each sample as processes listed. 
        * ``process instances`` keyword indicates an array of the process instances which rates are to be modified. If the ``process instances`` keyword is present, then the ``instance rates`` must also be present.
        * The ``instance rates`` object is a dicitonary with the sample index and an array of rates. The rates correspond to the process instances listed in the ``process instances`` array. There must be the same number of rates provided for each sample as processes instances listed.  

Note that the ``batch-rta`` solver type is the only solver that supports the ensemble section for now. The batch input should be used with ``kincat-batch.x`` executable. Also note that if a configuration file is provided in the ``sites`` object, the samples will be initialized from there rather than any fill ratios provided in the ``ensemble`` object. This configuration file may include either the same number of samples as the current batch, or may include only one, in which case all samples will be initialized with the same configuration.
   

Dictionary Input
----------------
      
This file contains the information needed to set up the lattice shape and define possible KMC events. It is generated as an output of KinCatPy and should not need to be modified by the user. This is an explanation of the information contained in the dictionary file. KinCatPy can be used to generate a dictionary file for three use-cases. KinCat uses so-called 1) 'Sets', 2) 'Uniconfig', and 3) full symmetry ('Fullsym') configurations specifying the event mechanism. The styles give statistically equivalent results, but use differing reliance on symmetry operations which changes how the processes and configurations are catalogued and retrieved. The Fullsym case specifies each possible process in a single configuration, while the sets and uniconfig cases use the symmetry of the lattice to reduce the size of the configuration(s) catalogued. The following example dictionary file specifies the sets case, which divides the processes by the sites they involve and define multiple configurations. 
 	
.. code-block:: javascript

        {"crystal": {
            "edge vectors": [[6.43, 0.0], [0.0, 3.12]], 
            "basis vectors": [[0.0, 0.0], [0.5, 0.0]], 
            "symmetry operations": {
                "shape": [4, 6], 
                    "data": [1, 0, 0, 1, 0.0, 0.0, -1, 0, 0, -1, 0.0, 0.0, -1, 0, 0, 1, 0.0, 0.0, 1, 0, 0, -1, 0.0, 0.0]
                }
            }, 
            "configurations": {
                "interaction range": [2, 2], 
                "site coordinates": [[0.0, 0.0], [0.5, 0.0], [0.0, 1.0], [0.5, 1.0]], 
                "variant orderings": [[0, 1, 2, 3], [2, 3, 0, 1]], 
                "symmetry orderings": {"0": [[0], [0]], "1": [[0], [0]], "2": [[0, 1], [1, 0]], ...}, 
                "configuration_sets": [0, 3, 6, 12, 18, 27], 
                "shape": [27, 4], 
                "data": [-1, 0, -1, -1, -1, 1, -1, -1, -1, 2, -1, ...]}, 
                "process dictionary": {
                "processes": ["CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_br_br", ...], 
                "process constraints": [[[1, 0, 2]], [[0, 0, 2]], [[1, 0, 1], [3, 0, 1]], [[0, 0, 1], ...], 
                "process symmetries": [[0], [0], [0], [0], [0, 1], [0], [0], [0], [0], [0, 1], [0, 1], ...], 
                "shape": [22, 3], 
                "data": [0, 2, 0, 2, 0, 5, 3, 5, 1, 5, 3, 6, 6, 9, 2, 7, 7, 14, 8, 8, ...]
            }
        }

* ``crystal`` object describes a unit cell structure and its symmetry operations.
    
    * ``edge vectors`` includes two vectors that form a unit cell (parallelogram) that is repeated to tile the lattice domain.
    * ``basis vectors`` represents the position of sites in the unit cell coordinates.
    * ``symmetry operations`` provides rotation matrices and translation vectors that produce equivalent symmetry configurations.

        * ``shape`` indicates ``[ # of symmetries, array size (4 entries for 2x2 matrix, 2 entries for 2x1 vector) ]`` and is used to interpret the ``data`` array. 

* ``configurations`` object includes a list of possible configurations. A configuration is defined as a unique arrangement of simulation species within the interaction range of a system process.

    * ``interaction range`` defines the range around the central site that needs to be recalculated after each event due to possible changes to processes and rates in that region. It is given in unit cells. It also affects the minimum domain size for parallel solvers. 
    * ``site coordinates`` includes the position of sites in the reference configuration. The position is given in units of the crystal edge vectors. 
    * ``variant orderings`` and ``symmetry orderings`` relate to the possible enumerations (or re-ordering of sites) forming symmetry equivalent configurations when the configuration is mapped to the lattice sites.
    * ``configuration_sets`` stores the divisions between different configuration definitions in the catalogue.
    * The list of configurations is stored as a 2D array ``[ # of configurations, configuration size ]`` where the entries of ``data`` are the species index. The sets dictionary style may result in some site occuptations being specified as '-1'. This indicates that this site is not important for that configuration definition. 


* ``process dictionary`` object describes process mapping from one configuration to the other configuration.

  * ``processes`` is a list of process labels. For the reduced symmetry configurations, the processes are unique.
  * ``process constraints`` is a rank-3 array ``[ # of processes, # of constraints, constraint size(3) ]``. For a corresponding process, a constraint ``[ site index, initial species, final species]`` represents the initial and final species of the specified site. The specified site in the initial and final configurations should match those given in the constraint. Otherwise, we do not consider it as a valid process instance. The initial and final species may be the same, indicating that the site is not changed by the process, but that it is a 'bystander' site important for the process definition but not changed by the process. 
  * ``process symmetries`` indicates a pattern index that should be accounted when computing valid events. For example, an absorption event can be counted multiple times in the reduced symmetry configurations. To prevent this multiple counting, the process symmetry information include ``[0]`` pattern index so that only the first symmetry pattern is used when searching for possible events. Without the process symmetry information, the KMC process will find multiple events that are equivalent (for the adsorption example, it would find four events: one for each symmetry operation). To avoid the duplicated event search, we can also use the Fullsym dictionary style as explained in the next section.
  * A process instance is specified as ``[initial configuration, final configuration, process index]`` and stored as rank-2 array ``[# of events, event size(3)]``.

A Fullsym input is shown below. Note that the event dictionary grows exponentially with the number of sites and the number of species. Using the Fullsym or even Uniconfig dictionary style might be prohibitive for a large reaction model.     
	
.. code-block:: javascript
   
   {
        "crystal": {
        "edge vectors": [[6.43, 0.0], [0.0, 3.12]], 
        "basis vectors": [[0.0, 0.0], [0.5, 0.0]], 
        "symmetry operations": {
            "shape": [1, 6], 
            "data": [1, 0, 0, 1, 0.0, 0.0]
        }
    }, 
    "configurations": {
        "interaction range": [2, 2], 
        "site coordinates": [[-0.5, 0.0], [0.0, -1.0], [0.5, -1.0], [0.0, 0.0], [0.5, 0.0], [0.0, 1.0], [0.5, 1.0], [1.0, 0.0]], 
        "variant orderings": [[0, 1, 2, 3, 4, 5, 6, 7]], 
        "symmetry orderings": {"0": [[0, 1, 2, 3, 4, 5, 6, 7]]}, 
        "configuration_sets": [0, 6561], 
        "shape": [6561, 8], 
        "data": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, ...]
    }, 
    "process dictionary": {
    "processes": ["CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_br_br", ...], 
    "process constraints": [[[4, 0, 2]], [[3, 0, 2]], [[4, 0, 1], [6, 0, 1]], ...], 
    "process symmetries": [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], ...], 
    "shape": [32076, 3], 
    "data": [0, 54, 0, 0, 162, 1, 0, 30, 2, 0, 90, 3, 0, 108, 4, 0, 2268, 5, 1, 55, 0, 1, ... }}
   }

Here, we explain the major differences from the Sets case.   

* The Fullsym style, which does not use symmetries, has an identity matrix for ``symmetry operations``.
* The ``variant`` ordering is trivial i.e., identity map.
* The ``configuration_sets`` is trivial. 
* As expected, the number of possible configurations and process intances is much bigger than the Sets case e.g., 6561 vs 27 and 32076 vs 22.
* ``events`` can have duplicated processes names e.g., same absorption process with different initial configurations.
* ``process constraints`` and ``process symmetries`` require inputs for the same number of ``processes`` array size.
* ``process symmetries`` is trivial and the Fullsym case does not have multiple event searches.

Rate Input
----------
  
The rate input is explained with the sample script below.

.. code-block:: javascript
   
    {
        "default rate": 0, 
        "process specific rates" : {
        "CO_ads_cus" : 2.04E+06 ,
        "CO_ads_br" : 2.04E+06,
        "O_ads_cus_cus" : 3.81E+03,
        "O_ads_br_br" : 3.81E+03,
        "O_ads_br_cus" : 3.81E+03,
        "CO_des_cus" : 1.82E+07,
        "CO_des_br" : 5.50E+04,
        ...
        },
        "event specific rates": {
        "5" : 1.85E+07,
        "11" : 2.00E+06,
        ...
        }
    }

* ``default rate`` is used when the rate is not otherwise specified. A negative default value can be used for error checking if the input file must specify all processes (or process instances).
* ``process specific rates`` includes rates for processes. 
    * All process instances with the same process will be set to the same rate.
    * Events with processes not specified will be set to the default rate.

* ``process instance specific rates`` includes rates for specific instances. 
    * Process instances not specified will be set to the process rate if one was specified, or the default rate otherwise. 

Dump Output
-----------

When ``dump`` is enabled from the main input, the code dumps the output of the sites in the following format. 

.. code-block:: javascript
		
   {
    "samples" : 4,
    "number of species": 3,
    "coordinates": {
        "shape": [ 800, 2 ], 
        "data": [ 0, 0, 3.215, 0, 0, 3.12, 3.215, 3.12, 0, 6.24, ...]
    },
    "sites": [
        { 
            "sample": 0,
            "time": 0,
            "data": [ 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 1, ...]
        }, 
        { 
            "sample": 1,
            "time": 0,
            "data": [ 1, 2, 2, 0, 2, 1, 0, 2, 1, 2, 0, 1, 2, 2, ...]
        }, 
        { 
            "sample": 2,
            "time": 0,
            "data": [ 1, 2, 1, 1, 0, 2, 1, 2, 2, 2, 1, 1, 0, 1, ...]
        }, 
        { 
            "sample": 3,
            "time": 0,
            "data": [ 1, 1, 1, 1, 1, 1, 0, 0, 1, 2, 1, 1, 1, 0, ...]
        }, 
        { 
            "sample": 0,
            "time": 5.00036e-07,
            "data": [ 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 0, ...]
        }, 
        { 
            "sample": 1,
            "time": 5.00061e-07,
            "data": [ 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 0, 2, 2, ...]
        }, 
        { 
            "sample": 2,
            "time": 5.00135e-07,
            "data": [ 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, ...]
        }, 
        { 
            "sample": 3,
            "time": 5.00028e-07,
            "data": [ 1, 1, 1, 2, 1, 2, 2, 0, 1, 2, 1, 1, 1, 2, ...]
        }, 
        ...
        ]
    }

A dump file can be used for post-processing and it includes 1) number of species, 2) coordinates, and 3) time series of sites information. The time series information includes the sample number (even if not running an ensemble), the time stamp, and the current occupation state at each site. Similar information is given with hdf5 style outputs. An example of post-processing is given with each example. Additionally, ``restart-sites.json`` and ``dump-batch-sites.json`` files are created to record the last site configurations when the code completes, which can be used for restarting the simulation. 

Stats Output
------------

When ``statistics`` is enabled from the main input, the code dumps the selected statistics in the following format.

.. code-block:: javascript

    {
        "lattice size" : [ 20, 20, 2 ], 
    "samples" : 4,
    "number of species" : 3,
    "number of processes" : 22,
    "processes" : [ "CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_br_br", ... ], 
    "readings" : [ 
        { 
            "sample" : 0,
            "time" : 0,
            "species coverage" : [ 0.2, 0.3, 0.5 ],
            "process counts" : [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            "site species coverage" : [[ 0.22, 0.29, 0.49 ], [ 0.18, 0.31, 0.51 ]]
        }, 
        { 
            "sample" : 1,
            "time" : 0,
            "species coverage" : [ 0.2, 0.4, 0.4 ],
            "process counts" : [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            "site species coverage" : [[ 0.22, 0.4175, 0.3625 ], [ 0.18, 0.3825, 0.4375 ]]
        }, 
        { 
            "sample" : 2,
            "time" : 0,
            "species coverage" : [ 0.2, 0.5, 0.3 ],
            "process counts" : [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            "site species coverage" : [[ 0.22, 0.48, 0.3 ], [ 0.18, 0.52, 0.3 ]]
        }, 
        { 
            "sample" : 3,
            "time" : 0,
            "species coverage" : [ 0.2, 0.6, 0.2 ],
            "process counts" : [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            "site species coverage" : [[ 0.19, 0.625, 0.185 ], [ 0.21, 0.575, 0.215 ]]
        }, 
        { 
            "sample" : 0,
            "time" : 4.00134e-07,
            "species coverage" : [ 0.0425, 0.21125, 0.74625 ],
            "process counts" : [ 2352, 160, 1, 2, 6, 2222, 4, 0, 0, 0, 0, 103, 0, 0, 0, 4, 0, 0, 11, 0, 0, 78 ],
            "site species coverage" : [[ 0, 0.315, 0.685 ], [ 0.085, 0.1075, 0.8075 ]]
        }, 
        { 
            "sample" : 1,
            "time" : 4.00335e-07,
            "species coverage" : [ 0.04125, 0.2925, 0.66625 ],
            "process counts" : [ 2153, 166, 3, 1, 4, 2001, 3, 0, 0, 0, 0, 94, 0, 0, 0, 8, 0, 0, 19, 0, 0, 83 ],
            "site species coverage" : [[ 0.005, 0.4325, 0.5625 ], [ 0.0775, 0.1525, 0.77 ]]
        }, 
        { 
            "sample" : 2,
            "time" : 4.00002e-07,
            "species coverage" : [ 0.0425, 0.37375, 0.58375 ],
            "process counts" : [ 1876, 196, 1, 0, 9, 1716, 8, 0, 0, 0, 0, 77, 0, 0, 0, 9, 0, 0, 11, 0, 1, 109 ],
            "site species coverage" : [[ 0.0025, 0.5, 0.4975 ], [ 0.0825, 0.2475, 0.67 ]]
        }, 
        { 
            "sample" : 3,
            "time" : 1.60017e-06,
            "species coverage" : [ 0.03375, 0.395, 0.57125 ],
            "process counts" : [ 8205, 202, 7, 1, 8, 7901, 13, 0, 0, 0, 0, 61, 0, 0, 0, 36, 0, 0, 73, 0, 0, 123 ],
            "site species coverage" : [[ 0, 0.65, 0.35 ], [ 0.0675, 0.14, 0.7925 ]]
        }, 
        ...
        ]
    }

In the ``readings`` object, the ``sample`` and ``time`` will always be present. However, the rest of the objects will depend on the selected options. The ``species coverage`` gives current fill fraction of each species over all sites. The ``site species coverage`` breaks this fill fraction down by lattice site. The ``process counts`` gives the count of each process since the start of the simulation, which can be used to calculate turn-over frequencies. 
    
.. autosummary::
   :toctree: generated
	     
