File Format
============

KinCat requires three input files specifying 1) solver, 2) dictionary, and 3) rates. The solver input can be specified by a user to control the KMC simulations while the dictionary input is auto-generated from a pre-processor, KinCatPy. The rates input is separated as we may want to use/test different set of rates that can be computed by multiple sources e.g., RMG.  

Main Input
----------

The following json input is the main input file extracted from the RuO2 example.

.. code-block:: javascript

   {
       "kincat": {
           "dictionary" : "${KINCAT_INSTALL_PATH}/example/RuO2/input-dictionary-RuO2.json",
           "rates" : "${KINCAT_INSTALL_PATH}/example/RuO2/input-rates-RuO2.json",
           "sites" : {
               "lattice" : [ 10, 10, 2 ],
               "type" : "random",
               "random seed" : 12345,
               "random fill ratio" : [0.1, 0.5],
               "filename" : "dump-sites.json"
           },
           "dump" : {
                "dump filename" : "dump.json",
                "dump interval" : 5E-7
             },
           "statistics" : {
                "stats filename" : "stats.json",
                "types" : ["species_coverage","process_counts"],
                "stats interval" : 4E-7
            },
           "solver" : {
               "type" : "sublattice",
               "domain" : [ 5, 5 ],
               "random seed" : 13542,
               "random pool size" : 10000,
               "max number of kmc kernel launches" : 100,
               "max number of kmc steps per kernel launch" : 50,
               "time range" : [ 0, 1000, 1E-7 ]
           }
       }
   }

* ``kincat`` object includes the following items:  
    * ``dictionary`` specifies the file path to the dictionary json file.
    * ``rates`` specifies the file path to the rate json file.
    * ``sites`` object includes lattice information and its initialization method.
    
        * ``lattice``: 2D lattice object. The array input implies ``[ x-length, y-length, n-basis ]``. The lengths are specified in unit cells, and must be integers.
        * ``type``: lattice initialization method.    
            * ``random`` value fills the lattice randomly with percentages specified by ``random fill ratio`` keyword; otherwise, ``random fill ratio`` is ignored. ``random fill ratio`` should be an array with as many entries as the number of species (excluding vacant sites), and in the same order as they were defined in KinCatPy. The sum of these fill ratios should not exceed 1.0. A ``random seed`` is given just for the initialization routine.
            * ``file`` value loads site specification from a file given by the ``filename``; otherwise, ``filename`` is ignored.
	
    * ``dump`` object relates to an output file with the full simulation state. 
        * ``dump filename`` object specifies the name of the dump output file. 
        * ``dump interval`` object indicates the minimum time between outputs of the full simulation state.  
    * ``statistics`` object relates to an output file of various desired statistics of the KMC simulation.
        * ``filename`` object specifies the name of the statistics output file.
        * ``types`` object is a list of the various statistics desired by the user. Further types may be added later, or added by the user by modifying the Stats Class. Currently implemented types are:
            * ``"species_coverage"`` outputs the fraction of the lattice sites filled by each species. 
            * ``"process_counts"`` outputs the number of times each process has occurred up to that point in the KMC simulation.
            * ``"site_species_coverage"`` outputs the fraction of each site in the lattice unit cell filled by each species. Thus it gives similar data to ``species coverage`` but broken out by site. Note that it breaks out the data by each site, not each unique site type. This is important if the unit cell has multiple symmetric sites. 
        * ``stats interval`` object indicates the minimum time between outputs of the statsitics listed in ``types``. 
    * ``solver`` object include the following items to specify the KMC simulation.
        * ``type`` specifies the algorithm selected i.e., ``serial-rta``, ``sublattice``, ``batch-rta``.
        * ``domain`` specifies the lattice size ``[x-length, y-length]`` processed by a single process (or a group of threads) in a parallel algorithm (``sublattice``).
      
        * ``serial-rta`` solver type requires that the domain size matches to the lattice size; thus, this domain information is ignored. This solver is based on the BKL or 1st Order KMC algorithm. The ``batch-rta`` solver implements this algorithm for multiple samples.
        * ``sublattice`` solver is based on the synchronous sublattice algorithm of Shim and Amar. It requires satisfing the following conditions for appropriate domain setup:
 	        * The lattice size should be an even integer multiple of sublattice domain size in both x and y dimensions.
 	        * The sublattice domain size should be bigger than the interaction range of the lattice model (included in the dictionary file) to avoid potential corruptions.
	    
    * ``random seed`` specifies the random seed used in the KMC simulations.
    * KinCat uses a random number generator and pre-generates an array of random numbers and ``random pool size`` indicates the array size. The minimum random number array is ``2*(# of subdomains)*(max number of KMC steps per kernel launch)``.
    * The simulation will complete either after it meets ``max number of kmc kernel launches`` or the simulation reaches ``time end``.
    * Each KMC kernel launch will complete either after the ``max number of kmc steps per kernel launch`` or the ``dt`` is reached. Note that if a kernel completes before reaching ``dt``, the next kernel will only proceed until the prior ``dt`` is reached. 
    * ``time range`` specifies the ``[ time begin, time end, dt (time-increment)]``.

Note that a single ``solver.advance`` function runs until it reaches ``dt`` time step or the maximum number of KMC steps per kernel launch. When ``dt`` is set zero, the code runs for the specified number of KMC steps. If a user wants to ensure the ``dt`` time step is reached for each kernel launch, then the number of KMC steps needs to be sufficiently large. If the number of steps is not large enough to reach the desired ``dt`` timestep, then another kernel will launch with the same limiting endtime as the original ``dt`` timestep. For the ``serial-rta``, running the code setting without the ``dt`` constraint does not cause any simulation errors. 
The ``sublattice`` algorithm requires synchronization among subdomains. A subset of subdomains are evolved simultaneously while the others are frozen. The algorithm rotates through subsets until all subdomains are synchronized at ``dt``, rejecting the final steps that would extend past ``dt``. The timestep ``dt`` should be sufficiently small so that the kinetics in the subdomain will not be significantly different at the end than at the beginning. Otherwise it will lead to significant errors. On the other hand, shorter timesteps lead to more numerous rejected events and reduced efficiency. The user is encouraged to carefully consider what parameters will balance errors and efficiency for their system.

Ensemble Input
--------------

To exploit ``kincat-batch.x``, it is required to append the following ensemble section to the above main input.  

.. code-block:: javascript
		
   {		
        "ensemble" : {
        "number of samples" : 4,
        "solver random number variations" :{
            "apply" : "enabled",
            },
        "sites random variations" : {
            "apply" : "disabled",
            "random fill ratio" : {
                "0" : [0.25 , 0.3 ],
                "1" : [0.35 , 0.4 ],
                "2" : [0.45 , 0.5 ],
                "3" : [0.1 , 0.15 ]
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
            "override filename" : "../example/RuO2/input-rates-override-RuO2.json"  
            }
        } 
   }


* ``ensemble`` object includes the following items:  
  * ``number of samples`` specifies the number of samples.
    
    * When ``solver random number variations`` is enabled, each sample uses a different sequence of random numbers selecting a KMC event. 
    * ``site random variations`` varies the site configurations of samples.
        * ``random fill ratio`` contains dictionary objects used to define the species fill ratios as in the sites object, but with each sample uniquely defined. If ``site random variations`` is enabled but the ``random fill ratio`` object is not present, each sample will initialize with a different initial configuration, but with the fill ratio set included in the 'sites' object.  
        * When this option is ``disabled``, samples use the same initial configuration, which is randomly configured.
    * ``process rates variations`` allows for samples to use different process rates when it is enabled.
        * ``type`` can be either ``inlined`` or ``file``.
	
	    * ``inlined`` looks for ``processes`` and ``process instances``. Either or both may be included.
	        * ``processes`` keyword indicates an array of the processes which rates are to be modified.
                * If the ``processes`` keyword is present, then the ``process rates`` must also be present. 
                * The ``process rates`` object is a dictonary with the sample index and an array of rates. The rates correspond to the processes listed in the ``processes`` array. 
	        * ``process instances`` keyword indicates an array of the process instances which rates are to be modified.
                * If the ``process instances`` keyword is present, then the ``instance rates`` must also be present.
                * The ``instance rates`` object is a dicitonary with the sample index and an array of rates. The rates correspond to the process instances listed in the ``process instances`` array. 
	   * ``file`` type value take ``override filename`` keyword to load the user specified rates for samples. The format of the file should be a json file with the ``processes`` and/or ``process instances`` arrays and their corresponding rates dictionaries.

Note that the ``batch-rta`` solver type is the only solver that supports the ensemble section for now. The batch input should be used with ``kincat-batch.x`` executable. 	  
   

Dictionary Input
----------------
      
This file contains the information needed to set up the lattice shape and define possible KMC events. It is generated as an output of KinCatPy and should not need to be modified by the user. This is an explanation of the information contained in the dictionary file. KinCatPy can be used to generate a dictionary file for two use-cases. KinCat uses so-called 1) reduced symmetry and 2) full symmetry configurations specifying the event mechanism. The full symmetry case specifies each possible process, while the reduced symmetry case uses the symmetry of the lattice to reduce the amount of information needed to be stored in the dictionary. The following example dictionary file specifies the reduced symmetry case. 
 	
.. code-block:: javascript
   
   { 
        "crystal": {
            "edge vectors": [[6.43, 0.0], [0.0, 3.12]], 
            "basis vectors": [[0.0, 0.0], [0.5, 0.0]], 
            "symmetry operations": {
                "shape": [4, 6], 
                "data": [-1, 0, 0, -1, 0.0, 0.0, -1, 0, 0, 1, 0.0, 0.0, ...]
            }
        }, 
        "configurations": {
            "site coordinates": [[0.0, 0.0], [0.5, 0.0], [0.0, 1.0], [0.5, 1.0]], 
            "variant orderings": [[2, 3, 0, 1], [0, 1, 2, 3]], 
            "interaction range": [3, 3], 
            "shape": [45, 4], 
            "data": [1, 1, 1, 1, 1, 2, 1, ...]
        }, 
        "process dictionary": {
            "processes": ["CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_br_br", ...], 
            "process constraints": [[[1, 0, 2]], [[0, 0, 2]], [[1, 0, 1], [3, 0, 1]], ...], 
            "process symmetries": [[0], [0], [0], [0, 1], [0, 2], [0], [0], [0, 1], ...], 
            "shape": [258, 3], 
            "data": [0, 17, 7, 0, 39, 8, 0, 8, 9, 1, 40, 8, 1, 16, 9, 1, 2, 5, ...]
        }
    }

* ``crystal`` object describes a unit cell structure and its symmetry operations.

  * ``edge vectors`` includes two vectors that form a unit cell (parallelogram) that is repeated to tile the lattice domain.
  * ``basis vectors`` represents the position of sites in the unit cell coordinates.
  * ``symmetry operations`` provides rotation matrices and translation vectors that produce equivalent symmetry configurations.
    
    * ``shape`` indicates ``[ # of symmetries, array size (4 entries for 2x2 matrix, 2 entries for 2x1 vector) ]`` and is used to interpret the ``data`` array.

* ``configurations`` object includes a list of possible configurations. A configuration is defined as a unique arrangement of simulation species within the interaction range of a system process. 

  * ``site coordinates`` includes the position of sites in the reference configuration. The position is given in units of the crystal edge vectors. 
  * ``variant orderings`` represents the possible enumerations (or re-ordering of sites) forming symmetry equivalent configurations when the configuration is mapped to the lattice sites.
  * ``interaction range`` defines the range around the central site that needs to be recalculated after each event due to possible changes to processes and rates in that region. It is given in unit cells. It also represents the minimum domain size for parallel solvers. 
  * The list of configurations is stored as a 2D array ``[ # of configurations, configuration size ]`` where the entries of ``data`` are the species index.


* ``process dictionary`` object describes process mapping from one configuration to the other configuration. 
  
  * ``processes`` is a list of process labels. For the reduced symmetry configurations, the processes are unique.
  * ``process constraints`` is a rank-3 array ``[ # of processes, # of constraints, constraint size(3) ]``. For a corresponding process, a constraint ``[ site index, initial species, final species]`` represents the initial and final species of the specified site. The specified site in the initial and final configurations should match those given in the constraint. Otherwise, we do not consider it as a valid process instance. The initial and final species may be the same, indicating that the site is not changed by the process, but that it is important for the process definition. 
  * ``process symmetries`` indicates a pattern index that should be accounted when computing valid events. For example, an absorption event can be counted multiple times in the reduced symmetry configurations. To prevent this multiple counting, the process symmetry information include ``[0]`` pattern index so that the first symmetry pattern is only used when searching for possible events. Without the process symmetry information, the KMC process will find multiple events that are equivalent (for the adsorption example, it would find four events: one for each symmetry operation). To avoid the duplicated event search, we can also use the full symmetry dictionary as explained in the next section.
  * A process instance is specified as ``[initial configuration, final configuration, process index]`` and stored as rank-2 array ``[# of events, event size(3)]``.

A full symmetry input is shown below. Note that the event dictionary grows exponentially with the number of sites and the number of species. Using a full symmetry dictionary might be prohibitive for a large reaction model.     
	
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
        "site coordinates": [[-0.5, 0.0], [0.0, -1.0], [0.5, -1.0], [0.0, 0.0], [0.5, 0.0], ...], 
        "variant orderings": [[0, 1, 2, 3, 4, 5, 6, 7]], 
        "interaction range": [3, 3], 
        "shape": [6561, 8], 
        "data": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, ...]
	   },
   "process dictionary": {
      "processes": ["CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_cus_cus", ...],
      "process constraints": [[[4, 0, 2]], [[3, 0, 2]], [[4, 0, 1], [6, 0, 1]], ...],
      "process symmetries": [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], ...],
      "shape": [34992, 3],
      "data": [0, 540, 9, 0, 60, 10, 0, 1620, 11, 0, 180, 12, 0, 4536, 13, ...]
      }
   }

Here, we only explain the major difference from the reduced symmetry case.   

* The full symmetry configuration has an identity matrix for ``symmetry operations``.
* The ``variant`` ordering is trivial i.e., identity map.
* As expected, the number of possible configurations and process intances is much bigger than the reduced configuration case e.g., 6561 vs 45 and 34992 vs 258.
* ``events`` can have duplicated processes names e.g., same absorption process with different initial configurations.
* ``process constraints`` and ``process symmetries`` require inputs for the same number of ``processes`` array size.
* ``process symmetries`` is trivial and the full symmetry case does not have the duplicated event search issue.

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
       "number of species": 3,
       "coordinates": {
            "shape": [ 200, 2 ],
            "data": [ 0, 0, 3.215, 0, 0, 3.12, 3.215, 3.12, 0, 6.24 ... ]
       },
       "sites": [
           {
	       "sample": 0,
	       "time": 0,
	       "data": [ 0, 0, 0, 2, 2, 0, 2, 0, 0, 1, 0, 0 ... ]
	   }
           {
	       "sample": 1,
	       "time": 0,
	       "data": [ 0, 1, 0, 1, 2, 0, 2, 0, 0, 1, 0, 0 ... ]
	   }
           {
	       "sample": 0,
	       "time": 0.1,
	       "data": [ 1, 0, 0, 2, 2, 0, 2, 0, 0, 1, 0, 0 ... ]
	   }
           {
	       "sample": 1,
	       "time": 0.15,
	       "data": [ 0, 2, 0, 2, 2, 0, 2, 2, 0, 1, 0, 0 ... ]
	   }
       ]
   }

A dump file can be used for post-processing and it includes 1) number of species, 2) coordinates, and 3) time series of sites information. An example of post-processing is illustrated in ``${KINCAT_INSTALL_PATH}/bin/plot-dump.ipynb``. Additionally, ``dump-sites.json`` and ``dump-batch-sites.json`` files are created to record the last site configurations when the code completes, which can be used for restarting the simulation.

Stats Output
------------

When ``statistics`` is enabled from the main input, the code dumps the selected statistics in the following format.

.. code-block:: javascript

    {
        "number of species": 3,
        "number of processes": 22,
        "processes": [ "CO_ads_cus", "CO_ads_br", "O_ads_cus_cus", "O_ads_br_br", ...], 
        "readings" : [ 
            { 
                "sample": 0,
                "time": 0,
                "species coverage": [ 1, 0, 0 ],
                "process counts": [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...]
            }, 
            { 
                "sample": 0,
                "time": 4.02101e-07,
                "species coverage": [ 0.705, 0, 0.295 ],
                "process counts": [ 75, 54, 0, 0, 0, 67, 3, 0, 0, 0, 0, 1751, 0, 0, 0, 0, 0, ...]
            }, 
            { 
                "sample": 0,
                "time": 8.00534e-07,
                "species coverage": [ 0.59, 0, 0.41 ],
                "process counts": [ 159, 74, 0, 0, 0, 146, 5, 0, 0, 0, 0, 4116, 0, 0, 0, 0, ...]
            }
            ...
        ]
    }

In the ``readings`` object, the ``sample`` and ``time`` will always be present. However, the rest of the objects will depend on the selected options. 
    
.. autosummary::
   :toctree: generated
	     
