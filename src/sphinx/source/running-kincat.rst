Running KinCat
==============

We can run the ``kincat.x`` executable in the ``${KINCAT_INSTALL_PATH}/core`` directory. The following command line options are available for the current version. There are three verbose options to show additional run-time information. The ``verbose-parse`` shows raw input from the input json file. Note that configurations and events specification can be different in the dictionary and within the code as KinCat requires a sorted list of configurations.The ``verbose`` option shows the object details used in KinCat. The ``verbose-iterate`` option shows the KMC algorithm details. Each verbosity type has four levels: 0) show nothing, 1) basic information, 2) reserved for users verbose comments, 3) lengthy details, 4) show everything.

.. code-block:: bash

   cd ${KINCAT_INSTALL_PATH}/bin
   ./kincat.x --help
    Usage: ./kincat.x [options]
      options:
      --echo-command-line           bool      Echo the command-line but continue as normal
      --help                        bool      Print this help message
      --input                       string    Input dictionary file name
                                              (default: --input=input.json)
      --verbose                     int       Verbosity level
                                              (default: --verbose=0)
      --verbose-iterate             int       Verbosity level in KMC iteration
                                              (default: --verbose-iterate=0)
      --verbose-parse               int       Verbosity level in parsing step (not yet sorted)
                                              (default: --verbose-parse=0)
    Description:
       KinCat Main

The source files include several examples. We show output for example 1 here. By default, we use Kokkos OpenMP as the execution space for CPU processors while this example uses the serial residential time algorithm. Thus, users may want to set ``export OMP_NUM_THREADS=1`` not to carry some Kokkos overhead of using multiple threads. However, this is unlikely to be significant for this limited run. 

.. code-block:: bash     

   cd ${KINCAT_INSTALL_PATH}/bin		
   ../../../../kincat_build/core/kincat.x --input='ex1-input.json' --verbose-iterate=1
    Kokkos::OpenMP::initialize WARNING: OMP_PROC_BIND environment variable not set
    In general, for best performance with OpenMP 4.0 or better set OMP_PROC_BIND=spread and OMP_PLACES=threads
    For best performance with OpenMP 3.1 set OMP_PROC_BIND=true
    For unit testing set OMP_PROC_BIND=false
    -- Kokkos::OpenMP is used
    solver_type = serial-rta
    -- Lattice : Lattice
      -- number of species : 3
      -- domain : [10, 10, 2]
      -- edge vectors : 
      [6.43,0]
      [0,3.12]
      -- basis : 
      [0,0]
      [0.5,0]

    Reset constraints (due to sets construction): 22, 2
      instance: 0, 0, [0, 0]
      instance: 1, 5, [0, 0]
      instance: 2, 1, [0, 0]
      instance: 3, 6, [0, 0]
      instance: 4, 2, [0, 0]
      instance: 5, 14, [-1, 0]
      instance: 6, 10, [-1, 0]
      instance: 7, 7, [0, 0]
      instance: 8, 18, [-1, 0]
      instance: 9, 3, [0, 0]
      instance: 10, 15, [-1, 0]
      instance: 11, 11, [-1, 0]
      instance: 12, 8, [0, 0]
      instance: 13, 19, [-1, 0]
      instance: 14, 4, [0, 0]
      instance: 15, 16, [0, 0]
      instance: 16, 12, [0, 0]
      instance: 17, 17, [0, 0]
      instance: 18, 9, [0, 0]
      instance: 19, 20, [0, 0]
      instance: 20, 13, [0, 0]
      instance: 21, 21, [0, 0]
    n_process_types found : 22
    -- ProcessDictionary : Dictionary
      -- variant ordering : (2, 20)
      -- configuration list : (27, 4)
      -- process instance list : (22, 3), constraints: (22, 2), rates : (22)
    -- ProcessDictionary : Instance details
      -- process instance list : (22), process list : (22)
    -- Solver : serial-rta
    -- Dump : Dump
      -- filename dump : "dump-ex1.json"
      -- filename restart: "restart_sites.json"
    creating json file : dump-ex1.json
    -- Stats : Stats
      -- filename stats : "stats-ex1.json"
    epoch = 0, t = 5.00028e-06
      -- # of events occured : 16889
    epoch = 1, t = 1.00014e-05
      -- # of events occured : 33363
    ...

This KMC runs and produces the output files ``dump-ex1.json`` and ``stats-ex1.json``. We can post-process these by using ``python plot-ex1.py``.

Solving for Multiple Samples
----------------------------

For small and medium problem sizes, it may be beneficial to use a batch parallel version of kincat i.e., ``kincat-batch.x``. Examples 4 and 5 demonstrate the use of this batch version, and more details of the batch parallelism use case will be explained later with the batch input file.

.. code-block:: bash

   cd ${KINCAT_INSTALL_PATH}/bin		
   ./kincat-batch.x --input="ex4-input.json"


Running on Weaver
-----------------------

To run the code with a GPU, we first allocate an interactive compute node with following command. A single node is dedicated for the user. The Power9 CPU has 40 cores and can utilize 160 threads with symmetric multi-processing (SMP4) accelerated with 4 GPUs. Since Kokkos does not support the multi GPU use case, a user explicitly maps the MPI processes to different GPUs by adding ``--kokkos-num-devices=4``. 

.. code-block:: bash

   [weaver11] bsub -gpu num=1 -Is -q rhel8 bash
   ***Forced exclusive execution
   Job <42355> is submitted to queue <rhel7W>.
   <<Waiting for dispatch ...>>
   <<Starting on weaver1>>
   [weaver1 ~]$ ./kincat.x --input='input.json' 

Note that this will use the default GPU, but there are four GPU's available. This can be seen through the following environment variables: 

.. code-block::bash

  $ env | grep CUDA
  CUDA_VISIBLE_DEVICES1=0,1,2,3
  CUDA_VISIBLE_DEVICES=0,1,2,3

By setting these environment variables to only include a single value, e.g. '1', then that specific GPU will be assigned the job, allowing the user to use all GPU's without interference from concurrent jobs. 

The same bsub command, except without the '-gpu num=1' arguments, is used to request a CPU node on Weaver. If desired, the number of threads used can be set by setting the environment variable ``OMP_NUM_THREADS``. Note that the KinCat build is specific to whether CPU or GPU is intended to be used. Two builds are required if the user wishes to be able to use both.


Restart Simulation
------------------

When a simulation completes before it reaches its steady state, the simulation can be restarted using ``restart-sites.json`` and ``dump-batch-sites.json`` output. These files contain the last snapshot of the simulation state and they are created at the end of the simulation or when the code catches an exception. See the lattice section of the input file format explained next. These filenames are currently hard-coded into KinCat. They are only produced if a dump-style output is enabled, and are produced with a json format even if HDF5 outputs are requested. 


   
.. autosummary::
   :toctree: generated
