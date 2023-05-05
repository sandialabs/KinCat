To Do
=====


1. Validation
   a. Serial RTA algorithm for realistic catalysis problems where data is public e.g., RuO2.
   
2. Kinetic Model Interface to KinCatPy
   a. For the problems of interest, find (or define) commonly used kinetic model format e.g., CHEMKIN, Cantera, RMG, etc.
   b. Automate the translation from existing input format to KinCatPy.

3. Parallel KMC
   a. Validate the sublattice algorithm against the serial RTA algorithm.
   b. Implement time warp algorithm and validate against serial RTA algorithm.
   c. Unit tests for all levels.

4. Documentation
   a. Refine this documentation.

5. UQ Interface
   a. Find problem sizes that is suited for single program multiple data (SPMD) parallelism case.
   b. Find problem sizes that runs better with ``kincat-batch.x``. 
      i. Find a good way to deal with massive data output from running KMC on many samples.
   c. Unify the basic data structure, i.e., lattice and dictionary, that can be used for the batch use case.
      i. Single problem interface is just a sub-case of batch use case where the number of samples is one.
      ii. Lattice - change site view from 1d to 2d, Dictionary - change rate view from 1d to 2d. Then, unify the interface including sample index.
	 
5. GPU simulation   
   a. Find good problem size that utilizes a GPU better.
   b. Find hotspots and optimize GPU performance on realistic application problems. 

6. Misc tasks
   a. Improve the restart option. The current restart option assumes a single simulation output without considering batch output. ``kincat-batch.x`` results ``dump-batch-site.json`` at the end of the simulation recording the last status of the lattice configuration of each sample and its time. The array of object should be properly parsed and store in the input object to assign it to lattice.
   b. Currently, I checked only RuO2 example inputs. Validate different combination of RuO2 examples and modify other example input files accordinlgy.     
      
.. autosummary::
   :toctree: generated

