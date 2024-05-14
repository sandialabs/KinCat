This directory has three examples. The first is of a run using the serial solver. It is run using the command:
./path_to_executable/kincat.x --input='ex1-input.json'

The second example is using the sublattice solver. It is run using the command: 
./path_to_executable/kincat.x --input='ex2-input.json'

Note that the sublattice solver requires the additional 'domain' keyword.

The third example is using hdf5 formatting for the output files. It is run using the command:
./path_to_executable/kincat.x --input='ex3-input.json'

Python files to plot results for these simulations are also included. They can be run by using the command:
python plot-ex*.py
