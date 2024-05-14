These files show how to run a batch ensemble. Example 4 is run by executing the command:
./path_to_executable/kincat-batch.x --input='ex4-input.json'

The details for the ensemble are contained in the ensemble object. It runs a batch of four simulations that have different initial fill ratios, and have some rates modified from those provided in For this example, the rates variations are sourced from the input-rates-override-RuO2.json file. It therefore ignores the rates given in the ex4-input.json file. However, they are in the correct format should the rates variations 'type' be changed from 'file' to 'inlined'. 

Example 5 is run by executing the command: 
./path_to_executable/kincat-batch.x --input='ex5-input.json'

This batch of simulations are initialized from the ex5-restart.json file. This file was automatically outputted from a prior run. 
