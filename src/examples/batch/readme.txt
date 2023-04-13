This example is of how to run a batch ensemble. It is run by executing the command:
./path_to_executable/kincat-batch.x --input='ex3-input.json'

The details for the ensemble are contained in the ensemble object. For this example, the rates variations are sourced from the input-rates-override-RuO2.json file. It therefore ignores the rates given in the ex3-input.json file. However, they are in the correct format should the rates variations 'type' be changed from 'file' to 'inlined'. 