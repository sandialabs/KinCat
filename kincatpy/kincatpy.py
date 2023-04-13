#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
======================================================================================
kincat version 1.0
Copyright (2023) NTESS
https://github.com/sandia/kincat
Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
This file is part of KinCat.
KinCat is open-source software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. You should have received a copy of the GNU Lesser General Public License along with KinCat. If not, see <http://www.gnu.org/licenses/>.

Questions? Contact Craig Daniels at <cjdanie@sandia.gov>

Sandia National Laboratories, Albuquerque, NM, USA
======================================================================================*/
N.B. This code is adapted from KineCluE (https://github.com/lukamessina/kineclue). 
Some unecessary portions may still remain, and some terminology or conventions may reflect that application.
These may include 'jump','event', or 'frequency' vs 'process' or 'process instance', and 'kinetic range (kira)'' vs 'range'
'''

import numpy as np
import sympy as sym
import time as tm
import kinepy as kp
import _pickle as pickle
import os
import sys
import logging
import datetime
import copy
from itertools import compress, product
from shutil import copyfile
from psutil import Process
import json
#import threading
import multiprocessing as mp
import functools
from tqdm.contrib.concurrent import thread_map
import subprocess

# Check for correct number of arguments
if len(sys.argv) != 2:
    raise Exception("ERROR! Incorrect call of {}. Correct usage: << ./kineclue_main input_file_path >>".format(sys.argv[0]))

sys.setrecursionlimit(kp.recursionlimit)  # for Pickle, else it causes error to save important data
start_time = tm.time() # Start measuring execution time
process_for_memory_tracking = Process(os.getpid()) # Save process id for memory usage tracking

# Get input file from first argument (#1)
myinput = sys.argv[1]

# Read input file
if os.path.exists(myinput):
    input_string = open(myinput, "r").read()
else:
    raise Exception ("ERROR! Input file {} not found.".format(myinput))

# Remove #-comments from input file (and remove new-line characters)
input_string = " ".join([string.split(sep='#')[0] for string in input_string.split(sep='\n') if len(string.split(sep='#')[0]) > 0])

# Definition of input_dataset dictionary
keywords = ['crystal', 'basis', 'directory', 'range', 'uniquepos', 'species',
            'procmech', 'fullsym']
input_dataset = {key: None for key in keywords}  # dictionary for reading user input
# Split input in keywords
input_list = input_string.split(sep='& ')
del input_list[0]  # delete all comments before first keyword
# Save input data into input_dataset
for ipt in input_list:
    keyword = ipt.split()[0].lower()  # get keyword
    input_dataset[keyword] = [i for i in ipt.split()[1:]]  # split all entries of each keyword in a list
del input_string, input_list, ipt

# Read output directory path from input
if input_dataset['directory'] is None:  # default directory is ./CALC/
    input_dataset['directory'] = ['./CALC/']
directory = input_dataset['directory'][0]
if directory == "cwd":
    directory = os.getcwd() + "/"
else:
    # Add a slash to the end of the directory path if it's missing
    if directory[-1] != "/":
        directory += "/"
# Save actual directory to input_dataset['directory'][0]
# the user-input string goes into input_dataset['directory'][1]
input_dataset['directory'].insert(0, directory)
# Create output directory
if not os.path.isdir(directory):
    os.makedirs(directory)

# Copying input file into directory
copy_input_path = directory + myinput.split("/")[-1]
if ".txt" not in myinput:
    copy_input_path += ".txt"
if not os.path.exists(copy_input_path):  # necessary when using the 'cwd' directory option
    copyfile(myinput, copy_input_path)

# Setting up logfile
logger = logging.getLogger('kinecluelog')
logger.setLevel(logging.INFO)
fh = logging.FileHandler(directory + 'kineclue_main.log', mode='w')
fh.setLevel(logging.INFO)
fh.setFormatter(logging.Formatter('%(message)s'))
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(fh)
logger.addHandler(ch)
logger.info(' __________________________________________________')
logger.info(' |                                                |')
logger.info(' |                   KinCatPy                     |')
logger.info(' |                 Craig Daniels                  |')
logger.info(' |                   Based on                     |')
logger.info(' |           KineCluE v{} - {}           |'.format(kp.version, kp.date))
logger.info(' |        T. Schuler, L. Messina, M. Nastar       |')
logger.info(' |________________________________________________|')
logger.info('')
kp.print_license_notice()
logger.info('')
logger.info('Calculation date: {}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
logger.info('Working in directory:  {}'.format(directory))

# Check that required information is in the input
if input_dataset['range'] is None:
    kp.produce_error_and_quit("!! You MUST provide a RANGE keyword.")
if input_dataset['species'] is None:
    kp.produce_error_and_quit("!! You MUST provide a SPECIES keyword.")
if input_dataset['procmech'] is None:
    kp.produce_error_and_quit("!! You MUST provide a PROCMECH keyword.")
if input_dataset['crystal'] is None:
    kp.produce_error_and_quit("!! You MUST provide a CRYSTAL keyword.")
if input_dataset['uniquepos'] is None:
    kp.produce_error_and_quit("!! You MUST provide a UNIQUEPOS keyword.")

# Convert kin range to float and round up to avoid numerical issues
KiRa = float(input_dataset['range'][0]) + 0.01
logger.info("Interaction range is set to {:.2f} lattice parameters".format(KiRa))

# Create crystal
crystal_input = input_dataset['crystal']  # !! Change of variable: crystal becomes an object
# SaveCheck for basis atoms (in a list of arrays)
if input_dataset['basis'] is None:
    logger.info("!! Bravais lattice because (BASIS keyword not found)")
    input_dataset['basis'] = ['s', '1', '0', '0', '0', '0']  # default crystal basis: Bravais lattice
n_basis_atoms =int(input_dataset["basis"][1])

# Check dimensionality of the crystal from input file
dim = int(np.sqrt(len(crystal_input)-1))
if n_basis_atoms > 1:
    logger.info('Creating {:.0f}D crystal {} with {:.0f} basis atoms'.format(dim, crystal_input[0], n_basis_atoms))
else:
    logger.info('Creating {:.0f}D crystal {}'.format(dim, crystal_input[0]))
if len(crystal_input) == 2: #1D crystal
    crystal_input += ['0.0', '0.0', '0.0', '6.5713', '0.0', '0.0', '0.0', '3.5187'] # adding small components in y and z directions
    kp.produce_error_and_quit("ERROR! KinCatPy is written for a 2D lattice/crystal")
elif len(crystal_input) == 5: #2D crystal, adding the third direction, perpendicular to the others
    crystal_input = crystal_input[0:3] + ['0.0'] + crystal_input[3:5] + ['0.0', '0.0', '0.0', '5.4618']
elif len(crystal_input) == 10: # 3D crystal, nothing to add
    kp.produce_error_and_quit("ERROR! KinCatPy is written for a 2D lattice/crystal")
    pass
else:
    kp.produce_error_and_quit("ERROR! Wrong definition of the crystal. Must be a string (name of the crystal), followed by 1, 4 or 9 numbers depending on the dimensionality of the system.")
# Create new crystal
crystal_input[1:10] = [kp.evalinput(i) for i in crystal_input[1:10]]
crystal = kp.Crystal(name=crystal_input[0], vec1=np.array(crystal_input[1:4]), vec2=np.array(crystal_input[4:7]),
                    vec3=np.array(crystal_input[7:10]), dim=dim)
# Set basis list and atomic volume
basis_list = []
for i in range(2, len(input_dataset["basis"])):
    input_dataset['basis'][i] = kp.evalinput(input_dataset['basis'][i])
for ba in range(n_basis_atoms):
    basis_vec = input_dataset["basis"][ba*3+2:ba*3+5]
    basis_vec.append(0.0)    
    basis_list.append(np.array(basis_vec, dtype=float))
    if input_dataset["basis"][0] == "o":  # convert position to supercell base
        basis_list[ba][1:4] = kp.change_base(arr=basis_list[ba][1:4], crystal=crystal)
crystal.set_basislist(basis_list)
crystal.set_atomicvolume()

# Find all valid symmetry operations
# crys_symop_list is a dictionary with following objects as keys: crystal, crystal_deformed
logger.info("Searching for symmetry operations in lattice...")
tmp_symop_list = kp.find_symmetry_operations(crystal=crystal)
json_dict={}
crystal_dict={}
edge_vectors=copy.copy(crystal.get_primvectors())
basis_vectors=copy.copy(crystal.get_basislist())
for i in range(len(basis_vectors)):
    basis_vectors[i]=basis_vectors[i].tolist()
    basis_vectors[i].pop(3) #remove z index of basis vector
    basis_vectors[i].pop(0) #remove index 0, which stores site type
crystal_dict.update({'edge vectors':[[edge_vectors[0][0],edge_vectors[1][0]],[edge_vectors[0][1],edge_vectors[1][1]]]})
crystal_dict.update({'basis vectors':basis_vectors})


# Saving (perfect) crystal object to file
pickle.dump([crystal, tmp_symop_list], open(directory + 'crystal_' + crystal.get_name() + '.pkl', 'wb'), -1)

symop_list = tmp_symop_list

del tmp_symop_list

if input_dataset['fullsym'] is None:
    full_sym_flag=False
else:
    full_sym_flag=True # Full symmetry output only if flagged in input file


crystal_file=directory+"KinCat_input.json"
short_symop_list=[]
if full_sym_flag: #Only pass identity symmetry operation
    short_symop_list.append(1)
    short_symop_list.append(0)
    short_symop_list.append(0)
    short_symop_list.append(1)
    short_symop_list.append(0.0)
    short_symop_list.append(0.0)
    n_output_symops=1
else: 
    for i in range(len(symop_list)):
        tmp=copy.copy(symop_list[i].get_rotation().tolist())
        tmp2=copy.copy(symop_list[i].get_translation()).tolist()
        tmp2.pop(2)
        tmp3=[tmp[0][0],tmp[0][1],tmp[1][0],tmp[1][1], tmp2[0], tmp2[1]]
        short_symop_list.append(tmp[0][0])
        short_symop_list.append(tmp[0][1])
        short_symop_list.append(tmp[1][0])
        short_symop_list.append(tmp[1][1])
        short_symop_list.append(tmp2[0])
        short_symop_list.append(tmp2[1])
    n_output_symops=len(symop_list)
sym_dict={}
sym_dict.update({"shape":[n_output_symops, 6]})
#json_dict.update({'number of symmetry operations':len(tmp_symop_list)})
sym_dict.update({'data':short_symop_list})
crystal_dict.update({"symmetry operations":sym_dict})
json_dict.update({"crystal":crystal_dict})

# Create list of defects from UNIQUEPOS
logger.info("Creating list of species and components, starting at {:.3f} s".format(tm.time()-start_time))
[defect_list, all_defect_list, doubledefects] = kp.dataset_to_defects(dataset=input_dataset, crystal=crystal, symop_list=symop_list)
n_defects = len(defect_list)
# defect_list contains all symmetry unique defects
# all_defect_list contains all defects
# doubledefects contains list of user-input defects that are symmetrically equivalent
# (they do not appear in defect_list/all_defect_list and need to be removed from the permission list in species - this is done in dataset_to_species)

# Create lists of species and components
[species_list, component_list] = kp.dataset_to_species(dataset=input_dataset, defect_list=defect_list, doubledefects=doubledefects)
n_species = len(species_list) - 1  # bulk is not included in the number of species
n_components = len(component_list)

# Creating list of component label permutations (for configuration names)
specs = []
for cp in component_list:
    specs.append(cp.get_species())  # list of species for each components
name_perms = kp.find_possible_permutations(specs=specs, species_list=species_list)
del cp, specs, doubledefects

# Create list of jumps (jump symmetries and constraints are analyzed in dataset_to_jumps)
# Jump_list contains complete list of jumps, including symmetry equivalents
logger.info("Creating list of processes, starting at {:.3f} s".format(tm.time()-start_time))
jump_list = kp.dataset_to_jumps_cat(dataset=input_dataset, crystal=crystal, symop_list=symop_list, all_defect_list=all_defect_list, species_list=species_list, sym_unique=True, component_list=component_list)
sym_jump_list, jump_sym_list = kp.dataset_to_jumps_cat(dataset=input_dataset, crystal=crystal, symop_list=symop_list, all_defect_list=all_defect_list, species_list=species_list, sym_unique=False, component_list=component_list)

if full_sym_flag:
    jump_list=sym_jump_list
    jump_sym_list=[[0] for _ in jump_list]

#print("\n Full Symmetries Jump List")
#print("length= ",len(sym_jump_list))
n_jumps = len(jump_list)  # total number of jump mechanisms (including symmetry equivalent ones)

#Create initial configuration 
#empty configuration assumes adsorbtion/desorption processes present
logger.info("Exploring Configuration Definition (starts at {:.3f} s)".format(tm.time()-start_time))

#Should this occur before finding symmetry equivalent jumps?
#Should there be any checking of equivalent coordinates by symmetry or translation?
    #Simpler to not
    #Checking symmetry can reduce number of configurations, so would be valuable...
event_coords = kp.find_event_coords(jump_list)
if not full_sym_flag: #If exploring reduced symmetry, still need to find interaction range of full symmetry for KinCat KMC. 
    full_event_coords = kp.find_event_coords(sym_jump_list)

logger.info("Found {:} process coordinates.".format(len(event_coords)))
#for cord in event_coords:
#    logger.info("{:}".format(cord))
logger.info("Exploring Configuration Relationships (starts at {:.3f} s)".format(tm.time()-start_time))

event_str_list=[]
for event in jump_list:
    event_str_list.append(str(event.get_name()))

#print("Finding configuration site_list:")
logger.info("Finding configuration sites.")

coords_list=kp.confsites(event_coords=event_coords, kira=KiRa, crystal=crystal, species=species_list) #finds all the sites that need to be specified
#logger.info("Initial guess at configuration sites: ")
#for cord in coords_list:
#    logger.info("{:}".format(cord))
if not full_sym_flag:
    full_coords_list=kp.confsites(event_coords=full_event_coords, kira=KiRa, crystal=crystal, species=species_list) #finds all sites needed to be specified in full symmetry case
    tmp=len(full_coords_list[0])
    x_max=0
    x_min=0
    y_max=0
    y_min=0
    for cord in full_coords_list:
        x_max=max(cord[0],x_max)
        x_min=min(cord[0],x_min)
        y_max=max(cord[1],y_max)
        y_min=min(cord[1],y_min)
    x_max=kp.index_round(x_max)
    x_min=kp.index_round(x_min)
    y_max=kp.index_round(y_max)
    y_min=kp.index_round(y_min)
    x_interaction=x_max-x_min
    y_interaction=y_max-y_min
else:
    tmp=len(coords_list[0])
    x_max=0
    x_min=0
    y_max=0
    y_min=0
    for cord in coords_list:
        x_max=max(cord[0],x_max)
        x_min=min(cord[0],x_min)
        y_max=max(cord[1],y_max)
        y_min=min(cord[1],y_min)
    x_max=kp.index_round(x_max)
    x_min=kp.index_round(x_min)
    y_max=kp.index_round(y_max)
    y_min=kp.index_round(y_min)
    x_interaction=x_max-x_min
    y_interaction=y_max-y_min
logger.info("Found {:} configuration coordinates.)".format(len(coords_list)))
#for cord in coords_list:
#    logger.info("{:}".format(cord))
if full_sym_flag:
        config_sym_relationships=[[cord for cord in range(len(coords_list))]]
        for p,op in enumerate(symop_list):
            if kp.are_equal_arrays(op.get_rotation(),np.identity(3)):
                temp_ind=copy.copy(p)
                config_symops=[temp_ind]
else:
    logger.info("Finding configuration symmetry relationships.")
    config_symops,config_sym_relationships= kp.config_sym_define(coords_list=coords_list, symop_list=symop_list)
#Translate event mechanisms into configuration index operations
config_event_list=kp.config_event_def(coords_list=coords_list, jump_list= jump_list) #may not need this anymore

n_coords=len(coords_list)
conf_symop_list=[]
for i in config_symops:
    conf_symop_list.append(kp.SymOp(symop_list[i].get_rotation()))
ini_defects=[]
for i in range(n_coords): #need to fix!!!???
    ini_defects.append([component_list[0].get_species().get_defects()[0].get_symeq()[0]]) #Assuming everything on main lattice (no basis or sublattice)

config_template=kp.CatConfTemplate(defects=ini_defects, translations=coords_list, config_symops=config_symops, symop_relations=config_sym_relationships, crystal=crystal, event_coords=event_coords)
config_template.set_events(config_event_list)
config_json_dict={}
json_coords=[]
for i in range(len(coords_list)):
    tmp=coords_list[i].tolist()
    tmp.pop(2)
    json_coords.append(tmp)
config_json_dict.update({"site coordinates":json_coords})
config_json_dict.update({"variant orderings":config_sym_relationships})
config_json_dict.update({"interaction range":[int(x_interaction+1),int(y_interaction+1)]})

# END READING INPUT AND FORMATTING DATA-------------------------

# Creating configuration space starting from iniconf and applying successive jumps to symmetry unique configurations
logger.info("Creating configuration space (starts at {:.3f} s)".format(tm.time() - start_time))
logger.info("Peak memory usage up to now: {:.3f} MB.".format(process_for_memory_tracking.memory_info()[0]*1e-6))  # peak memory usage

freq_list = {}  # list of jump frequencies (w)
config_list = {}  # comprehensive configuration list

##Configuration permutations:
def permute_species_configs(ix,perm) -> None:
    if ix<=n_coords:
        for i in range(n_species+1):
            new_perm=perm.copy()
            new_perm.append(i)
            permute_species_configs(ix+1,new_perm)
    else:
        ini_defects=[]
        ini_species=[]
        
        for i in range(n_coords):
            ini_defects.append([component_list[0].get_species().get_defects()[0].get_symeq()[0]]) #Assuming everything on main lattice (no basis or sublattice)
            ini_species.append(kp.Species(species_list[perm[i]].get_index(),species_list[perm[i]].get_name(),species_list[perm[i]].get_defects(),species_list[perm[i]].get_permissions(),species_list[perm[i]].get_radius()))
            #NB!!! [j.get_index() for j in ini_species] is NOT THE SAME as PERM! (0 index species is last in species list, not first). 
            #test_spec_count[species_list[perm[i]].get_index()]+=1 #Do I need this anymore?

        #Check if configurations match system definition
        if (kp.check_subconfiguration_consistence(species_list=ini_species, position_list=coords_list, all_defect_list=all_defect_list,crystal=crystal,component_list=[])):
            if (kp.check_configuration_exclusion(species_list=ini_species, config_template=config_template)):
                ini_species_ind=[j.get_index() for j in ini_species]
                kp.sym_config_add_unique(ini_species_ind=ini_species_ind, config_template=config_template, config_dict=config_dict, config_species_lists=config_species_lists) # adds configuration if not already found


config_dict={}
config_species_lists=[]
permute_species_configs(ix=1,perm=[])

config_json_dict.update({"shape": [len(config_species_lists), len(config_species_lists[0])]})
config_data=np.array(config_species_lists)
config_data=np.reshape(config_data,(len(config_species_lists)*len(config_species_lists[0])))
config_json_dict.update({"data": config_data.tolist()})
json_dict.update({"configurations":config_json_dict})
logger.info("Found {:} configurations.".format(len(config_species_lists)))
logger.info("Peak memory usage up to now: {:.3f} MB.".format(process_for_memory_tracking.memory_info()[0]*1e-6))  


logger.info("Finding process instances (starts at {:.3f} s)".format(tm.time() - start_time))

freq_list=[]
events_set= set()
event_list=[]
jump_list_str=[jump.get_name() for jump in jump_list]
event_json_dict={}
event_json_dict.update({"processes": jump_list_str})
event_json_dict.update({"process constraints":config_event_list})
event_json_dict.update({"process symmetries": jump_sym_list})



#Begin Parallel Block
#Uses threading/multiprocessor to try to speed up finding events (process instances)
full_event_list=[]
full_event_set=set()
symmetry_orders=config_template.get_relations()
coord_dict=config_template.get_coord_dict()
non_sym_event_list=[]


if (len(config_species_lists)<100000): #Serial is fast enough, little to no benefit to parallelizing due to overhead
    full_event_list=kp.find_all_events_config(allevents=config_event_list,config_species_lists=config_species_lists, config_template=config_template, config_dict=config_dict)

else:
    ## Multiprocessing Attempt
    if __name__=="__main__":
        #print("In multi")

        ##Using separate script, single output file. Works, but only achieves speedup of 3x, benefits maximxize with 8 processors with 8-core machine
        nproc=8
        logger.info("Conducting parallel process instance search.")
        logger.info("Using {} processors.".format(nproc))
        pickle.dump([config_species_lists, config_template, config_dict, nproc],open(directory + 'events_calc.pkl', 'wb'), -1)
        bashCommand=('python ../kincat_event_calc.py '+(directory))
        #print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout= subprocess.PIPE)
        output, error = process.communicate()
        #print(output)
        full_event_list=pickle.load(open((directory+'events_list.pkl'), 'rb'))

#End Parallel Block

event_json_dict.update({"shape": [len(full_event_list), 3]})
event_data=np.array(full_event_list)

event_data=np.reshape(event_data,(len(full_event_list)*3)).tolist()
event_json_dict.update({"data": event_data})
json_dict.update({"process dictionary": event_json_dict})
#print(len(full_event_list)," event frequencies found.")

logger.info("Found {:} process instances.".format(len(full_event_list)))

# Writing configuration file 
logger.info("Writing configuration values file (starts at {:.3f} s)".format(tm.time() - start_time))
conf_file=directory+"configurations.json"


with open(crystal_file,"w") as write_file:
    json.dump(json_dict,write_file)#, indent=2)
conf_file = directory + "configurations.txt"
disfmt="{:^"+str(2*(len(species_list)-1)*len(component_list)+2)+"s}"
with open(conf_file, 'w') as output:
    output.writelines("Configuration Site Positions (supercell coordinates)\n")
    for i in coords_list:
        output.writelines(str(i))
    output.writelines("\n1) Configuration class; 2) Entropy prefactor (no units); 3) Binding energy (eV, >0 means attraction);"+
                      "4) Species List\n")
    for i,conf in enumerate(config_species_lists):
        tofile = []
        tofile += ["{:<5.0f}".format(i)]#conf.get_label())] #configuration thermodynamic interaction
        tofile += ["{:5.3f}".format(-1)]  # default prefactor
        tofile += ["{:4.3f}".format(0)] #default binding energy
        output.writelines("{:s}  {:s}  {:s} ".format(*tofile))
        tofile=[]
        tofile += conf
        output.writelines(["%s" % item + ' ' for item in tofile])
        output.writelines("%s\n" % "")


# Writing jump frequency file (ThRa)
logger.info("Writing process instances file (starts at {:.3f} s)".format(tm.time() - start_time))
freq_file = directory + 'process_instances.txt'

with open(freq_file, 'w') as output:
    output.writelines("1) Jump frequency number; 2) Jump prefactor (no units) 3) Saddle-point energy (eV);"+
        " 4) Initial Configuration class; 5) Final Configuration class; 6) Jump mechanism name;\n")
    for i,freq in enumerate(full_event_list):
        tofile = []
        tofile += ["{:<5.0f}".format(i)] #frequency label
        tofile += ["{:5.3f}".format(-1)] #default prefactor
        tofile += ["{:4.3f}".format(0)] #default binding energy
        tofile += ["{:^5.0f}".format(freq[0])] #initial configuration label
        tofile += ["{:^5.0f}".format(freq[1])] #final configuration label
        tofile += ["{}".format(freq[2])] #event name
        output.writelines("{:s}  {:s}  {:s} {:s} {:s} {:s}\n".format(*tofile))

# Stop measuring execution time and print elapsed time
stop_time = tm.time()
logger.info("Execution time: {:.3f} s.".format(stop_time - start_time))
logger.info("Peak memory usage: {:.3f} MB.".format(process_for_memory_tracking.memory_info()[0]*1e-6))  # peak memory usage


# END CODE-----------------------------------
