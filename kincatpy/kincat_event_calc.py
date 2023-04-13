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
This script is intended to allow only the findevents_proc function to be run on multiple processors.
It is called from kincatpy.py .  
The resulting events are combined and loaded back into kincatpy.py
'''

import numpy as np
import _pickle as pickle
import subprocess
import sys
import multiprocessing as mp
import copy


#This method utilizes a single input file for all configurations, and uses a configuration-based description of mechanisms to find events.
directory= sys.argv[1]
myinput = (directory +'events_calc.pkl')
#print(myinput)
[config_species_lists, config_template, config_dict, nproc]=pickle.load(open(myinput, 'rb'))
#print("IN EVENT CALC")
#print(config_dict)
#divide configurations into nproc sections
config_per_proc=int(np.ceil(float(len(config_species_lists))/nproc))
#print(" config_per_proc= ",config_per_proc,' ')

def findevents_proc(proc) -> list:
	#print(" In findevents_proc , proc= ",proc)
	#Checks list of configuration species lists. Finds all events and returns list of event information formatted: [initial, final, event_type_index].
	full_event_list=[]
	full_event_set=set()
	symmetry_orders=config_template.get_relations()
	coord_dict=config_template.get_coord_dict()
	sym_species_list=[config_species_lists[0] for _ in range(len(symmetry_orders))] # dummy values, overwrite later
	jump_list=config_template.get_events()
	#Only find events for configurations in one section of list
	config_start=proc*(config_per_proc)
	config_end= (proc+1)*config_per_proc
	if config_end>len(config_species_lists):
		config_end=len(config_species_lists)
	#print(" Config start/end= ",config_start, ' ', config_end,' ')
	for c in range(config_start,config_end):
	#print(" c= ",c)
		#print(config_species_lists[c])
		for i in range(len(symmetry_orders)):
			sym_species=[config_species_lists[c][j] for j in symmetry_orders[i]]
			sym_species_list[i]=sym_species
		for sym in sym_species_list:
			for e, event in enumerate(jump_list): #loop over each event
				found_event_flag=True
				finspec=copy.copy(sym)
				for con in event:
					if sym[con[0]]==con[1]:
						finspec[con[0]]=con[2]
					else:
						found_event_flag=False 
						break
				if found_event_flag: #All constraint initial conditions met.
					finconfig_ind = -1
					for i in range(len(symmetry_orders)):
						sym_finspec=[finspec[j] for j in symmetry_orders[i]]
						finspec_str=np.array_str(np.array(sym_finspec))
						if finspec_str in config_dict:
							finconfig_ind= int(config_dict[finspec_str])
							break
					if (finconfig_ind == -1):
						#Constraints for initial configuration met, but no matching final configuration due to size exclusions. Event disallowed.
						continue
					unique_event=[c, finconfig_ind, e]
					unique_event_str=np.array_str(np.array(unique_event))
					if not unique_event_str in full_event_set:
						full_event_set.add(unique_event_str)
						full_event_list.append(unique_event)
	return full_event_list

print(' NProc= ',nproc,' ')
complete_event_list=[]
if __name__ == '__main__':
	with mp.Pool() as parpool:
		results=[]
		for proc in range(nproc):
			results.append( parpool.apply_async( findevents_proc, (proc,)))
		#print("Par. Ran")
		for res in results:
			answer = res.get()
			for event in answer:
				complete_event_list.append(event)
#print("Number of Events found = ",len(complete_event_list))
#for event in complete_event_list:
#	print(' ',event,' ')
pickle.dump(complete_event_list, open(directory + 'events_list.pkl', 'wb'), -1)
