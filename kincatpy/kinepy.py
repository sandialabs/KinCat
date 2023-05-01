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
import copy
import os
import _pickle as pickle
import mpmath as mp
import logging
import datetime
from scipy.sparse import linalg
from itertools import product, permutations, combinations
from shutil import copyfile


# Definition of variables
version = "1.0"
date = "01/04/2023"
tol = 1e-12  # numeric tolerance on the comparisons (machine error)
tol_strain = 1e-6
db_displ = 0.0003  #  fictitious displacement of each atom in a dumbbell (in lattice parameter units, with respect to dumbbell center)
maxdbdist = 0 # CHECK FOR POSSIBLE ROUNDING ERRORS FOR DUMBBELLS! HARD CODED '4' valid for 100,110,111 dumbbells.
              # The correct value should be sqrt(or.or) where or is the dumbbells orientation vector
              # Other possibility is to normalize dumbbells d distances to always have the same dumbbell bond length
tol_pbc = 0.001  # tolerance when applying periodic boundary conditions (to exclude dumbbell from automatic translations)
recursionlimit = 10000  # for Pickle, else it causes error to save important data
e = sym.Symbol('e')  # symbolic variable for strain
strainval = 0.00435427869  # numerical value for strain if needed
# Special variables to treat bulkspecies (symmetric species that lead to a double counting of the symmetry equivalent configurations - e.g. pure dumbbells)
bulkspecies_list = {}  # list of special "bulkspecies" (list of indices corresponding to the bulkspecies components)
bulkspecies_symdefects = {}  # dictionary where for each key (each defect in all_defect_list) corresponds the defect with the opposite sublattice, so that the first component is positive

# Definition of logger
logger = logging.getLogger('kincatpylog')

# Definition of classes
class Crystal:
    """The crystal object contains the primitive vectors (as numpy arrays) and the methods to perform a base change"""

    def __init__(self, name, vec1, vec2, vec3, dim):
        self.__name = name
        self.__primvectors = np.transpose(np.array([vec1, vec2, vec3],dtype='float'))   # Primitive vectors (vec1, vec2, vec3 in columns)
        self.__invvectors = np.linalg.inv(self.__primvectors)  # Compute inverse matrix of prim. vectors for base changes
        self.__strain = np.zeros((3,3), dtype=float)  # strain in the cell
        self.__orthostrain = np.zeros((3,3), dtype=float)  # strain in orthonormal coordinates
        self.__atomicvolume = 0
        self.__basislist = []
        self.__dimensionality = dim # 1D, 2D or 3D
        self.__primvectorsnum = self.get_primvectors()  # with symbolic strain, primvectors become object array, not float array
        self.__invvectorsnum = self.get_invvectors()  # with symbolic strain, invvectors become object array, not float array

    def get_name(self):
        return self.__name

    def get_primvectors(self):
        return self.__primvectors

    def get_primvectorsnum(self):
        return self.__primvectorsnum

    def get_invvectors(self):
        return self.__invvectors

    def get_invvectorsnum(self):
        return self.__invvectorsnum

    def get_strain(self):
        return self.__strain

    def get_orthostrain(self):
        return self.__orthostrain

    def get_basislist(self):
        return self.__basislist

    def get_atomicvolume(self):
        return self.__atomicvolume

    def get_dim(self):
        return self.__dimensionality

    def set_basislist(self, a):
        self.__basislist = a

    def set_atomicvolume(self):
        # self.__atomicvolume = np.abs(np.linalg.det(self.__primvectors))/len(self.__basislist)
        self.__atomicvolume = np.abs(sym.simplify(sym.Matrix(self.__primvectors).det())) / len(self.__basislist)

    def set_strain(self, strain: np.ndarray):  # Apply user-input strain
        self.__strain = strain
        self.__orthostrain = change_base(arr=strain, crystal=self, inv=True, matrix=True)
        self.__primvectors = np.dot(self.__primvectors, (np.identity(3) + self.__strain))   # Check this formula!
        self.__primvectorsnum = sym.matrix2numpy(sym.Matrix(self.__primvectors).subs(e, strainval), dtype='float')
        #self.__invvectors = np.linalg.inv(self.__primvectors)  # Could be written without computing linalg.inv
        self.__invvectors = sym.simplify(sym.Matrix(self.__primvectors).inv())
        self.__invvectorsnum = sym.matrix2numpy(self.__invvectors.subs(e, strainval), dtype='float')
        self.__invvectors = sym.matrix2numpy(self.__invvectors, dtype='object')


class Component:
    """A component belongs to a species which defines the possible defects (sublattice/orientation) available"""

    def __init__(self, species, index):
        self.__species = species   # species that the component belong to
        self.__index = index  # component index

    def get_species(self):
        return self.__species

    def get_sublattices(self):
        return tuple([o.get_sublattice() for o in self.__species.get_defects()])

    def get_defects(self):
        return flat_list([o.get_symeq() for o in self.__species.get_defects()])

    def get_info(self):
        logger.info("    Component {} is a {}". format(self.__index, self.__species.get_name()))

class CatConfTemplate:
    """Since all configurations have the same generic coordinates, this establishes the defects and translations lists once, 
    so that doesn't have to be stored with each CatConf.
    Also stores relationships between sites and symmetry operations, for easier checking of symmetry equivalence."""
    def __init__(self,defects, translations, config_symops, symop_relations, crystal, event_coords):
        self.__defects = defects #List of defects (a defect for each site)
        self.__translations = translations #Indices of the supercell where each species is sitting, besides the first one, which is always in [0, 0, 0].
        self.__config_symops = config_symops
        self.__symop_relations = symop_relations
        self.__n_sites = len(translations)
        self.__event_site_indicies = []

        for i in event_coords:
            for j in range(self.__n_sites):
                if are_equal_arrays(np.array(i),np.array(translations[j])):
                    self.__event_site_indicies.append(j)

        positions = []
        self.__distances = [] #distances between sites, orthonormal, not supercell coordinates
        for i in translations:
            positions.append((change_base(i.tolist(),crystal, inv=True, matrix=False)).tolist())
        for i in range(self.__n_sites):
            dist_list=[]
            for j in range(self.__n_sites):
                temp=np.array(positions[i])-np.array(positions[j])
                temp2=np.sqrt((temp[0]**2)+(temp[1]**2)+(temp[2]**2))
                dist_list.append(temp2)
            self.__distances.append(dist_list)

        self.__coord_dict={}
        for i in range(self.__n_sites):
            self.__coord_dict.update({np.array_str(np.array(translations[i])):i})

    def get_relations(self):
        return self.__symop_relations

    def get_translations(self):
        return self.__translations

    def get_defects(self):
        return self.__defects 

    def get_n_sites(self):
        return self.__n_sites

    def get_distances(self): 
        return self.__distances

    def get_event_sites(self):
        return self.__event_site_indicies

    def get_coord_dict(self):
        return self.__coord_dict

    def set_events(self, config_event_list:list):
        self.__event_list = config_event_list

    def get_events(self):
        return self.__event_list


class Configuration:
    """A Configuration contains defect types and position for each component"""

    def __init__(self, defects, translations, label=None):
        self.__defects = defects  # List of defects (a defect for each component)
        self.__translations = translations   # Indices of the supercell where each component is sitting, besides the first one, which is always in [0, 0, 0].
        # Position is given by self.__translations * crystal_primvectors + self__defects
        self.__th0 = None  # Thermodynamic interaction
        self.__nu1 = None  # Kinetic interaction
        self.__jumplist = []  # List of list of accessible configurations and corresponding jump mechanism
        self.__label = label # name of the configuration in c[defects_][translations] or d[defects_][translations] if kinetically dissociated
        self.__subconfigurations = [] # list of subconfigurations to apply thermodynamic range reduction
        self.__species = [] # list of species, needed to compare subconfigurations
        self.__beyond = True  # true if this configuration is thermodynamically but not kinetically dissociated
        self.__kineticinter = [] # list of kinetic interaction indexes associated with this configuration
        self.__posthreplist = None
        self.__subname = None


    def get_defects(self):
        return self.__defects

    def get_translations(self):
        return self.__translations

    def get_defect_position(self, def_idx: int) -> np.ndarray:
        """Return defect position as sublattice+translation (defect order is that defined in list self.__defects)"""
        return np.array(self.get_defects()[def_idx].get_sublattice(), dtype=float) + np.array(self.get_translations()[3*def_idx:3*def_idx+3], dtype=float)

    def get_thermoint(self):
        return self.__th0

    def get_kinint1(self):
        return self.__nu1

    def get_jumplist(self):
        return self.__jumplist

    def get_label(self):
        return self.__label

    def get_kineticinter(self):
        return self.__kineticinter

    def get_posthreplist(self):
        return self.__posthreplist

    def set_posthreplist(self, a):
        self.__posthreplist = a

    def get_subname(self):
        return self.__subname

    def set_translations(self, a):
        self.__translations = a

    def set_kineticinter(self, a):
        if isinstance(a, list):
            self.__kineticinter.extend(a)
        else:
            self.__kineticinter.append(a)

    def set_thermoint(self, a):
        self.__th0 = a

    def set_kinint1(self, a):
        self.__nu1 = a

    def set_label(self, a):
        self.__label = a

    def set_jumpfrequency(self, pos, freq):
        self.__jumplist[pos][4]=freq

    def get_beyond(self):
        return self.__beyond

    def set_beyond(self):
        self.__beyond = False


    def cleanconfig(self):
        #self.__beyond = None
        self.__nu1 = None
        #self.__subconfigurations = None
        self.__species = None
        self.__jumplist = None

    def get_subconfigurations(self):
        return self.__subconfigurations

    def get_species(self):
        return self.__species

    def set_subconfiguration(self, a):
        self.__subconfigurations.append(a)

    def set_species(self, a):
        self.__species += a

    def replace_jumplist(self, old: str, new: str):
        for i, k in enumerate(self.__jumplist):
            if k[0] == old:
                self.__jumplist[i][0] = new

    def replace_disso_jumplist(self, toreplace: list):
        n = -1
        for i, k in enumerate(self.__jumplist):
            if k[0][0] == 'd': #disso jump; remove it and replace by another name from toreplace list
                n += 1
                self.__jumplist[i][0] = toreplace[n]
    

class Defect:
    """A defect object is defined by a sub-lattice"""
    def __init__(self, crystal, sublattice, index, sublattice_nod):
        self.__crystal = crystal  # Crystal object in which the defect is defined
        self.__sublattice = sublattice  # position of the defect in the primitive cell
        self.__sublattice_nod = sublattice_nod # same but removing the dumbbell displacement if any
        self.__index = index  # index labeling each defect
        # List of symmetry equivalent defect objects in the supercell (e.g. there are 12 for a mixed dumbbell)
        # This is a list of defect objects (object in the object!) where, for each of them, the list __symeqdefects is empty,
        # sublattice and direction are obtained with method find_symeq(symop), and oriented and index are the same as the mother class.
        self.__symeqdefects = []
        # symindexlist is a list of lists showing the resulting defect from a given symmetry operation
        # first item in the sublist is the index of the defect resulting from the symmetry operation
        # next three items in the sublists are translations to account for translations of defect positions
        # because defects are always defined inside the base supercell (coordinates between 0 and 1)
        # For instance an inversion on +0.5,0,0 gives -0.5,0,0 which is in fact described by +0.5,0,0
        # and a translation of -1,0,0 supercells such that the sublist would be [a,-1,0,0] if this defect is indexed a
        self.__symindexlist = []
        self.__sudefect = None # symmetry unique defect from which the current defect was obtained by symmetry
        self.__align_dist = None # distance between two "aligned" defects to make sure that dumbbell configurations are correct
        self.__siteinter = None # dictionary relating species at this defect position to site interactions

    def get_crystal(self):
        return self.__crystal

    def get_sublattice(self):
        return self.__sublattice

    def get_sublattice_nod(self):
        return self.__sublattice_nod

    def get_index(self):
        return self.__index

    def get_symeq(self):
        return self.__symeqdefects

    def get_symindexlist(self):
        return self.__symindexlist

    def append_symindexlist(self, indexitem):
        self.__symindexlist.append(indexitem)

    def get_info(self):
        logger.info("  Defect {} on sublattice {} has {} symmetry equivalent.".format(self.__index, self.__sublattice, len(self.__symeqdefects)))

    def get_sudefect(self):
        return self.__sudefect

    def set_sudefect(self, a):
        self.__sudefect = a

    def set_align_dist(self): # to check distance between atmoms sharing a site if it is a dumbbell
        self.__align_dist = np.max([distance(vect1=self.__symeqdefects[a].get_sublattice()-self.__symeqdefects[a].get_sublattice_nod(),
                                             vect2=self.__symeqdefects[b].get_sublattice()-self.__symeqdefects[b].get_sublattice_nod(),
                                             crystal=self.__crystal) for a in range(len(self.__symeqdefects)) for b in range(a,len(self.__symeqdefects))])

    def get_align_dist(self):
        return self.__align_dist

    # Find list of symmetry equivalent defects, based on list of symmetry operations symop.
    def find_symeq(self, symop):
        # Apply symmetry operations to sublattice (defect position in the primitive cell) and direction (orientation)
        symeq = apply_symmetry_operations([self.__sublattice, self.__sublattice_nod], symop)
        # Apply "periodic boundary conditions"
        for eq in symeq:
            for sub in range(2): # with and without dumbbell displacement
                eq[sub] = apply_pbc(np.array(eq[sub], dtype=float))
            # Check that the new found defect was not already in the list.
            found = False
            for r in range(0, len(self.__symeqdefects)):
                refdefect = np.array(self.get_symeq()[r].get_sublattice(), dtype=float)
                if are_equal_arrays(refdefect, eq[0]):
                    found = True
                    break
            if not found:
                self.__symeqdefects.append(Defect(crystal=self.__crystal, sublattice=eq[0], sublattice_nod=eq[1], index=self.__index+len(self.__symeqdefects)))
        # Find the defect indices for each symmetry operation to speed up the identification of symmetry equivalent configurations
        for test_defect in self.__symeqdefects:
            test_defect.set_sudefect(self) # set symmetry unique defect
            for operator in symop:
                eq = apply_symmetry_operations(vector_list=[test_defect.get_sublattice()], symop_list=[operator], unique=False, applytrans=True)[0][0]
                trans = [0., 0., 0.]
                for r in range(3):
                    if eq[r] < (0 - tol_pbc):
                        trans[r] = -np.ceil(-eq[r])  # we need to get this translation correct that is why we are not using apply_pbc function
                        eq[r] += np.ceil(-eq[r])
                    elif eq[r] >= (1 - tol_pbc):
                        trans[r] = +np.floor(eq[r]+tol_pbc)  # we need to get this translation correct that is why we are not using apply_pbc function
                        eq[r] += -np.floor(eq[r]+tol_pbc)
                for r in range(0, len(self.__symeqdefects)):
                    refdefect = np.array(self.get_symeq()[r].get_sublattice(), dtype=float)
                    if are_equal_arrays(refdefect, np.array(eq, dtype=float)):
                        test_defect.append_symindexlist(flat_list([[self.get_symeq()[r].get_index()],trans]))
                        break
            if len(test_defect.get_symindexlist()) != len(symop):
                produce_error_and_quit("In Defect.find_symeq - Some defect symmetry equivalent indices of defect {} were not found.".format(self.__get_sublattice()))

    #def set_siteinter(self, species_list, species, siteinter):
    #    if self.__siteinter is None:
    #        self.__siteinter={x: 0 for x in species_list} # site inter
    #    self.__siteinter[species] = siteinter # interaction 0 is for bulk species or forbidden defects or a given species

    def get_siteinter(self): # requires a species object, not bulk species
        return self.__siteinter

    def set_siteinter(self, a):
        self.__siteinter = a

    def find_siteinter_classes(self, species_list: list, index: int, symopcpg: list, all_defect_list: list): # must be applied to a su_defect
        tmp = copy.copy(self.__symeqdefects) # list containing all symmetry equivalent defects of this Symmetry-Unique defect
        while len(tmp) > 0: # the first element of tmp is the defect that is currently being looked into
            if tmp[0].get_siteinter() is None:
                found = False
                for symeq,symop in zip(tmp[0].get_symindexlist()[:len(symopcpg)], symopcpg):
                    if symeq[0] == tmp[0].get_index() and symop.get_cpgsign() == -1: # there is a symmetry operation that reverses the CPG while preserving this defect. Site interaction is not necessary
                        found =True
                        break
                if found: # site interaction is not needed
                    for a in set([x[0] for x in tmp[0].get_symindexlist()[:len(symopcpg)]]):
                        all_defect_list[a].set_siteinter(a={x: 0 for x in species_list})  # site inter set to zero for all species because no site interaction is required
                        if all_defect_list[a] in tmp:
                            tmp.remove(all_defect_list[a]) # this defect has already been assigned a site interaction so it does not need to be studied further on
                else: # site interaction is needed
                    siteinterdict = {x: 0 for x in species_list}
                    minussiteinterdict = {x: 0 for x in species_list} # dict for symmetry equivalents reversing the CPG direction
                    for sp in species_list: # checking if this species can occupy this type of defect
                        if self in sp.get_defects(): # species sp can be in this type of defect (non zero permission)
                            index += 1 # new interaction
                            siteinterdict[sp] = index
                            minussiteinterdict[sp] = -index
                    # assigning the same siteinterdict dictionary to all symmetry equivalents
                    for symeq in set([(x[0],y.get_cpgsign()) for x,y in zip(tmp[0].get_symindexlist()[:len(symopcpg)], symopcpg)]):
                        if all_defect_list[symeq[0]].get_siteinter() is not None:
                            produce_error_and_quit("In find_siteinter_classes - defect has already been assigned a site interaction...this is odd.")
                        if symeq[1] == 1:
                            all_defect_list[symeq[0]].set_siteinter(a=siteinterdict)
                        elif symeq[1] == -1:
                            all_defect_list[symeq[0]].set_siteinter(a=minussiteinterdict)
                        else:
                            produce_error_and_quit("In find_siteinter_classes - CPG conserving symmetry operation is neither +1 or -1...this is odd")
                        if all_defect_list[symeq[0]] in tmp:
                            tmp.remove(all_defect_list[symeq[0]])
            else:
                produce_error_and_quit("In find_siteinter_classes - defect has already been assigned a site interaction...this is odd.")
        return index


class JumpConstraint:
    """A JumpConstraint defines the initial position and orientation of a given species
    for a jump mechanism to occur. Note that a jump will only be possible if all constraints
    are verified"""

    def __init__(self, inispecies, iniposition, iniposition_nod, finspecies, finposition, finposition_nod):
        # String label (for comparison purposes)
        self.__label = str(inispecies.get_index()) + "_" + ' '.join(['{:.6f}'.format(np.round(a+1e-8, 6)) for a in iniposition]) +\
                    "_" + str(finspecies.get_index()) + "_" + ' '.join(['{:.6f}'.format(np.round(a+1e-8, 6)) for a in finposition])
        # "_nod" positions are juste the same as position but the eventual dumbbell displacement is set 0 to have the exact jump vector.
        # Attributes before the jump
        self.__inispecies = inispecies  # species type
        self.__iniposition = np.array(iniposition)  # position of the defect
        self.__iniposition_nod = np.array(iniposition_nod)  # position of the defect without the dumbbell displacement from site
        self.__inidefect = None  # defect type
        self.__initrans = None  # integer supercell translation
        # Attributes after the jump
        self.__finspecies = finspecies
        self.__finposition = np.array(finposition)
        self.__finposition_nod = np.array(finposition_nod)
        self.__findefect = None
        self.__fintrans = None

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            complist = []
            complist.append(self.__inispecies == other.get_inispecies())
            complist.append(self.__finspecies == other.get_finspecies())
            complist.append(are_equal_arrays(self.__iniposition, other.get_iniposition()))
            complist.append(are_equal_arrays(self.__finposition, other.get_finposition()))
            return np.all(complist)
        else:
            return False

    def get_inispecies(self):
        return self.__inispecies

    def get_iniposition(self):
        return self.__iniposition

    def get_iniposition_nod(self):
        return self.__iniposition_nod

    def get_finspecies(self):
        return self.__finspecies

    def get_finposition(self):
        return self.__finposition

    def get_finposition_nod(self):
        return self.__finposition_nod

    def get_inidefect(self):
        return self.__inidefect

    def get_initrans(self):
        return self.__initrans

    def get_findefect(self):
        return self.__findefect

    def get_fintrans(self):
        return self.__fintrans

    def get_label(self):
        return self.__label

    def get_info(self):
        return "From species {} as defect {}  in position {} to species {} as defect {}  in position {}"\
            .format(str(self.get_inispecies()), str(self.get_inidefect()), str(self.get_iniposition()),\
                    str(self.get_finspecies()), str(self.get_findefect()), str(self.get_finposition()),)

    def vect2deftrans_constraint(self, all_defect_list):
        """Assign defect and translation attributes to the initial and final configuration of the constraint"""
        # (third item in mylist is a defect object)
        mylist = vect2deftrans(vec=copy.copy(self.__iniposition), all_defect_list=all_defect_list)
        self.__inidefect = mylist[0]
        self.__initrans = mylist[1]
        mylist = vect2deftrans(vec=copy.copy(self.__finposition), all_defect_list=all_defect_list)
        self.__findefect = mylist[0]
        self.__fintrans = mylist[1]

    def check_sublattices(self, all_defect_list: list) -> bool:
        """Check that all species in the constraint are sitting in the correct sublattices"""
        self.vect2deftrans_constraint(all_defect_list=all_defect_list)
        check_ini = check_subconfiguration_consistence(species_list=[self.__inispecies], position_list=[self.__iniposition], all_defect_list=all_defect_list)
        check_fin = check_subconfiguration_consistence(species_list=[self.__finspecies], position_list=[self.__finposition], all_defect_list=all_defect_list)
        return (check_ini and check_fin)


class JumpMech:
    """"A Jump is read from the input and contains a number of constraints"""

    def __init__(self, index, name):
        self.__index = index   # index identifying the jump mechanism
        self.__label = name  # label containing jump name and constraint labels (for comparison purposes)
        self.__nconstraints = 0   # number of constraints in the jump mechanism
        self.__name = name  # name of this jump mechanism
        self.__constraints = []  # list of constraint objects
        self.__netdisp = []  # list of net displacements for each species ordered by species index
        self.__netdisp_cpg = []  # list of net displacements for each species by species index projected on cpg and perpendicular directions
        self.__symeqs = [] # list of symmetric equivalents of this jump
        self.__summednetdisp = None  # average of net displacements over all symmetry equivalent
        self.__cpgsummednetdisp = None  # average over symmetry equivalents of net displacements projected on cpg and perp. directions

    def __eq__(self, other):
        #if isinstance(other, self.__class__)  and  self.get_name() == other.get_name():
        if isinstance(other, self.__class__):
            complist = [constr1 == constr2 for constr1, constr2 in zip(self.__constraints, other.get_constraints())]
            return np.all(complist)
        else:
            return False

    def get_nconstraints(self):
        return self.__nconstraints

    def get_name(self):
        return self.__name

    def get_constraints(self):
        return self.__constraints

    def get_index(self):
        return self.__index

    def get_label(self):
        return self.__label

    def get_netdisp(self):
        return self.__netdisp

    def get_netdisp_cpg(self):
        return self.__netdisp_cpg

    def get_symeqs(self):
        return self.__symeqs

    def get_summednetdisp(self):
        return self.__summednetdisp

    def get_cpgsummednetdisp(self):
        return self.__cpgsummednetdisp

    def add_constraint(self, c):  # c is a JumpConstraint object
        self.__constraints.append(c)
        self.__nconstraints = len(self.__constraints)
        self.__label += "_" + c.get_label()

    def set_symeqs(self, symops: list, jump_catalogue: list, all_defect_list: list):
        self.__symeqs = self.find_symeqs(symops=symops, old_jumps=jump_catalogue, unique=False, addreverse=False, all_defect_list=all_defect_list)

    def set_name(self, new_name: str):
        self.__name = new_name

    def set_summednetdisp(self, a):
        if self.__summednetdisp is None:
            self.__summednetdisp = a

    def set_cpgsummednetdisp(self, a):
        if self.__cpgsummednetdisp is None:
            self.__cpgsummednetdisp = a

    def set_net_species_displacement(self, strcrystal: Crystal, directions: list, spec: list):
        species_disp = [np.array([0, 0, 0], dtype=float) for _ in range(0, len(spec))] # list of displacement vectors for each species
        #comp_disp = [np.array([0, 0, 0], dtype=float) for _ in range(0, len(spec))] # list of displacement vectors for each component
        for cons in self.__constraints:
            species_disp[cons.get_inispecies().get_index()] = species_disp[cons.get_inispecies().get_index()] - np.array(cons.get_iniposition_nod(), dtype=float)
            species_disp[cons.get_finspecies().get_index()] = species_disp[cons.get_finspecies().get_index()] + np.array(cons.get_finposition_nod(), dtype=float)
        self.__netdisp = species_disp
        for species in species_disp:
            # Projection of species displacements on three specific diffusion directions
            self.__netdisp_cpg.append(np.array(proj(vect=species, crystal=strcrystal, direc=directions), dtype=object))

    
    def find_symeqs(self, symops: list, all_defect_list: list, old_jumps: dict = {}, unique: bool = True, addreverse: bool = True, print: bool = True, symop_ind: bool = False) -> list:
        """Applies a set (symops) of symmetry operations to a jump mechanism to find all symmetry equivalent jumps"""
        jump_list = []  # list of jump mechanism objects
        vectors = []  # list of vectors containing initial positions for each constraint, and then final positions for each constraint
        reversejump = []  # to make sure reverse jumps are included we also try to reverse initial and final positions
        symop_indices_list= []
        # Build vectors containing all initial and final positions for each constraint
        for jcons in self.__constraints:
            vectors.append(jcons.get_iniposition())
            reversejump.append(jcons.get_finposition())
        for jcons in self.__constraints:
            vectors.append(jcons.get_finposition())
            reversejump.append(jcons.get_iniposition())
        # Apply symmetry operations on initial and final positions for each constraint simultaneously
        tmpsymeqjumps, tmpsymopidx = apply_symmetry_operations(vector_list=vectors, symop_list=symops, unique=unique, symop_indices=True)  # symmetry equivalent jump vectors
        # Apply periodic boundary conditions to make sure the initial position of the first constraint is in the base supercell
        # Is this a problem? Have removed this step from initial jump_mech definition. 
        symeqjumps = []
        symopidx = []
        for idx, jump in enumerate(tmpsymeqjumps):
            tvec = -np.floor(jump[0] + tol_pbc * np.array([1, 1, 1]))
            jump = [np.array(i) + tvec for i in jump]  # translates each position
            # Modify dumbbell coordinates if needed (bulkspecies), first for initial configuration then for final
            if len(bulkspecies_list) > 0:
                jump[:self.__nconstraints] = apply_bulkspecies_positions(positions=jump[:self.__nconstraints], all_defect_list=all_defect_list, species=[c.get_inispecies() for c in self.__constraints])
                jump[self.__nconstraints:] = apply_bulkspecies_positions(positions=jump[self.__nconstraints:], all_defect_list=all_defect_list, species=[c.get_inispecies() for c in self.__constraints])
            # Remove jump duplicates
            found = False
            for comp in symeqjumps:
                if are_equal_arrays(np.array(jump, dtype=float), np.array(comp, dtype=float)):
                    found = True
                    break
            if not unique:
                found = False
            if not found:
                symeqjumps.append(jump)
                symopidx.append(tmpsymopidx[idx])

        del tmpsymeqjumps, tmpsymopidx

        # Add reverse jump if it is not already part of the symmetry equivalents (only for symmetry-unique list)
        # Check if reverse jump is already part of symetry equivalent jumps
        isreverseincluded = True
        revprt = ""
        #addreverse=True
        #logger.info(' addreverse ={}'. format(addreverse))

        if addreverse:
            # Checking species; if there are transmutations the reverse jump will not be part of symmetrical equivalents
            for jcons in self.__constraints:
                if jcons.get_inispecies() != jcons.get_finspecies():
                    isreverseincluded = False
                    break
            if isreverseincluded:
                isreverseincluded = False
                # There are no transmutations so reverse jump might be among symmetry equivalents, so let's check positions
                # Apply translation to make sure the atom of the first constraint is initially in the first supercell
                tvec = -np.floor(reversejump[0]+ tol_pbc*np.array([1, 1, 1]))
                for i in range(0, len(reversejump)):
                    reversejump[i] = reversejump[i] + tvec # translates each position
                # Modify dumbbell coordinates if needed (bulkspecies), first for initial configuration then for final
                if len(bulkspecies_list) > 0:
                    reversejump[:self.__nconstraints] = apply_bulkspecies_positions(positions=reversejump[:self.__nconstraints], all_defect_list=all_defect_list, species=[c.get_finspecies() for c in self.__constraints])
                    reversejump[self.__nconstraints:] = apply_bulkspecies_positions(positions=reversejump[self.__nconstraints:], all_defect_list=all_defect_list, species=[c.get_inispecies() for c in self.__constraints])
                for symeq in symeqjumps: # compare reverse jump with each symmetry equivalent jump
                    if are_equal_arrays(np.array(flat_list(symeq), dtype=float), np.array(flat_list(reversejump), dtype=float)):
                        isreverseincluded = True
                        break
            if not isreverseincluded:  # reverse jump was not found among symmetry equivalent and must be added
                
                revprt = "(!! Reverse jump added !!)"
                #logger.info('{}'. format(revprt))
                tmpsymeqjumps, tmpsymopidx = apply_symmetry_operations(vector_list=reversejump, symop_indices=True, symop_list=symops, unique=unique)
                for idx, jump in enumerate(tmpsymeqjumps):
                    tvec = -np.floor(jump[0] + tol_pbc * np.array([1, 1, 1]))
                    jump = [np.array(i) + tvec for i in jump]  # translates each position
                    # Modify dumbbell coordinates if needed (bulkspecies), first for initial configuration then for final
                    if len(bulkspecies_list) > 0:
                        jump[:self.__nconstraints] = apply_bulkspecies_positions(positions=jump[:self.__nconstraints], all_defect_list=all_defect_list, species=[c.get_finspecies() for c in self.__constraints])
                        jump[self.__nconstraints:] = apply_bulkspecies_positions(positions=jump[self.__nconstraints:], all_defect_list=all_defect_list, species=[c.get_inispecies() for c in self.__constraints])
                    # Remove jump duplicates
                    found = False
                    for comp in symeqjumps:
                        if are_equal_arrays(np.array(jump, dtype=float), np.array(comp, dtype=float)):
                            found = True
                            break
                    if not unique:
                        found = False
                    if not found:
                        symeqjumps.append(jump)
                        symopidx.append(tmpsymopidx[idx])
        #if unique and print:
            #logger.info('    Found {} symmetry equivalent jump mechanisms for {} {}'. format(len(symeqjumps), self.__name, revprt))
        # Create new jumps for all symmetry equivalents (including reverse jumps)
        for isym in range(0, len(symeqjumps)):
            symeq = symeqjumps[isym]
            species = []
            if (not isreverseincluded) and (isym >= 0.5*len(symeqjumps)):
                # Final species become initial species when jump is reversed
                species.append([e.get_finspecies() for e in self.get_constraints()])
                species.append([e.get_inispecies() for e in self.get_constraints()])
            else:
                species.append([e.get_inispecies() for e in self.get_constraints()])
                species.append([e.get_finspecies() for e in self.get_constraints()])
            new = JumpMech(self.__index, self.__name)
            # Add constraints to this JumpMech object
            for cons in range(0, self.__nconstraints):
                if (not isreverseincluded) and (isym >= 0.5 * len(symeqjumps)):
                    inipos_nod = apply_symmetry_operations(vector_list=[np.array(self.__constraints[cons].get_finposition_nod())], symop_list=[symops[symopidx[isym]]])[0][0]
                    finpos_nod = apply_symmetry_operations(vector_list=[np.array(self.__constraints[cons].get_iniposition_nod())], symop_list=[symops[symopidx[isym]]])[0][0]
                else:
                    inipos_nod = apply_symmetry_operations(vector_list=[np.array(self.__constraints[cons].get_iniposition_nod())], symop_list=[symops[symopidx[isym]]])[0][0]
                    finpos_nod = apply_symmetry_operations(vector_list=[np.array(self.__constraints[cons].get_finposition_nod())], symop_list=[symops[symopidx[isym]]])[0][0]
                if cons == 0:  # actualize the translation vector with the initial position of the first constraint
                    tvec = -np.floor(np.array(inipos_nod) + tol_pbc * np.array([1, 1, 1]))
                inipos_nod += tvec
                finpos_nod += tvec
                new.add_constraint(JumpConstraint(inispecies=species[0][cons], iniposition=symeq[cons], iniposition_nod=inipos_nod, finspecies=species[1][cons], finposition=symeq[cons+self.__nconstraints], finposition_nod=finpos_nod))
            # Retrieve the new jump from the already existing catalogue, or create a new one
            if old_jumps.get(new.get_label()) is None:
                jump_list.append(new)
            else:
                jump_list.append(old_jumps[new.get_label()])
        if symop_ind:
            return jump_list, symopidx
        else:
            return jump_list

    def check_jump(self, all_defect_list: list, component_list: list, crystal: Crystal):
        """Check that the jump constraints are fully consistent with each other, in terms of proper sublattices and no overlapping components"""
        ini_species = [c.get_inispecies() for c in self.get_constraints()]
        fin_species = [c.get_finspecies() for c in self.get_constraints() if c.get_finspecies().get_name() != 'bulk']
        ini_positions = [c.get_iniposition() for c in self.get_constraints()]
        fin_positions = [c.get_finposition() for c in self.get_constraints() if c.get_finspecies().get_name() != 'bulk']
        check_ini = check_subconfiguration_consistence(species_list=ini_species, position_list=ini_positions, all_defect_list=all_defect_list, component_list=component_list, crystal=crystal)
        check_fin = check_subconfiguration_consistence(species_list=fin_species, position_list=fin_positions, all_defect_list=all_defect_list, component_list=component_list, crystal=crystal)
        return (check_ini and check_fin)


class Species:
    """Defect species are provided by the user in the input file and are assigned a short and long name"""

    def __init__(self, index, name, defects, permissions, radius):
        self.__index = index  # species index (goes from 0 to n, there are n user-input species (from 0 to n-1), plus the bulk (n))
        self.__name = name  # species written name (user-input)
        self.__bulk = False  # flag marking if species resembles a bulk atom (e.g. a pure dumbbell)
        self.__defects = defects  # list of defects where the species can sit
        # For each allowed defect, self.__defect_permission is a list that gives 1 if the species can be in this defect by itself,
        # or -1 if it needs to have another defect in the same site (useful to treat mixed dumbbells) (default is 1)
        self.__permissions = permissions
        self.__radius = radius #Exclusion radius. No other species can be within the sum of this radius and the other species radius. 

    def set_bulk(self):
        self.__bulk = True

    def get_bulk(self):
        return self.__bulk

    def get_index(self):
        return self.__index

    def get_name(self):
        return self.__name

    def get_defects(self):
        defecttuple = self.__defects
        if not(type(defecttuple) is tuple):
            defecttuple = tuple(defecttuple)
        return defecttuple

    def get_permissions(self):
        permissions = self.__permissions
        if not(type(permissions) is tuple):
            permissions = tuple(permissions)
        return permissions

    def get_radius(self):
        return self.__radius

    def get_info(self):
        logger.info("  Species {} is a {} and can be in sublattices {}".format(self.__index, self.__name, [i.get_index() for i in self.__defects]))


class SymOp:
    """A Symop object is a 3x3 matrix (and a translation vector) representing a symmetry operation that is valid for the supercell"""

    def __init__(self, rotation):
        self.__rotation = rotation  # rotation matrix
        self.__translation = None  # 1-D translation vector
        self.__cpgsign = None  # +/-1 depending  on the transformation of CPG with this symetry operation

    def get_rotation(self):
        return self.__rotation

    def get_translation(self):
        return self.__translation

    def get_cpgsign(self):
        return self.__cpgsign

    def set_rotation(self, new_rot):
        self.__rotation = new_rot

    def set_translation(self, new_trans):
        self.__translation = new_trans

    def set_cpgsign(self, cpgsign):
        self.__cpgsign = cpgsign

    def apply_symop(self, vect):
        return np.dot(self.__rotation, vect) + self.__translation


# Definition of functions

def apply_pbc(vec: np.ndarray) -> np.ndarray:
    # applies periodic boundary conditions to numpy.array vec
    for r in range(0, 3):
        if vec[r] < (0 - tol_pbc):
            vec[r] += np.ceil(-vec[r])
        elif vec[r] >= (1 - tol_pbc):
            vec[r] += -np.floor(vec[r]+tol_pbc)
    return vec

def dupe_sort(t_list: list) -> list:
    n_dupes = len(t_list)
    if (len(t_list)<2):
        return t_list
    #Sorts list of list pairings by last entry. Intended for sorting duplicate jump index pairings
    r_list = []
    r_list.append(t_list[0])
    for t in range(1, len(t_list)):
        if (t_list[t][1] < r_list[0][1]):
            r_list.insert(0, t_list[t])
        else: 
            if (len(r_list) == 1):
                r_list.append(t_list[t])
            else:
                insert_flag = False 
                for r in range(1, len(r_list)):
                    if (t_list[t][1] < r_list[r][1]):
                        r_list.insert(r, t_list[t])
                        insert_flag = True
                        break
                if not insert_flag:
                    r_list.append(t_list[t])
    if (n_dupes != len(r_list)):
        "PROBLEM IN dupe_sort()!!!"
    return r_list


def apply_symmetry_operations(vector_list: list, symop_list: list, unique: bool = True, applytrans: bool = True, symop_indices: bool = False) -> list:
    """Apply a set (SymOp object list) of symmetry operations to a list of vectors, giving as output, for each symmetry operation,
    the list of equivalent vectors. For instance, if vector_list has 2 components and there are 48 symmetries, the output is a list
    of 48 elements, where each element is a list of 2 vectors. If unique=True, the number of symmetries is reduced if different symmetry operations
    yield the same result."""
    n_vec = len(vector_list)
    # Transform vector list in matrix (so to perform each sym operation once on the whole set of vectors)
    vector_matrix = np.column_stack(tuple(vector_list[i] for i in range(0,n_vec)))
    # Loop over list of symmetry operations
    symm_vec_matrix_list = []
    symm_vec_list_list = []
    symop_idx_list = []
    for idx, symop in enumerate(symop_list):
        # Apply symmetry operation
        if applytrans:
            # Transform translation vector into a 3xn matrix where the translation vector is repeated n times
            new_vector_matrix = np.dot(symop.get_rotation(), vector_matrix) + np.column_stack(tuple(symop.get_translation() for _ in range(0, n_vec)))
        else:
            new_vector_matrix = np.dot(symop.get_rotation(), vector_matrix)
        # Check whether this symmetry equivalent has already been found
        if unique:
            found = False
            for previous_found_matrix in symm_vec_matrix_list:
                if are_equal_arrays(previous_found_matrix, new_vector_matrix):
                    found = True
                    break
            if not found:
                symm_vec_matrix_list.append(new_vector_matrix)  # save new matrix in the database, for comparison in the next iterations
                symm_vec_list_list.append(new_vector_matrix.transpose().tolist())  # save the obtained list of vectors in the output
                symop_idx_list.append(idx)
        else:
            symm_vec_list_list.append(new_vector_matrix.transpose().tolist())  # save the obtained list of vectors in the output
            symop_idx_list.append(idx)
    if symop_indices:
        return symm_vec_list_list, symop_idx_list
    else:
        return symm_vec_list_list


def are_equal_arrays(A: np.ndarray, B: np.ndarray ) -> bool:
    """Check if two arrays (of any shape) have all equal elements (including tolerance)"""
    subtraction_array = np.fabs ( B - A )
    if np.all ( subtraction_array < tol ):
        return True
    else:
        return False

def are_equal_arrays_test(A: np.ndarray, B: np.ndarray ) -> bool:
    """Check if two arrays (of any shape) have all equal elements (including tolerance)"""
    subtraction_array = np.fabs ( B - A )
    if np.all ( subtraction_array < tol ):
        return True
    else:
        return False

def are_equals(a: float, b: float) -> bool:
    """Check if two numbers are equal based on tolerance defined in kinepy"""
    return np.fabs(a-b) < tol

def change_base(arr: np.ndarray, crystal: Crystal, inv: bool = False, matrix: bool = False) -> np.ndarray:
    """Change vector or matrix from orthonormal base to supercell base or viceversa"""
    # Select correct crystal matrix depending on "inv" variable (if True, pass from crystal to orthonormal base with invvectors)
    if inv:
        M = crystal.get_invvectors()
        Minv = crystal.get_primvectors()
    else:
        M = crystal.get_primvectors()
        Minv = crystal.get_invvectors()
    # Check if passed array is a vector or a matrix
    if matrix:
        return np.dot(Minv, np.dot(arr, M))
    else:
        return np.dot(Minv, arr)

def check_configuration_consistence(label: str, component_list: list, all_defect_list: list, crystal: Crystal) -> bool:
    """Check consistence of a configuration, given its label. The consistence criteria are the ones defined in check_subconfiguration_consistence"""
    spec_list = [comp.get_species() for comp in component_list]
    config = name2conf(name=label, all_defect_list=all_defect_list, species=spec_list)
    position_list = [config.get_defect_position(d) for d in range(0, len(component_list))]
    return check_subconfiguration_consistence(species_list=spec_list, position_list=position_list, all_defect_list=all_defect_list, crystal=crystal)  # do not pass component_list because we don't need to check for correct number of components


def check_subconfiguration_consistence(species_list: list, position_list: list, all_defect_list: list, crystal: Crystal, component_list: list = []) -> bool:
    """Check that list of species and corresponding positions is consistent with the system defined in the input.
    This entails checking on: 1) correct number of components (only if component_list has been passed) 2) sublattices,
    3) no overlapping on same sites, 4) special case of "-1" permission."""
    # Check consistency of input (length of lists)
    if len(species_list) != len(position_list):
        produce_error_and_quit("In check_subconfiguration_consistence - Species list and position list have different lengths.")
    # Check if we have the correct number of components (only if the full component list has been passed)
    if len(component_list) > 0:
        component_species_list = [comp.get_species().get_name() for comp in component_list]
        for spec in species_list:
            if spec.get_name() != 'bulk':
                if spec.get_name() in component_species_list:
                    component_species_list.remove(spec.get_name())  # remove used component from list
                else:
                    return False
    # Check if species sit on the correct sublattices
    def_list = [] # list of defects occupied in this configuration
    perm_list = []  # list of permissions for each found defect
    for spec, pos in zip(species_list, position_list):
        if spec.get_name() != 'bulk':  # skip bulk species, assuming it is not a problem if user inputs a wrong lattice site
            # Transform position into a defect+translation pair (but translation is not needed)
            defect, _ = vect2deftrans(vec=pos, all_defect_list=all_defect_list)
            if defect.get_sudefect() not in spec.get_defects():
                return False
            else:  # locate and add corresponding permission
                perm_list.append((defect.get_sudefect(), spec.get_permissions()[spec.get_defects().index(defect.get_sudefect())]))
        else:  # for bulk species, add an artificial permission (just to maintain the correct species order in perm_list)
            perm_list.append((None, 1))
    if len(position_list) == 1:  # simply check that there are no negative permissions
        if not perm_list[0][1] > 0:
            return False
    
    # If all checks were successful at this point, the configuration is acceptable!
    return True

def dataset_to_defects(dataset: dict, crystal: Crystal, symop_list: list) -> list:
    n_defects = int(dataset['uniquepos'][0])
    defect_list = []  # list of all symmetry unique defects
    all_defect_list = []  # list of all defects, including symmetric ones
    dindex = 0
    doubledefects = []
    for i in range(n_defects):
        uniquepos_vec = dataset['uniquepos'][i * 3 + 2:i * 3 + 4]
        uniquepos_vec.append('0')
        sublattice = np.array([evalinput(e) for e in uniquepos_vec], dtype=float)
        sublattice_nod = np.array([evalinputd0(e) for e in uniquepos_vec], dtype=float)
        if dataset['uniquepos'][i * 3 + 1] == 'o':  # convert positions to supercell base if necessary
            sublattice = change_base(arr=sublattice, crystal=crystal)
            sublattice_nod = change_base(arr=sublattice_nod, crystal=crystal)
        sublattice = apply_pbc(vec=sublattice)  # apply boundary conditions
        sublattice_nod = apply_pbc(vec=sublattice_nod)
        defect_list.append(Defect(crystal=crystal, sublattice=sublattice, index=dindex, sublattice_nod=sublattice_nod))  # add to list of symmetry-unique defects
        defect_list[i].find_symeq(symop_list)  # find all symmetry equivalents of this defect
        defect_list[i].set_align_dist()
        doubles = False
        # Check that this defect is not equivalent to a previously found defect
        for new in defect_list[i].get_symeq():
            for iold, old in enumerate(defect_list[:-1]):
                if (np.allclose(new.get_sublattice(), old.get_sublattice(), atol=tol, rtol=tol)):
                    doubles = True
                    logger.info('!! Defect {} is ignored because symmetrically equivalent to defect {}.'.format(i,old.get_index()))
                    doubledefects.append(i)
                    break
            if doubles:
                break
        if not doubles:
            all_defect_list += defect_list[i].get_symeq()
            dindex += len(defect_list[i].get_symeq())
            #defect_list[i].get_info()
    doubledefects=list(reversed(sorted(doubledefects)))
    for i in doubledefects:
        del defect_list[i]
    return [defect_list, all_defect_list, doubledefects]


def dataset_to_jumps_cat(dataset: dict, crystal: Crystal, symop_list: list, all_defect_list: list, species_list: list, sym_unique: bool, component_list: list) -> list:
    su_jump_list = []  # list of symmetry unique jump objects
    jump_list = []  # list of all possible jump objects
    jump_symop_ind_list = [] # list of symmetry operations valid for each jump object
    jumpmech2 = (''.join(str(e) + ' ' for e in dataset['procmech'])).split(sep='%%')[1:]  # '%%' symbol separates process mechanisms
    #print('jumpmech2:',jumpmech2)
    for i in range(len(jumpmech2)-1,-1,-1):  # Checking the number of species for this jump
        tmpcomp = [a.get_species().get_index() for a in component_list] # non-bulk species in cluster
        tmp_component_list=[] #Collects both initial and final species for all constraints.
        tmplist=[a.split()[0] for a in jumpmech2[i].split(sep='>')[1:]]
        for a in tmplist:
            tmp_component_list.append(a)
        tmplist=[a.split()[-1] for a in jumpmech2[i].split(sep='>')[0:(len(jumpmech2[i].split(sep='>'))-1)]]
        for a in tmplist:
            tmp_component_list.append(a)
        for sp in [int(b) for b in tmp_component_list if b !='0']: # non bulk species in jump constraints
            if not sp in tmpcomp: #just consider if species in list, not if number matches.
                logger.info("!! Jump {} is not taken into account because it requires components that are not part of the cluster".format(jumpmech2[i].split(sep=' ')[2]))
                del jumpmech2[i]
                break
    n_jumps = len(jumpmech2)
    #logger.info("  Found {} process mechanisms".format(n_jumps))
    for i in range(0, n_jumps):  # Loop on jump mechanisms
        jumpmech3 = jumpmech2[i].split(sep=' ')[1:]
        su_jump_list.append(JumpMech(index=i, name=jumpmech3[1]))  # List of symmetry unique jump mechanisms
        # Check user-input constraints for this jump
        tvec = np.array([0, 0, 0], dtype=float)
        tmplist = []  # temporary constraint list
        #print('jumpmech3:', jumpmech3)
        for j in range(0, int(jumpmech3[0])):  # Loop on constraints for this jump
            # Check for correct constraint format in input file
            #print('jumpmech3 slice',jumpmech3[6 + j * 6])
            if jumpmech3[6 + j * 6] != '>':
                produce_error_and_quit("In dataset_to_jumps - Problem in the format of constraints {} of jump {}".format(j+1, i+1))
            # If everything looks fine, add this constraint to the temporary list
            tmplist.append(jumpmech3[2 + j * 6:8 + j * 6])
        #print('tmplist:',tmplist)
        jumpmech3 = jumpmech3[0:2] + flat_list(sorted(tmplist, key=lambda x: -int(x[3])))  # sort constraints by species in decreasing order (from highest index to 0)
        # Add constraints to jump object
        for j in range(0, int(jumpmech3[0])):  # Loop on constraints for this jump
            # Initial species and position
            inispec = species_list[int(jumpmech3[5 + j * 6]) - 1]
            #print('inispec:', inispec.get_name())
            pos_vec = jumpmech3[3 + j * 6:5 + j * 6]
            pos_vec.append('0')
            #print('inipos', pos_vec)
            inipos = np.array([evalinput(e) for e in pos_vec], dtype=float)
            inipos_nod = np.array([evalinputd0(e) for e in pos_vec], dtype=float)
            if jumpmech3[2 + j * 6] == 'o':  # convert to supercell basis if necessary
                inipos = change_base(arr=inipos, crystal=crystal)
                inipos_nod = change_base(arr=inipos_nod, crystal=crystal)
            
            ## Translate every atom so that the initial position of the atom in the first constraint is in the initial supercell
            #N.B. No longer checking if in initial supercell. Trusting user to intelligently define constraints.
            
            # Final species and position
            finspec = species_list[int(jumpmech3[7 + j * 6]) - 1]
            #print('finspec:', finspec.get_name())
            finpos = np.array([evalinput(e) for e in pos_vec], dtype=float)
            #print('finpos', finpos)
            finpos_nod = np.array([evalinputd0(e) for e in pos_vec], dtype=float)
            if jumpmech3[2 + j * 6] == 'o':  # convert to supercell basis if necessary
                finpos = change_base(arr=finpos, crystal=crystal)
                finpos_nod = change_base(arr=finpos_nod, crystal=crystal)
            finpos += tvec  # translating (as above)
            finpos_nod += tvec  # translating (as above)
            # Add constraint object to this jump object
            su_jump_list[i].add_constraint(JumpConstraint(inispecies=inispec, iniposition=copy.copy(inipos), iniposition_nod=copy.copy(inipos_nod), finspecies=finspec, finposition=copy.copy(finpos), finposition_nod=copy.copy(finpos_nod)))
        # Check jump, making sure that species sit on correct sublattices and there are no overlapping components on the same site (unless permission -1 is present)
        if not su_jump_list[i].check_jump(all_defect_list=all_defect_list, component_list=component_list, crystal=crystal):
            produce_error_and_quit("In dataset_to_jumps - Jump {} is not consistent with sublattice permissions and/or site occupancy.".format(su_jump_list[i].get_name()))
        # Search and assign defect type and translation to each constraint
        for cons in su_jump_list[i].get_constraints():
            cons.vect2deftrans_constraint(all_defect_list=all_defect_list)
        # Check that this jump is not symmetrically equivalent to other jumps (but do not remove it!)
        #for poteq in jump_list:
        #    if su_jump_list[i] == poteq:
        #        logger.info("!! Jumps {} and {} are symmetrically equivalent.".format(str(i), str(poteq.get_index())))
        #        logger.info("Ignore for now. Likely not symmetrically equivalent, just due to altered event mechanism input format.")
        
        tmp_jump_list, jump_symop_ind=su_jump_list[i].find_symeqs(symops=symop_list, all_defect_list=all_defect_list, addreverse=False ,symop_ind=True)
        jump_list += tmp_jump_list
        jump_symop_ind_list.append(jump_symop_ind)
    if sym_unique:
        return su_jump_list
    else:
        
        #This part is to ensure that the full list of events does not include jumps that are symmetric with a different jump from a different cell (no double counting)
        n_spec=len(species_list)
        duplicate_jumps = []
        for j, jump in enumerate(jump_list):
            dupe_flag = False
            for dupe_ind in range(len(duplicate_jumps)):
                if duplicate_jumps[dupe_ind][1] == j:
                    dupe_flag = True
                    break
            if dupe_flag:
                #This search & skip also prevents multiple references to the same jumps in duplicate_jumps list (occurs if triplicate or more jumps)
                continue
            for k in range(j+1,len(jump_list)):
                if jump.get_name()==jump_list[k].get_name():
                    temp_cons1=jump.get_constraints()
                    temp_cons2=jump_list[k].get_constraints()
                    inispec_coord_lists1=[[] for _ in range(n_spec)]
                    finspec_coord_lists1=[[] for _ in range(n_spec)]
                    inispec_coord_lists2=[[] for _ in range(n_spec)]
                    finspec_coord_lists2=[[] for _ in range(n_spec)]
                    for con in temp_cons1:
                        spec=con.get_inispecies().get_index()
                        inispec_coord_lists1[spec].append(con.get_iniposition().tolist())
                        spec=con.get_finspecies().get_index()
                        finspec_coord_lists1[spec].append(con.get_finposition().tolist())
                    for spec in range(n_spec):
                        inispec_coord_lists1[spec]=sortcoords(inispec_coord_lists1[spec])
                        finspec_coord_lists1[spec]=sortcoords(finspec_coord_lists1[spec]) 
                    for con in temp_cons2:
                        spec= con.get_inispecies().get_index()
                        inispec_coord_lists2[spec].append(con.get_iniposition().tolist())
                        spec=con.get_finspecies().get_index()
                        finspec_coord_lists2[spec].append(con.get_finposition().tolist())
                    for spec in range(n_spec):
                        inispec_coord_lists2[spec]=sortcoords(inispec_coord_lists2[spec])
                        finspec_coord_lists2[spec]=sortcoords(finspec_coord_lists2[spec])
                        
                    break_flag = False
                    prior_x_shift = -10
                    prior_y_shift = -10
                    for spec in range(n_spec):
                        if break_flag:
                            continue
                        ini_shift_x=0.0
                        ini_shift_y=0.0
                        fin_shift_x=0.0
                        fin_shift_y=0.0
                        if len(inispec_coord_lists1[spec]) == 0: #-10 becomes marker for no constraints for that species & ini/fin (e.g. not present)
                            ini_shift_x = -10
                            ini_shift_y = -10
                        if len(finspec_coord_lists1[spec]) == 0:
                            fin_shift_x = -10
                            fin_shift_y = -10
                        for c in range(len(inispec_coord_lists1[spec])):#Assumes that constraints have same number of initial and final species coordinates. 
                            ini_shift_x+=((inispec_coord_lists1[spec][c][0]-inispec_coord_lists2[spec][c][0])/(len(inispec_coord_lists1[spec])))
                            ini_shift_y+=((inispec_coord_lists1[spec][c][1]-inispec_coord_lists2[spec][c][1])/(len(inispec_coord_lists1[spec])))
                        for c in range(len(finspec_coord_lists1[spec])):
                            fin_shift_x+=((finspec_coord_lists1[spec][c][0]-finspec_coord_lists2[spec][c][0])/(len(finspec_coord_lists1[spec])))
                            fin_shift_y+=((finspec_coord_lists1[spec][c][1]-finspec_coord_lists2[spec][c][1])/(len(finspec_coord_lists1[spec])))
                        #NB: Shifts only refers to center of mass for species
                        temp=(ini_shift_x%1)
                        if not are_equals(temp,0.0): #Only checking for whole unit cell shifts.
                            break_flag = True
                            continue
                        temp=(ini_shift_y%1)
                        if not are_equals(temp,0.0):
                            break_flag = True
                            continue
                        temp=(fin_shift_x%1)
                        if not are_equals(temp,0.0):
                            break_flag = True
                            continue
                        temp=(fin_shift_y%1)
                        if not are_equals(temp,0.0):
                            break_flag = True
                            continue
                        #If to this point, only whole cell shifts
                        if not ((ini_shift_x == fin_shift_x) or (ini_shift_x == -10 or fin_shift_x == -10)):
                            break_flag = True
                        if not ((ini_shift_y == fin_shift_y) or (ini_shift_x == -10 or fin_shift_x == -10)):
                            break_flag = True
                        if not break_flag: #This section ensures that any whole lattice shifts are the same between initial and final and for all species
                            if ini_shift_x != -10:
                                if prior_x_shift == -10:
                                    prior_x_shift = ini_shift_x
                                else:
                                    if (prior_x_shift != ini_shift_x):
                                        break_flag = True
                            if fin_shift_x != -10:
                                if prior_x_shift == -10:
                                    prior_x_shift = fin_shift_x
                                else:
                                    if (prior_x_shift != fin_shift_x):
                                        break_flag = True
                            if ini_shift_y != -10:
                                if prior_y_shift == -10:
                                    prior_y_shift = ini_shift_y
                                else:
                                    if (prior_y_shift != ini_shift_y):
                                        break_flag = True
                            if fin_shift_y != -10:
                                if prior_y_shift == -10:
                                    prior_y_shift = fin_shift_y
                                else:
                                    if (prior_y_shift != fin_shift_y):
                                        break_flag = True 
                        if not break_flag: #This section ensures that shifts do result in same jumps, not just in jumps with shifted centers of mass (e.g. up-right vs. down-right diagonals)
                            for c in range(len(inispec_coord_lists1[spec])):
                                if (inispec_coord_lists1[spec][c][0] != (inispec_coord_lists2[spec][c][0] + ini_shift_x) 
                                    or (inispec_coord_lists1[spec][c][1] != (inispec_coord_lists2[spec][c][1] + ini_shift_y))):
                                        break_flag = True 
                                        continue
                            for c in range(len(finspec_coord_lists1[spec])):
                                if ((finspec_coord_lists1[spec][c][0] != (finspec_coord_lists2[spec][c][0] + fin_shift_x))
                                    or (finspec_coord_lists1[spec][c][1] != (finspec_coord_lists2[spec][c][1] + fin_shift_y))):
                                        break_flag = True 
                                        continue 
                                
                    if not break_flag:
                        duplicate_jumps.append([j,k]) # j is 1st jump, k is duplicate jump to be removed
        
        duplicate_jumps = dupe_sort(duplicate_jumps) #Sort to ensure that jumps to be removed are in increasing order
        for i in range(len(duplicate_jumps)):
            jump_to_remove = duplicate_jumps[len(duplicate_jumps)-i-1][1]
            jump_list.remove(jump_list[jump_to_remove])
            #Removed jump from list, need to remove corresponding symop_ind
            symop_count = 0
            for j in range(len(jump_symop_ind_list)):
                if jump_to_remove < (symop_count+ len(jump_symop_ind_list[j])):
                    jump_symop_ind_list[j].remove(jump_symop_ind_list[j][jump_to_remove-symop_count])
                    break
                else:
                    symop_count += len(jump_symop_ind_list[j])
        return jump_list, jump_symop_ind_list


def dataset_to_species(dataset: dict, defect_list: list, doubledefects: list = []) -> list:
    n_species = int(dataset['species'][0])
    n_defects = int(dataset['uniquepos'][0])  # n_defects is not taken from len(defect_list) because doubledefects have already been removed from that list
    n_components = 0
    species_list = []  # list of species
    component_list = []  # list of components
    for i in range(0, n_species):
        # For each species, create the list of defects that are allowed, based on user input
        tmp_permissions = dataset['species'][i * (3 + n_defects) + 2:i * (3 + n_defects) + (2 + n_defects)]
        if len(tmp_permissions) != (n_defects):  # wrong list of permissions in user input file!
            produce_error_and_quit("In dataset_to_species - The number of permissions for each species does not match the number of defects in UNIQUEPOS!\n" \
                                   "Correct input for each species: << X  0/1 0/1 ... spec_name >>,\n" \
                                   "where X is the number of components for that species, and 0/1 marks if the species can occupy the corresponding position in UNIQUEPOS.\n" \
                                   "Check that the amount of 0/1 flags matches the amount of defects defined in UNIQUEPOS, \n" \
                                   "Check that species exclusion radius is entered (last float, after name).")
        for idef in doubledefects:
            del tmp_permissions[idef]
        tmp_bool = np.array(tmp_permissions, dtype=int).astype(bool)  # OBS! anything that is not 0 returns True!
        tmp_defectlist = [k for (k, v) in zip(defect_list, tmp_bool) if v]  # selection of items in defect_list corresponding to a "True" value in tmp_bool
        tmp_reduced_permissions = np.array([k for (k, v) in zip(tmp_permissions, tmp_bool) if v], dtype=int)
        tmp_name = dataset['species'][i * (3 + n_defects) + (2 + n_defects)].lower()
        tmp_radius = np.float(dataset['species'][i * (3 + n_defects) + (2 + n_defects) +1])

        # Check that there is at least one positive permission, otherwise give warning
        if len(tmp_reduced_permissions) == 0:
            logging.info("WARNING! Species {} is not permitted on any sublattice. Is this what you meant doing?".format(tmp_name))
        # Check species name (name 'bulk' is forbidden because it is assigned to the matrix atoms)
        if tmp_name == 'bulk':
            produce_error_and_quit("Name 'bulk' for a species is forbidden. Choose another name, and keep in mind that you don't need to define the matrix atoms as a species in the input file.")
        # Append new species to species list
        species_list.append(Species(index=i + 1, name=tmp_name, defects=tmp_defectlist, permissions=[int(p) for p in tmp_reduced_permissions],radius=tmp_radius))
        #species_list[-1].get_info()
        # Check negative permissions (-n, where n is the species that has to be on the same site in order to activate permission)
        # Create adequate number of components for each species
        for j in range(0, int(dataset['species'][i * (3 + n_defects) + 1])):
            component_list.append(Component(species=species_list[i], index=n_components))
            #component_list[n_components].get_info()
            n_components += 1

    # Create bulk (matrix) species (for jump constraints)
    species_list.append(Species(index=0, name="bulk", defects=[], permissions=[],radius=0))
    return [species_list, component_list]

def distance(vect1: np.ndarray, vect2: np.ndarray, crystal: Crystal) -> float:
    """Returns the distance between vectors 1 and 2"""
    # First the difference between both vectors is transformed into the orthonormal basis
    ortho_dist = change_base(arr=vect2 - vect1, crystal=crystal, inv=True)
    # In strain calculations, ortho_dist contains the 'e' variable, but we take the distance at e=0
    if ortho_dist.dtype == 'O':
        ortho_dist = np.array([x.subs('e',0) for x in ortho_dist], dtype=float)
    return np.linalg.norm(ortho_dist)

def evalinputS(number: str) -> object: # for symbolic strain
    return eval(number.replace('sqrt', 'sym.sqrt'))


def evalinputd0(number: str) -> float:
    return float(eval(number.replace("d", str(0)).replace('sqrt', 'np.sqrt')))


def evalinput(number: str) -> float:
    return float(eval(number.replace("d", str(db_displ)).replace('sqrt', 'np.sqrt')))

def find_symmetry_operations(crystal: Crystal) -> list:
    """Gives a list of all valid symmetry operations (SymOp class) of a given crystal"""
    # Search for all possible rotation matrices is done with a recursive function that creates all possible 3x3 matrices with (-1,0,1)
    # by updating the list of symop as a function attribute
    def build_unitary_matrices(ix, args_list) -> None:
        if ix <= 9:  # Call function recursively until 9 indices are defined
            for i in range(-1,2):
                new_args_list = args_list.copy()
                new_args_list.append(i)  # a new index is added to args_list
                build_unitary_matrices(ix+1, new_args_list)
        else:  # Now that we have 9 indices, check if matrix can be a symmetry operation (determinant = 1 and 6-idempotent)
            Nmat = np.array(args_list, dtype=int).reshape((3,3))  # Building a 3x3 matrix from the 9 indices
            determ = np.linalg.det(Nmat)
            if are_equals(np.fabs(determ), 1):  # if determinant is 1 or -1
                for j in range(1, 7):  # check if 6-idempotent
                    Nmat2 = np.linalg.matrix_power(Nmat, j)    # Nmat2 = np.matrix(Nmat) ** j
                    if are_equal_arrays(Nmat2, np.identity(3)):
                        build_unitary_matrices.symop_list.append(SymOp(rotation=Nmat))  # Store found matrix

    # Here the build_unitary_matrices recursive function is called
    build_unitary_matrices.symop_list = []  # initialization of the symop_list as function attribute, which is updated inside the function
    build_unitary_matrices(ix=1, args_list=[])
    symop_list = build_unitary_matrices.symop_list

    # Now check if At * A = id, where A = C * M * C^-1, and keep only matrices fulfilling this criterion
    # (loop is done backward so that matrices can be removed safely without interfering with the loop)
    C = crystal.get_primvectorsnum()
    Cinv = crystal.get_invvectorsnum()
    for j in range(len(symop_list)-1, -1, -1 ):
        M = symop_list[j].get_rotation()
        A = np.dot(C, np.dot(M, Cinv))  # C*M*C^-1
        if not are_equal_arrays(np.dot(np.transpose(A), A), np.identity(3)):  # Is At*A - id == 0 ?
            del symop_list[j]   # matrix is removed

    # Modify matrices for 2D or 1D crystals
    if crystal.get_dim() <= 2:
        for j in range(len(symop_list) - 1, -1, -1):
            A = symop_list[j].get_rotation()
            A[2,:] = np.array([0,0,1])
            A[:,2] = np.array([0,0,1])
            symop_list[j].set_rotation(A)
    if crystal.get_dim() == 1:
        for j in range(len(symop_list) - 1, -1, -1):
            A = symop_list[j].get_rotation()
            A[1,:] = np.array([0,1,0])
            A[:,1] = np.array([0,1,0])
            symop_list[j].set_rotation(A)

    # Finally, check for identical matrices and remove them (backward loop for the same reason as above)
    for j in range(len(symop_list)-1, -1, -1):
        for k in range(len(symop_list)-1, j, -1):
            if are_equal_arrays(symop_list[k].get_rotation(), symop_list[j].get_rotation()):
                del symop_list[k]
    
    # We have symmetry matrices for the lattice but not yet for the crystal. Need to look for fractional translations
    SGop = 0  # number of space group operations
    PGop = 0 # number of point group operations
    basis = crystal.get_basislist()

    if len(basis) > 1:
        basis = np.array(sorted([a for a in basis], key=lambda x: [x[0], x[1], x[2], x[3]]))
        eqbasis = apply_symmetry_operations(vector_list=[a[1:4] for a in basis], symop_list=symop_list, unique=False, applytrans=False)
        symtoremove=[] # indexes of symmetry operations to be removed if any
        for idx, eqpositions in enumerate(eqbasis):  # loop on symmetry operations
            for id2, pos in enumerate(eqpositions):  # loop on basis atoms
                #apply periodic boundary conditions
                eqpositions[id2] = apply_pbc(np.array(pos))
            # find possible translations on the symmetry equivalent of the first atom
            # this symmetry equivalent of the first basis atom must coincide with one of the basis atoms
            # hence the possible translations; same goes for each and every symmetry operations
            possible_trans = []
            for pos in basis:
                if pos[0]==basis[0][0]:
                    possible_trans.append(pos[1:4]-eqpositions[0])
            found = False
            for trans in possible_trans:
                translated_symeq = []
                for pos in eqpositions:
                    translated_symeq.append(apply_pbc(pos+trans))
                      
                temp_array=np.array(sorted([np.array([basis[i][0]]+list(a)) for i,a in enumerate(translated_symeq)], key=lambda x:[x[0], x[1], x[2], x[3]]))
                temp_array2= sorted([np.array([basis[i][0]]+list(a)) for i,a in enumerate(translated_symeq)], key=lambda x:[x[0], x[1], x[2], x[3]])
                temp_array3=[]
                for temp_a in temp_array2:
                    temp_a = temp_a.tolist()
                #Correctly sorts 1st (maybe 2nd) terms, but not guaranteed beyond. Need to resort.
                temp_array = []
                temp_array.append(temp_array2[0][1:].tolist())
                for a in range(1,len(temp_array2)):
                    if are_equals(temp_array2[a][0], temp_array2[a-1][0]):
                        temp_array.append(temp_array2[a][1:].tolist())
                    else:
                        if (len(temp_array) == 1):
                            sort_array = temp_array
                        else:
                            sort_array = sortcoords(temp_array)
                        for b_int in range(len(sort_array)):
                            temp_array3.append([temp_array2[a-1][0], sort_array[b_int][0], sort_array[b_int][1], sort_array[b_int][2]])
                        temp_array = []
                        temp_array.append(temp_array2[a][1:].tolist())
                if (len(temp_array) == 1):
                    sort_array = temp_array
                    temp_array3.append([temp_array2[a][0], sort_array[b_int][0], sort_array[b_int][1], sort_array[b_int][2]])
                else:
                    sort_array = sortcoords(temp_array)
                    for b_int in range(len(sort_array)):
                        temp_array3.append([temp_array2[a-1][0], sort_array[b_int][0], sort_array[b_int][1], sort_array[b_int][2]])

                if are_equal_arrays(np.array(temp_array3), basis):
                    found = True
                    #print("found")
                    if symop_list[idx].get_translation() is None:
                        symop_list[idx].set_translation(trans)

                    else:
                        symop_list.append(SymOp(rotation=symop_list[idx].get_rotation()))
                        symop_list[-1].set_translation(trans)
                    if are_equal_arrays(trans, np.array([0,0,0])):
                        PGop += 1
                    else:
                        SGop += 1
            if not found:
                symtoremove.append(idx)
        symtoremove.sort(reverse=True)
        for idx in symtoremove:
            del symop_list[idx]
    else:
        for idx in range(len(symop_list)):
            PGop += 1
            symop_list[idx].set_translation(np.array([0,0,0]))
    logger.info("  Found {} symmetry operations ({} point group op. and {} space groups op.)".format(len(symop_list),PGop,SGop))
    #print("  Found {} symmetry operations ({} point group op. and {} space groups op.)".format(len(symop_list),PGop,SGop))
    return symop_list


def find_possible_permutations(specs: list, species_list):
    # Creating permutations list "name_perms" for configuration names
    # specs is a list of species for each component
    # species_list is the list of all species in the system, the last one being bulk
    tmp = []
    for isp, sp in enumerate(species_list[:-1]):  # loop over each species, except the last one which is bulk
        indices = []
        for i, x in enumerate(specs):
            if x == sp:
                indices.append(i)
        tmp.append(list(permutations(indices)))
    name_perms = list(product(*tmp))
    for id, variation in enumerate(name_perms):
        name_perms[id] = flat_list(variation)
    return name_perms


def flat_list(listoflist) -> list:
    """Takes a list of lists and turns it into a simple list"""
    flat_list = [inner for outer in listoflist for inner in outer]
    return flat_list


def name2conf(name: str, all_defect_list: list, species: list) -> Configuration:
    tmpname = name.replace('c', '').replace('d', '').split(sep='|')[0].split(sep='_')
    translations=[0, 0, 0] + [int(a) for a in tmpname[-1].replace('m','p-').split(sep='p')[1:]]
    defects = [all_defect_list[int(a)] for a in tmpname[:-1]]
    conf = Configuration(defects=copy.copy(defects), translations=copy.copy(translations), label=name)
    conf.set_species(species)
    return conf


def confsites(event_coords: list, kira: float, crystal: Crystal, species: list)->list: #Finds all permutations of species on lattice, saves them in config_list
    n = int(10) #estimate of maximum range needed to explore lattice.
    limitp = [n, n, 0]
    limitm = [-n, -n, 0]
    coords_list=[]
    for superx in range(limitm[0], limitp[0]+1):
        for supery in range(limitm[1], limitp[1]+1):
            trans = np.array([superx, supery, 0])
            for bulk in crystal.get_basislist():
                tmp_coord=bulk[1:4]+trans
                for coord in event_coords:
                    if distance(tmp_coord,np.array(coord),crystal=crystal)<=kira:
                        coords_list.append(bulk[1:4]+trans)
                        break
    #print("found ",len(coords_list),"lattice sites within range") #Need to check how sublattices/basis atoms work
    return coords_list


def printxyz(direc: str, input: list, thconfig_list: dict, freq_list: dict, n_components: int, crystal: Crystal, component_list: list, name_perms: list) -> None:
    if not os.path.isdir(direc+'CONFIGURATIONS'):
        os.makedirs(direc+'CONFIGURATIONS')
    k = 0
    for conf in thconfig_list:
        with open(direc+'CONFIGURATIONS/conf_{}.dat'.format(thconfig_list[conf].get_thermoint()), 'w') as output:
            limitp = [-1, -1, -1]
            limitm = [0, 0, 0]
            if len(input) == 2:
                if input[1] == 'wrap':
                    for b in range(3):
                        limitp[b] = 1+int(np.max(a=[thconfig_list[conf].get_translations()[3*e+b] for e in range(n_components)]))
                        limitm[b] = -1+int(np.min(a=[thconfig_list[conf].get_translations()[3*e+b] for e in range(n_components)]))
                else:
                    try:
                        n = int(float(input[1]))
                        limitp = [n, n, n]
                        limitm = [-n, -n, -n]
                    except ValueError:
                        logger.info("!! Second item in PRINTXYZ tag must be 'wrap' or an integer")
            output.writelines('{:.0f}\n'.format(n_components+(limitp[0]+1-limitm[0])*(limitp[1]+1-limitm[1])*(limitp[2]+1-limitm[2])*len(crystal.get_basislist())))
            if input[0] == 'o':  # orthonormal basis
                cell = "".join('{:8.4f} {:8.4f} {:8.4f}\n {:8.4f} {:8.4f} {:8.4f}\n {:8.4f} {:8.4f} {:8.4f}'.format(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
            elif input[0] == 's':  # supercell basis
                cell = "".join(['{:8.4f} {:8.4f} {:8.4f}\n '.format(*(crystal.get_primvectors()[b]*(limitp[b]-2*limitm[b]))) for b in range(3)])
                output.writelines(cell)
            else:
                logger.info('!! Problem in the PRINTXYZ tag; must be s (supercell basis) or o (orthonormal basis)')
            #output.writelines("Lattice=\"{}\" Properties=species:S:1:pos:R:3 Time=0.0\n".format(cell))
            output.writelines("\n")
            if input[0] == 'o':  # orthonormal basis
                for idx in range(n_components):
                    output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}\n'.format(component_list[idx].get_species().get_name(),*change_base(arr=thconfig_list[conf].get_defect_position(def_idx=idx),crystal=crystal, inv=True)))
            elif input[0] == 's':  # supercell basis
                for idx in range(n_components):
                    output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}\n'.format(component_list[idx].get_species().get_name(),*thconfig_list[conf].get_defect_position(def_idx=idx)))
            for superx in range(limitm[0], limitp[0]+1):
                for supery in range(limitm[1], limitp[1]+1):
                    for superz in range(limitm[2], limitp[2]+1):
                        trans = np.array([superx, supery, superz])
                        for bulk in crystal.get_basislist():
                            if input[0] == 'o':  # orthonormal basis
                                output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}\n'.format('Bulk{:.0f}'.format(bulk[0]),*change_base(arr=(bulk[1:4]+trans), crystal=crystal, inv=True)))
                            elif input[0] == 's':  # supercell basis
                                output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}\n'.format('Bulk{:.0f}'.format(bulk[0]),*(bulk[1:4]+trans)))
        if k == 0:
            k = 1
            if not os.path.isdir(direc+'JUMP_FREQUENCIES'):
                os.makedirs(direc+'JUMP_FREQUENCIES')
            for freq0 in freq_list:
                for freq in  freq_list[freq0]:
                    with open(direc+'JUMP_FREQUENCIES/freq_{}.dat'.format(freq.get_number()), 'w') as output:
                        if input[0] == 's': # print crystal vector only when using supercell basis
                            for b in range(3):
                                output.writelines('{:8.4f} {:8.4f} {:8.4f}\n'.format(*crystal.get_primvectors()[b]))
                        limitp = [-1, -1, -1]
                        limitm = [0, 0, 0]
                        ini_position, fin_position = untranslatedfinal(jumpmech=freq.get_jump(), component_list=component_list, name_perms=name_perms,
                                                    finconf=freq.get_config_fin(), iniconf=freq.get_config_ini(), fnum=freq.get_number())
                        if len(input) == 2:
                            if input[1] == 'wrap':  # must be changed if we change the output
                                for b in range(3):
                                    limitp[b] = int(np.max([np.floor(flat_list(flat_list([ini_position, fin_position]))[3*e+b]+tol_pbc) for e in range(2*n_components)]))
                                    limitm[b] = int(np.min([np.floor(flat_list(flat_list([ini_position, fin_position]))[3*e+b]+tol_pbc) for e in range(2*n_components)]))
                            else:
                                try:
                                    n = int(float(input[1]))
                                    limitp = [n, n, n]
                                    limitm = [-n, -n, -n]
                                except ValueError:
                                    logger.info("!! Second item in PRINTXYZ tag must be 'wrap' or an integer")
                        output.writelines('{:.0f}\n'.format(n_components + (limitp[0]+1-limitm[0])*(limitp[1]+1-limitm[1])*(limitp[2]+1-limitm[2])*len(crystal.get_basislist())))
                        output.writelines('\n')
                        # Print coordinates
                        if input[0] == 'o':  # orthonormal basis
                            for idx in range(n_components):
                                output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}  > {:8.4f} {:8.4f} {:8.4f}\n'.format(component_list[idx].get_species().get_name(),
                                    *change_base(arr=ini_position[idx], crystal=crystal,inv=True),
                                    *change_base(arr=fin_position[idx], crystal=crystal, inv=True)))
                        elif input[0] == 's': # supercell basis
                            for idx in range(n_components):
                                output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}  > {:8.4f} {:8.4f} {:8.4f}\n'.format(component_list[idx].get_species().get_name(),
                                    *ini_position[idx], *fin_position[idx]))
                        for superx in range(limitm[0], limitp[0]+1):
                            for supery in range(limitm[1], limitp[1]+1):
                                for superz in range(limitm[2], limitp[2]+1):
                                    trans = np.array([superx, supery, superz])
                                    for bulk in crystal.get_basislist():
                                        if input[0] == 'o':  # orthonormal basis
                                            output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}  > {:8.4f} {:8.4f} {:8.4f}\n'.format('Bulk{:.0f}'.format(bulk[0]),*change_base(arr=(bulk[1:4]+trans), crystal=crystal, inv=True), *change_base(arr=(bulk[1:4]+trans), crystal=crystal, inv=True)))
                                        elif input[0] == 's':  # supercell basis
                                            output.writelines('{:12s} {:8.4f} {:8.4f} {:8.4f}  > {:8.4f} {:8.4f} {:8.4f}\n'.format('Bulk{:.0f}'.format(bulk[0]),*(bulk[1:4]+trans), *(bulk[1:4]+trans)))


def print_license_notice() -> None:
    """Print license notice at the beginning of each log file"""
    license_notice = "KineCluE - Copyright (C) 2018 CEA, École Nationale Supérieure des Mines de Saint-Étienne\n" \
                     "This program comes with ABSOLUTELY NO WARRANTY; for details, see the license terms.\n" \
                     "This is free software, and you are welcome to redistribute it\n" \
                     "under certain conditions; see the license terms for details.\n" \
                     "You are required to cite the following paper when using KineCluE or part of KineCluE:\n" \
                     "T. Schuler, L. Messina and M. Nastar, Computational Materials Science (2019) [doi: https://doi.org/10.1016.j.commatsci.2019.109191]"
    logger.info(license_notice)


def produce_error_and_quit(error_message: str) -> None:
    logger.info("ERROR! " + error_message)
    raise Exception("ERROR! " + error_message)


def proj(vect: np.ndarray, crystal: Crystal, direc: list) -> list:
    vect = change_base(arr=vect, crystal=crystal, inv=True) # vector in the orthonormal base
    return [np.dot(e, vect) for e in direc]


def rotate_matrix(matrix: np.ndarray, symop: SymOp) -> np.ndarray:
    """Apply a symmetry rotation to a matrix (for instance, for elastic dipole symmetries)"""
    R = symop.get_rotation()
    Rinv = np.linalg.inv(R)
    return np.dot(np.dot(R, matrix), Rinv)


def untranslatedfinal(jumpmech, component_list: list, iniconf, finconf, name_perms: list, fnum=0, init=False):
    ini_position = np.array([iniconf.get_defect_position(def_idx=idx) for idx in range(len(component_list))], dtype=float)
    if (not init) and (iniconf.get_beyond() or finconf.get_beyond()):
        return ini_position, np.array([finconf.get_defect_position(def_idx=idx) for idx in range(len(component_list))], dtype=float)

    jumpvec = []
    specs = [k.get_species().get_index() for k in component_list]
    for cons in jumpmech.get_constraints():
        if cons.get_inispecies().get_index() != 0:
            jumpvec.append([cons.get_inispecies().get_index(), cons.get_finposition() - cons.get_iniposition()])
            for i, s in enumerate(specs):
                if s == cons.get_inispecies().get_index():
                    del specs[i]
                    break
    for s in specs:
        jumpvec.append([s, np.array([0., 0., 0.])])
    jumpvec = sorted(jumpvec, key=lambda x:x[0])

    ljumpvec = [[jumpvec[a] for a in b] for b in name_perms]
    ldiff = [np.array([finconf.get_defect_position(a) for a in b])-ini_position for b in name_perms]
    for diff in ldiff:
        for j in ljumpvec:
            tall = diff - np.array([a[1] for a in j])
            for k in range(len(tall) - 1, -1, -1):
                tall[k] -= tall[0]
            if sum(sum(np.abs(tall))) < 1e-8:
                return ini_position, ini_position+np.array([a[1] for a in j])
    produce_error_and_quit("Could not find the untranslated final position for jump frequency {}".format(fnum))


def vect2deftrans(vec: np.ndarray, all_defect_list: list) -> list:
    """ Takes a vector and transforms it into a defect index and an integer translation vector, written as a 4-item list."""
    trans = np.floor(vec + tol_pbc)
    vec2 = vec - trans  # do not modify vec, otherwise it will mess up the values in the upper-level functions
    found = False
    for defect in all_defect_list:
        if are_equal_arrays(defect.get_sublattice(), vec2):
            return [defect, trans]
    if not found:
        produce_error_and_quit("In vect2deftrans, problem in identifying defect for vector {}. Please check coordinates in the input file.".format(vec))

def config_sym_define(coords_list: list, symop_list: list) -> list:
    config_symops=[]
    symop_relations=[]
    test_com= np.array([0,0,0]) #Find center of mass
    n_coords=len(coords_list)
    #coords_list is a list of arrays, list_coords is a list of lists representing the coordinates
    list_coords=[]
    for i in range(n_coords):
        test_com= np.array(test_com)+np.array(coords_list[i])
        list_coords.append(coords_list[i].tolist())
    test_com= test_com/float(n_coords)
    
    sorted_coords=sortcoords(list_coords)

    n_symops=len(symop_list)
    rotated_coords_lists=apply_symmetry_operations(vector_list=list_coords, symop_list=symop_list, unique=False)
    for i in range(n_symops):
        com_=[0,0,0]
        for j in range(n_coords): #calculate 'center of mass' for rotated coordinates
            com_= np.array(com_)+np.array(rotated_coords_lists[i][j])
        com_=com_/float(n_coords)
        com_=test_com-com_ # now the shift in center of mass between initial and rotated configurations. 
        
        com_flag=True 
        for comp in com_:
            temp=abs(comp) % 1.0 #Need to evaluate separately to prevent error in are_equals()
            if not (are_equals(temp,0) or are_equals(temp,1)): #Only allow symmetries with whole unit cell shifts. How will this play into lattices with multiple basis atoms but same type?
                com_flag=False 
        if not com_flag:
            continue
        
        for j in range(n_coords): #shift coordinates so have same center of mass
            rotated_coords_lists[i][j]=np.array(rotated_coords_lists[i][j])+np.array(com_)
            rotated_coords_lists[i][j]=rotated_coords_lists[i][j].tolist()
            #coordinates now have been rotated and shifted back to same center of mass
        
        test_sorting_coords=copy.copy(rotated_coords_lists[i])
        test_sorting_coords=sorted(test_sorting_coords)
        unique_x=[0]# stores indicies of new x-coordinate values.
        for j in range(1,len(test_sorting_coords)):
            if not (are_equals(test_sorting_coords[j-1][0],test_sorting_coords[j][0])):
                unique_x.append(j)
        for j in range(len(unique_x)-1):
            if ((unique_x[j+1]-unique_x[j])>1): #multiple coordinates with the same x-value
                test_sorting_coords[unique_x[j]:unique_x[j+1]]=sorted(test_sorting_coords[unique_x[j]:unique_x[j+1]],key=sortSecond)

        if (are_equal_arrays(np.array(sorted_coords),np.array(test_sorting_coords))):
            #Symop valid for configuration
            config_symops.append(i)
            relations=[]
            for j in range(n_coords):
                for k in range(n_coords):
                    if are_equal_arrays(np.array(coords_list[j]), np.array(rotated_coords_lists[i][k])):
                        relations.append(k)
                        break
            symop_relations.append(relations)

    return config_symops, symop_relations

def sym_config_add_unique(ini_species_ind:list, config_template: CatConfTemplate, config_dict: dict, config_species_lists: list) -> None:
    symmetry_orders=config_template.get_relations()
    sym_list=[]
    flag= True
    config_index=len(config_dict)
    for i in range(len(symmetry_orders)):
        sym_species=[ini_species_ind[j] for j in symmetry_orders[i]]
        sym_list.append(sym_species)
        if (flag and (test_config(sym_species,config_dict))):
            flag=False
            break

    if flag: #Unique, so add to dictionary and list. 
        config_dict.update({np.array_str(np.array(sym_list[0])):str(config_index)}) #adding first of symmetry variants to dictionary. 
        config_species_lists.append(sym_list[0])


def test_config(test_conf,config_dict):
    return np.array_str(np.array(test_conf)) in config_dict.keys()


def check_configuration_exclusion(species_list: list, config_template: CatConfTemplate) -> bool:
    #Checks if species arrangement in configuration is impossible to size exclusion.
    #Returns True if arrangement is possible. 
    flag= True
    n_sites=config_template.get_n_sites()
    distances=config_template.get_distances()

    for i in range(n_sites):
        r1=species_list[i].get_radius()
        if (r1>0.0): #non-vacant site
            for j in range(i,n_sites):
                r2=species_list[j].get_radius()
                if (r2>0.0 and distances[i][j]>0.0):
                    temp=(r1+r2)-distances[i][j]
                    if (temp>0.0):
                        flag= False 
                if not flag:
                    break
        if not flag:
            break
    return flag

def find_all_events(allevents:list, config_species_lists:list, config_template: CatConfTemplate, config_dict: dict) -> list:
    #Checks list of configuration species lists. Finds all events and returns list of event information formatted: [initial, final, event_type_index].
    full_event_list=[]
    full_event_set=set()
    symmetry_orders=config_template.get_relations()
    coord_dict=config_template.get_coord_dict()
    sym_species_list=[config_species_lists[0] for _ in range(len(symmetry_orders))] # dummy values, overwrite later
    for c, conf in enumerate(config_species_lists):
        for i in range(len(symmetry_orders)):
            sym_species=[conf[j] for j in symmetry_orders[i]]
            sym_species_list[i]=sym_species
        for sym in sym_species_list:
            for event in allevents: #loop over each event
                for cons in event.get_constraints():
                    found_inicon= False 
                    ini_ind=int(coord_dict[np.array_str(np.array(cons.get_iniposition()))])
                    if cons.get_inispecies().get_index() == sym[ini_ind]:
                        found_inicon= True 
                    if not found_inicon:
                        break 
                if found_inicon: #All constraint initial conditions met.
                    #find what the final configuration is
                    finspec=copy.copy(sym)
                    for cons in event.get_constraints():
                        finindex= int(coord_dict[np.array_str(np.array(cons.get_finposition()))])
                        finspec[finindex]=cons.get_finspecies().get_index()
                    tmp_flag=True 
                    for i in range(len(symmetry_orders)):
                        sym_finspec=[finspec[j] for j in symmetry_orders[i]]
                        finspec_str=np.array_str(np.array(sym_finspec))
                        if finspec_str in config_dict:
                            tmp_flag= False 
                            finconfig_ind= int(config_dict[finspec_str])
                            break
                    if tmp_flag:
                        print("PROBLEMS!!!")
                    unique_event=[c, finconfig_ind, event.get_index()]
                    unique_event_str=np.array_str(np.array(unique_event))
                    if not unique_event_str in full_event_set:
                        full_event_set.add(unique_event_str)
                        full_event_list.append(unique_event)
    return full_event_list

def config_event_def(coords_list: list, jump_list: list) -> list:
    #translates each event mechanism and related constraints into a set of constraints defined in terms of configuration site indicies.
    config_coords_dict={}
    for i,cord in enumerate(coords_list):
        config_coords_dict.update({np.array_str(np.array(cord)):str(i)})
    config_jump_list=[]
    for event in jump_list:
        config_cons=[]
        for cons in event.get_constraints():
            ini_ind=int(config_coords_dict[np.array_str(np.array(cons.get_iniposition()))])
            fin_ind=int(config_coords_dict[np.array_str(np.array(cons.get_finposition()))])
            if ini_ind!=fin_ind:
                produce_error_and_quit("Please input constraints as relating to a single coordinate at a time. (Event {})".format(event.get_name()))
            config_cons.append([ini_ind, cons.get_inispecies().get_index(), cons.get_finspecies().get_index()])
        config_jump_list.append(config_cons) 
    return config_jump_list #each entry relates to one event. Each event is a list of constraints defined within the configuration template. Each constraint is the index of of the relevant configuration site, and the initial and final species at that site. 


def find_all_events_config(allevents:list ,config_species_lists: list, config_template: CatConfTemplate, config_dict:dict) -> list:
    full_event_list=[]
    full_event_set=set()
    symmetry_orders=config_template.get_relations()
    coord_dict=config_template.get_coord_dict()
    sym_species_list=[config_species_lists[0] for _ in range(len(symmetry_orders))] # dummy values, overwrite later
    for c in range(len(config_species_lists)):
        for i in range(len(symmetry_orders)):
            sym_species=[config_species_lists[c][j] for j in symmetry_orders[i]]
            sym_species_list[i]=sym_species
        for sym in sym_species_list:
            for e, event in enumerate(allevents): #loop over each event
                found_event_flag=True
                finspec=copy.copy(sym)
                for con in event:
                    if sym[con[0]]==con[1]:
                        finspec[con[0]]=con[2]
                    else:
                        found_event_flag=False 
                        break
                if found_event_flag:
                    finconfig_ind = -1
                    for i in range(len(symmetry_orders)):
                        sym_finspec=[finspec[j] for j in symmetry_orders[i]]
                        finspec_str=np.array_str(np.array(sym_finspec))
                        if finspec_str in config_dict:
                            finconfig_ind= int(config_dict[finspec_str])
                            break
                    if (finconfig_ind == -1):
                        #Constraints for intial configuration met, but no matching final configuration due to size exclusions. Event disallowed.
                        continue
                    unique_event=[c, finconfig_ind, e]
                    unique_event_str=np.array_str(np.array(unique_event))
                    if not unique_event_str in full_event_set:
                        full_event_set.add(unique_event_str)
                        full_event_list.append(unique_event)
    return full_event_list

def sortSecondtwo(vecs:list) -> list:
    #intended to return list of two vectors (in list form) sorted by second entry
    if (are_equals(vecs[0][1], vecs[1][1])): #same second entry, don't sort
        return vecs 
    if (vecs[0][1] < vecs[1][1]):
        return vecs 
    if (vecs[0][1] > vecs[1][1]):
        new_vecs = []
        new_vecs.append(vecs[1])
        new_vecs.append(vecs[0])
        return new_vecs

def sortSecond(val):
    return val[1]

def sortcoords(vecs:list) -> list:
    if len(vecs) <= 1:
        return vecs 
    else:
        new_vecs = sortcoordsinternal(vecs)
        return new_vecs

def sortcoordsinternal(coords:list) -> list:
    #convert coordinates from arrays to lists before passing to this function
    #Should NOT call for a single coordinate, should now be taken care of by outer sortcoords function
    coords=sorted(coords)#sorts by 1st term
    new_coords = []
    temp_array = []
    temp_array.append(coords[0])
    array_len = 1
    for i in range(1, len(coords)):
        if (are_equals(coords[i-1][0], coords[i][0])):
            temp_array.append(coords[i])
            array_len = array_len + 1
        else: 
            if (array_len == 1):
                sort_array = temp_array
            if (array_len == 2):
                sort_array = sortSecondtwo(temp_array)
            if (array_len > 2):
                sort_array = sorted(temp_array, key= sortSecond) #sorts by second term
            for b_int in range(len(sort_array)):
                new_coords.append(sort_array[b_int])
            temp_array = []
            temp_array.append(coords[i])
            array_len = 1
    if (array_len == 1):
        sort_array = temp_array
    if (array_len == 2):
        sort_array = sortSecondtwo(temp_array)
    if (array_len > 2):
        sort_array = sorted(temp_array, key=sortSecond)
    for b_int in range(len(sort_array)):
        new_coords.append(sort_array[b_int])
    return new_coords

def find_event_coords(jump_list:list) -> list: #finds all unique coordinates changed by any process/jump
    event_coords = []
    n_jumps = len(jump_list)
    for event in range(n_jumps):
        cons = jump_list[event].get_constraints()
        for con in range(len(cons)):
            #Assuming that initial and final positions are the same, All events coordinates described in terms of transmutations at those sites. 
            if cons[con].get_inispecies().get_name()!=cons[con].get_finspecies().get_name():
                cord_check = []
                cord_check.append(cons[con].get_iniposition()) 
                cord_check.append(cons[con].get_finposition())

                for k in range(2):
                    found_flag = False
                    for cord in range(len(event_coords)):
                        if are_equal_arrays(event_coords[cord], cord_check[k]):
                           found_flag = True
                           break
                    if not found_flag:
                        event_coords.append(cord_check[k])
    #print("found ",len(event_coords), "event coordinates")
    return event_coords

def index_round(x:float) -> int:
    r_x=round(x)
    if x<0:
        if x<r_x:
            x=int(r_x-1)
        else:
            x=int(r_x)
    else:
        if x>r_x:
            x=int(r_x+1)
        else:
            x=int(r_x)
    return x
