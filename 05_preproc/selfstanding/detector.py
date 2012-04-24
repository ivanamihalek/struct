#!/usr/bin/python
import sys
from Bio.PDB.PDBParser import PDBParser

# parameters
pattern_vector_alpha = [5.43, 5.05, 6.20]
pattern_vector_beta = [6.5, 10, 13]
pattern_vector_310 = [5.14, 6, 8.63]
alpha_count_min = 4
beta_count_min = 3

####################################################
# function that prints output to the screen
####################################################

def print_ss(sec_structures):
    for key in sorted(sec_structures.iterkeys()):
        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(sec_structures[key][0], key, sec_structures[key][1], sec_structures[key][2][0], sec_structures[key][2][1], sec_structures[key][2][2])


####################################################
# function that tracks back and change secondary structure information to "C"
# in sec_structure dictionary
####################################################

def trace_back(count, count_min, res_id_old, sec_structures):
    if (count < count_min):
        for i in range(count):
            sec_structures[res_id_old -i][1] = "C"

####################################################
# function that removes isolated beta strands and alpha helicies 
# from sec_structure dictionary in accordance with alpha_count_min i beta_count_min parameters. 
# It changes secondary structure information
# to "C"
####################################################

def clean_short_structures(sec_structures):
    state_old = "S"
    count = 0
    res_id_old = 0
    for res_id in sorted(sec_structures.iterkeys()):
        state = sec_structures[res_id][1]
        if (state == "C"):
            if (state_old == "S" or state_old == "C"):
                pass
            elif (state_old == "H"):
                trace_back(count, alpha_count_min, res_id_old, sec_structures)
            elif (state_old == "E"):
                trace_back(count, beta_count_min, res_id_old, sec_structures)
            count = 0
            state_old = "C"
            
        elif (state == "H"):
            if (state_old == "S" or state_old == "C"):
                count = 1
            elif (state_old == "E"):
                trace_back(count, beta_count_min, res_id_old, sec_structures)
            elif (state_old == "H"):
                if (res_id_old + 1 == res_id):
                    count = count +1
                else:
                    trace_back(count, alpha_count_min, res_id_old, sec_structures)
                    count  = 0
            res_id_old = res_id
            state_old = "H"
            
        elif (state == "E"):
            if (state_old == "S" or state_old == "C"):
                count = 1
            elif (state_old == "H"):
                trace_back(count, alpha_count_min, res_id_old, sec_structures)
            elif (state_old == "E"):
                if (res_id_old + 1 == res_id):
                    count = count +1
                else:
                    trace_back(count, beta_count_min, res_id_old, sec_structures)
                    count  = 0
            res_id_old = res_id
            state_old = "E"

####################################################
# function that finds beta strands and in accordance 
# chain the secondary structure values for particular residue.
# Values are stored in sec_structures dictionary
####################################################

def find_beta(distances, sec_structures):
    threshold = 2
    for rec in distances:
        if (rec[1] == None):
            continue
        if distance(rec[1:len(rec)+1], pattern_vector_beta) < threshold:
            sec_structures[rec[0]][1] = "E"

####################################################
# function that finds alpha helicies and in accordance 
# chain the secondary structure values for particular residue.
# Values are stored in sec_structures dictionary
####################################################

def find_alpha(distances, sec_structures):
    threshold = 1.5
    for rec in distances:
        if (rec[1] == None):
            continue
        if distance(rec[1:len(rec)+1], pattern_vector_alpha) < threshold:
            sec_structures[rec[0]][1] = "H"

####################################################
# function that finds 310 helicies and in accordance 
# chain the secondary structure values for particular residue.
# Values are stored in sec_structures dictionary
####################################################

def find_310(distances, sec_structures):
    threshold = 1.8
    for rec in distances:
        if (rec[1] == None):
            continue
        if distance(rec[1:len(rec)+1], pattern_vector_310) < threshold:
            sec_structures[rec[0]][1] = "H"



####################################################
# function that returns a dictionary of the secundary structures
# - key - residue number
# - value a list of (chain_id, secondary_structure, coordinates)
####################################################

def find_ss(distances, coord_list):
    sec_structures = {}
    for i in range(len(coord_list)):
        res_id = coord_list[i].get_full_id()[3][1]
        chain = coord_list[i].get_full_id()[2]
        rec = [chain, "C", distances[i][1:4]]
        sec_structures[res_id] = rec
        
    find_alpha(distances, sec_structures)
    find_beta(distances, sec_structures)
#    find_310(distances, sec_structures)
    clean_short_structures(sec_structures)
    return sec_structures

####################################################
# function that returns eucledian distance between two vectors
####################################################

def distance(vector1, vector2):
    dist = sum((x - y) ** 2 for (x, y) in zip(vector1, vector2)) ** 0.5
    return dist

####################################################
# function that returns a list of coordinates (a list with x,y,z) of atoms
####################################################

def determine_ca_coordinates(chain):
    coord_list = []
    for residue in chain.get_list():
        if residue.has_id("CA"):
            ca = residue["CA"]
            coord_list.append(ca)
    return coord_list

####################################################
# function that returns a list of ca atom objects
####################################################

def determine_ca_distances(coord_list):
    distances = []
    for i in range(len(coord_list)):
        res_id = coord_list[i].get_full_id()[3][1]
        if (i + 4 > len(coord_list) -1): # last element in list
            rec = res_id, None, None, None
            distances.append(rec)
            continue
        res_id_plus_4 = coord_list[i+4].get_full_id()[3][1]
        coord = coord_list[i].get_coord()
        if ((res_id + 4) != res_id_plus_4):
            rec = res_id, None, None, None
        else:
            
            coord_plus_2 = coord_list[i+2].get_coord()
            coord_plus_3 = coord_list[i+3].get_coord()
            coord_plus_4 = coord_list[i+4].get_coord()

            dist2 = distance(coord, coord_plus_2)
            dist3 = distance(coord, coord_plus_3)
            dist4 = distance(coord, coord_plus_4)
            
            rec = res_id, dist2, dist3, dist4
        distances.append(rec)
    return distances

####################################################
# function that calculates the secondary structure for each residue by chain
####################################################

def calculate_ss(model):
    for chain in model.get_list():
        coord_list = determine_ca_coordinates(chain)
        distances = determine_ca_distances(coord_list)
        sec_struct = find_ss(distances, coord_list)
        print_ss(sec_struct)
    


def main():
	if len(sys.argv) < 2:
		sys.exit('Usage: %s input_pdb_file' % sys.argv[0])
	pdb_name = sys.argv[1]
	
	parser=PDBParser(PERMISSIVE=1)
	structure_id = "temp"
	structure = parser.get_structure(structure_id, pdb_name)
	model = structure[0]
	
	calculate_ss(model)
	
	
if __name__ == "__main__":
    main()
