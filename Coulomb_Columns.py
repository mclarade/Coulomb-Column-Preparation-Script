from __future__ import division
import os, math, random
import numpy as np
import argparse

from scipy.spatial.distance import pdist, squareform


def file_list_builder():
    '''
    This function takes in the names of all of the files to be used and stores them in lists
    '''
    InputList = []
    for file in os.listdir(os.curdir):
        if file.endswith(args.InputExtension):
            InputList.append(file)
    InputList.sort()
    return InputList


def initialize_numpy_bins():
    '''Creates numpy bins for each atom type'''
    counter_dict = {}
    for atom_type in args.AtomInputList.keys():
        counter_dict.setdefault(atom_type, 0)
    wavefunction_and_file_dict = {}
    InputList = file_list_builder()
    for atomfilename in InputList:
        wavefunction = atomfilename.split(args.WaveFunctionExtension)
        if wavefunction[0] not in wavefunction_and_file_dict:
            wavefunction_and_file_dict[wavefunction[0]] = [atomfilename]
        else:
            #This may be redundant
            atomfilelist = wavefunction_and_file_dict[wavefunction[0]]
            atomfilelist.append(atomfilename)
            wavefunction_and_file_dict[wavefunction[0]] = atomfilelist
        for atom_type in args.AtomInputList.keys():
            string_check = args.WaveFunctionExtension + atom_type
            if string_check in atomfilename:
                counter_dict[atom_type] += 1
                break
    array_dict = {}
    for atom_type in args.AtomInputList.keys():
        if counter_dict[atom_type] == 0:
            pass
        else:
            #look here for output dimensions
            dimension0 = counter_dict[atom_type]
            array_dict[atom_type] = np.zeros((dimension0, args.size_of_wavefunction))
    keylist = wavefunction_and_file_dict.values()
    keylist.sort()
    return keylist, array_dict


def retrieve_coordinates(wavefunction):
    with open(wavefunction, 'r') as waveinput:
        wavelines = waveinput.readlines()
        x_list = []
        y_list = []
        z_list = []
        label_list = []
        for line in wavelines:
            linesplit = line.split()
            if (line == wavelines[0]) or (line == wavelines[1]):
                pass
            elif line.startswith("CENTRE ASSIGNMENTS"): 
                break
            else:
                label_list.append((str(linesplit[0]) + str(linesplit[1])))
                x_list.append(float(linesplit[4]))
                y_list.append(float(linesplit[5]))
                z_list.append(float(linesplit[6]))
    position_matrix = np.stack((x_list,y_list,z_list),axis=-1)
    distance_list = pdist(position_matrix)
    distance_matrix = squareform(distance_list)
    distance_and_identity_matrix = distance_matrix  +  np.identity(len(distance_matrix))
    inverse_distance_and_identity_matrix = 1/distance_and_identity_matrix
    inverse_distance_matrix = inverse_distance_and_identity_matrix - np.identity(len(distance_matrix))
    try:
        inverse_distance_matrix[inverse_distance_matrix <= 1/args.cutoff] = 0
    except ValueError:
        pass
    labels = np.asarray(label_list)
    charges = []
    for element in labels:
        charges.append(label_to_charge(element))
    hcharges = np.asarray(charges)
    vcharges = hcharges.reshape(len(hcharges),-1)
    coulomb_in_progress = np.multiply(inverse_distance_matrix,hcharges)
    coulomb_matrix = np.multiply(coulomb_in_progress,vcharges)
    return labels, coulomb_matrix

def label_to_charge(label):
    atom_label = label[0]
    if atom_label == 'H':
        return 1
    if atom_label == 'C':
        return 6
    if atom_label == 'O':
        return 8

def gen_energy_list(int_file):
    with open(int_file,'r') as atomic_file:
        atomic_lines = atomic_file.readlines()
        for line in atomic_lines:
            if line == atomic_lines[6]:
                atom_label = line[0]
            if line.startswith('              K'):
                floatenergy = float(line.split()[3])
    # if atom_label == 'C':
    #     for i in range(0, args.randomizations):
    #         CEnergyList.append(floatenergy)
    # if atom_label == 'H':
    #     for i in range(0, args.randomizations):
    #         HEnergyList.append(floatenergy)
    # if atom_label == 'O':
    #     for i in range(0, args.randomizations):
    #         OEnergyList.append(floatenergy)

def gen_coulomb_column(matrix, labels):
    for i in range(0, len(matrix)):
        #global CarbonWritten, OxygenWritten, HydrogenWritten
        coulomb_column = []
        label_list = []
        for element in matrix[i]:
            if element > 1e-7:
                coulomb_column.append(element)
        coulomb_column.sort()
        label_list.append(labels[i][0])
        for atom in args.AtomInputList:
            if labels[i] == atom:
                for i in range(0, args.randomizations):
                    while (len(coulomb_column) != args.size_of_wavefunction):
                        coulomb_column.append(0)
                        Array[i] = coulomb_column
        # if atom_label == 'H':
        #     for i in range(0, args.randomizations):
        #         # random.shuffle(coulomb_column)
        #         while (len(coulomb_column) != args.size_of_wavefunction):
        #             coulomb_column.append(0)
        #         HydrogenArray[HydrogenWritten] = coulomb_column
        #         HydrogenWritten += 1
        # if atom_label == 'O':
        #     for i in range(0, args.randomizations):
        #         # random.shuffle(coulomb_column)
        #         while (len(coulomb_column) != args.size_of_wavefunction):
        #             coulomb_column.append(0)
        #         OxygenArray[OxygenWritten] = coulomb_column
        #         OxygenWritten += 1
    


def generate_random_mutations(Data_Array):
    print "Random Folding Initiated"
    counter = 0
    Copy_Array = np.copy(Data_Array)
    Dimension0 = int(Copy_Array.shape[0])
    Dimension1 = int(Copy_Array.shape[1])
    Data_Array = np.zeros((Dimension0,Dimension1))
    for line in Copy_Array:
        copy_list = line.tolist()
        np.random.shuffle(copy_list)
        Data_Array[(counter)] = copy_list
        counter += 1
    return Data_Array


def trim_zero_columns(Array):
    Array = Array[:, (Array == 0).sum(axis=0) != Array.shape[0]]
    return Array

def main(args):
    filelist = os.listdir(os.curdir)
    filelist.sort()
    keylist, array_dict_sorted_by_atom = initialize_numpy_bins()
    # HydrogenArray = np.zeros([HydrogenCounter, args.size_of_wavefunction])
    # CarbonArray = np.zeros([CarbonCounter, args.size_of_wavefunction])
    # OxygenArray = np.zeros([OxygenCounter, args.size_of_wavefunction])

    # HydrogenWritten = 0
    # CarbonWritten = 0
    # OxygenWritten = 0
    for wavefunction in keylist:
        print wavefunction
        labels, distance_matrix = retrieve_coordinates(wavefunction)
        gen_coulomb_column(distance_matrix, labels)
        for intfile in master_dict[wavefunction]:
            gen_energy_list(intfile)
        

    # HydrogenEnergyOut = np.asarray(HEnergyList)
    # CarbonEnergyOut = np.asarray(CEnergyList)
    # OxygenEnergyOut = np.asarray(OEnergyList)

    # HydrogenArray = trim_zero_columns(HydrogenArray)
    # CarbonArray = trim_zero_columns(CarbonArray)
    # OxygenArray = trim_zero_columns(OxygenArray)

    # HydrogenArray = generate_random_mutations(HydrogenArray)
    # CarbonArray = generate_random_mutations(CarbonArray)
    # OxygenArray = generate_random_mutations(OxygenArray)

    # print "Hydrogen", HydrogenArray.shape, HydrogenEnergyOut.shape
    # print "Carbon", CarbonArray.shape, CarbonEnergyOut.shape
    # print "Oxygen", OxygenArray.shape, OxygenEnergyOut.shape

    # np.save(('H_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), HydrogenArray)
    # np.save(('C_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), CarbonArray)
    # np.save(('O_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), OxygenArray)
    # np.save(('Energy_H_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), HydrogenEnergyOut)
    # np.save(('Energy_C_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), CarbonEnergyOut)
    # np.save(('Energy_O_Out_Cutoff' + str(args.cutoff) + 'Randomzations' + str(args.randomizations)), OxygenEnergyOut)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script automates the conversion of wfn and int files into coulomb columns for use with neural network inputs. Note that each wavefunction must have a unique name, or behavior may become unpredictable')
    parser.add_argument('-x', '--extension',
                        dest='InputExtension',
                        help='Extension of aimpac output scripts',
                        default='.int')
    parser.add_argument('-w', '--wavefunction',
                        dest='WaveFunctionExtension',
                        help='Extension of wavefunction files',
                        default='.wfn')
    parser.add_argument('-c', '--cutoff',
                        dest='cutoff',
                        help='Set cutoff distance in Bohr, default = none, suggested values are between 2.0 and 11.0 Bohr',
                        type=float)
    parser.add_argument('-r', '--randomizations',
                        dest='lambda_value',
                        help='Set number of randomizations to perform when preparing coulomb columns, default = 1',
                        type=int,
                        default=-1)
    parser.add_argument('-a', '--atoms',
                        dest='AtomInputList',
                        help='Add to list of atoms to be inspected, takes input in the form Symbol:Name (eg, H:Hydrogen)',
                        type=dict,
                        default={'H': 'Hydrogen', 'C': 'Carbon', 'N': 'Nitgrogen', 'O': 'Oxygen'})
    parser.add_argument('-s', '--size_of_wavefunction',
                        dest='size_of_wavefunction',
                        help='set this value equal to the largest number of atoms in a wavefunction, less one, default = 18',
                        type=int,
                        default=18)
    args = parser.parse_args()
    main(args)
