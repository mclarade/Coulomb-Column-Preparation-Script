from __future__ import division
import os, math, random
import numpy as np
import argparse
import json

from scipy.spatial.distance import pdist, squareform


class Atomic_Data:
    def __init__(self, atom_type, atom_info):
        self.atom_type = atom_type
        self.atom_name = atom_info[0]
        self.charge = atom_info[1]
        self.position_in_array = 0
        self.energies = []
        self.coulomb_column_array = None

    def increment_position(self):
        self.position_in_array += args.randomizations


    def add_energies(self, int_file):
        with open(int_file, 'r') as atomic_file:
            atomic_lines = atomic_file.readlines()
            #make this compatable with the dictonary code
            for line in atomic_lines:
                if line == atomic_lines[6]:
                    #Fix this to work with two character symbols
                    atom_label = line[0]
                if line.startswith('              K'):
                    floatenergy = float(line.split()[3])
            for i in range(0, args.randomizations):
                self.energies.append(floatenergy)

    def initialize_numpy_bins(self, wavefunction_and_files):
        counter = 0
        string_check = str(args.WaveFunctionExtension + self.atom_type)
        for intfile in wavefunction_and_files.values():
            for file in intfile:
                if string_check in str(file):
                    counter += 1
        dimension0 = counter * args.randomizations
        if dimension0 == 0:
            pass
        else:
            self.coulomb_column_array = np.zeros([dimension0, args.length_of_wavefunction])
    
    def save_out_data(self):
        if args.cutoff:
            np.save((self.atom_name + '_Cutoff' + str(args.cutoff) + 'Randomzations' +
                     str(args.randomizations)), self.coulomb_column_array)
            np.save((self.atom_name + '_Cutoff' + str(args.cutoff) + 'Randomzations' +
                     str(args.randomizations)), self.energies)
        else:
            np.save((self.atom_name + 'Randomzations' +
                     str(args.randomizations)), self.coulomb_column_array)
            np.save((self.atom_name + 'Randomzations' +
                     str(args.randomizations))+'Energies', self.energies)


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


def read_in_data(AtomInputList):
    wavefunction_and_file_dict = {}
    InputList = file_list_builder()
    for atomfilename in InputList:
        wavefunction = atomfilename.split(args.WaveFunctionExtension)
        if wavefunction[0] not in wavefunction_and_file_dict:
            wavefunction_and_file_dict[wavefunction[0]] = [atomfilename]
        else:
            atomfilelist = wavefunction_and_file_dict[wavefunction[0]]
            atomfilelist.append(atomfilename)
            wavefunction_and_file_dict[wavefunction[0]] = atomfilelist
    return wavefunction_and_file_dict


def generate_coulomb_matrix(wavefunction, AtomInputList):
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
    if args.cutoff:
        inverse_distance_matrix[inverse_distance_matrix <= 1/args.cutoff] = 0
    labels = np.asarray(label_list)
    charges = []
    for element in labels:
        element = ''.join([i for i in element if not i.isdigit()])
        charges.append(label_to_charge(element, AtomInputList))
    hcharges = np.asarray(charges)
    vcharges = hcharges.reshape(len(hcharges),-1)
    coulomb_in_progress = np.multiply(inverse_distance_matrix,hcharges)
    coulomb_matrix = np.multiply(coulomb_in_progress,vcharges)
    return labels, coulomb_matrix


def label_to_charge(label, AtomInputList):
    if label in AtomInputList.keys():
        return AtomInputList[label][1]


def gen_coulomb_column(matrix, label):
    coulomb_column = []
    label_list = []
    coulomb_list = []
    for element in matrix:
        if element > 1e-7:
            coulomb_column.append(element)
    coulomb_column.sort()
    label_list.append(label[0])
    for atom in args.AtomInputList:
        if label == atom:
            for i in range(0, args.randomizations):
                while (len(coulomb_column) != args.length_of_wavefunction):
                    coulomb_column.append(0)
                coulomb_list.append(coulomb_column)



def shuffle_coulomb_columns(Data_Array):
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
    with open('Atom_Dict.json') as AtomsIn:
        AtomInputList = json.loads(AtomsIn.read())
    Atomic_Data_Dict = {}
    wavefunction_and_file_dict = read_in_data(AtomInputList)
    for Atom_Type in AtomInputList:
        New_Class = Atomic_Data(Atom_Type, AtomInputList[Atom_Type])
        Atomic_Data_Dict[Atom_Type] = New_Class
        New_Class.initialize_numpy_bins(wavefunction_and_file_dict)
    for wavefunction in wavefunction_and_file_dict.keys():
        print wavefunction
        labels, distance_matrix = generate_coulomb_matrix(wavefunction+'.wfn', AtomInputList)
        for i in range(0, len(distance_matrix)):
            gen_coulomb_column(distance_matrix[i], labels[i])
            atom_label = ''.join([j for j in labels[i] if not j.isdigit()])
            for intfile in wavefunction_and_file_dict[wavefunction]:
                #fix this
                Atomic_Data_Dict[atom_label].add_energies(intfile)
    for key in Atomic_Data_Dict.keys():
        Atomic_Data_Dict[key].save_out_data()

                if Atom_Class == atom_label:
                    for intfile in wavefunction_and_file_dict[wavefunction]:
                        Atomic_Data_Dict[Atom_Type].add_energies(intfile)

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
                        help='Set cutoff distance in Bohr, default = None, suggested values are between 2.0 and 11.0 Bohr',
                        type=float)
    parser.add_argument('-r', '--randomizations',
                        dest='randomizations',
                        help='Set number of randomizations to perform when preparing coulomb columns, default = 1',
                        type=int,
                        default=1)
    parser.add_argument('-a', '--atoms',
                        dest='AtomInputList',
                        help='Add to list of atoms to be inspected, takes input in the form Symbol:Name (eg, H:Hydrogen)',
                        type=dict,
                        default={'H': 'Hydrogen', 'C': 'Carbon', 'N': 'Nitgrogen', 'O': 'Oxygen'})
    parser.add_argument('-l', '--length_of_wavefunction',
                        dest='length_of_wavefunction',
                        help='set this value equal to the largest number of atoms in a wavefunction, less one, default = 18',
                        type=int,
                        default=18)
    args = parser.parse_args()
    main(args)
