from __future__ import division
from sklearn import preprocessing
from scipy.spatial.distance import pdist, squareform
import os
import math
import numpy as np
import argparse
import json



class Atomic_Data:
    """
    This class contains all of the information for each atom type being processed
    """
    def __init__(self, atom_type, atom_info):
        self.atom_type = atom_type
        self.atom_name = atom_info[0]
        self.charge = atom_info[1]
        self.position_in_array = 0
        self.energies = []
        self.coulomb_column_array = None

    def add_energies(self, int_file):
        """
        Add_energies opens the appropriate input file, reads in the energy, and appends it 
        to the energy list once for each time the coulomb column is randmized.
        """
        with open(int_file, 'r') as atomic_file:
            atomic_lines = atomic_file.readlines()
            for line in atomic_lines:
                if line.startswith('              K'):
                    floatenergy = float(line.split()[3])
            for _ in range(0, args.randomizations):
                self.energies.append(floatenergy)


    def take_coulomb_column(self, coulomb_column):
        """
        Take_coulomb_column places each coulomb column into the appropriate array once
        for each time it is to be randomized.
        """
        try: 
            for _ in range(0, args.randomizations):
                self.coulomb_column_array[self.position_in_array] = coulomb_column
                self.position_in_array += 1
        except:
            pass


    def initialize_numpy_array(self, wavefunction_and_files):
        """
        Initialize_numpy_array creates an empty array to store coulomb columns
        """
        counter = 0
        string_check = str(args.WaveFunctionExtension + self.atom_type)
        #TODO check functionality
        for intfile in wavefunction_and_files.values():
            for file in intfile:
                if string_check in str(file):
                    counter += 1
        dimension0 = counter * args.randomizations
        if dimension0 == 0:
            pass
        else:
            self.coulomb_column_array = np.zeros([dimension0, args.length_of_wavefunction])


    def trim_zero_columns(self):
        """
        Trim_zero_columns removes all of the columns that contain only zeroes before
        randomization happens.
        """
        try:
            self.coulomb_column_array = self.coulomb_column_array[:, (self.coulomb_column_array == 0).sum(axis=0) != self.coulomb_column_array.shape[0]]
        except AttributeError:
            pass
#            raise AttributeError('Please remove unused atom type from Atom_Dict.json, ' + self.atom_name)


    def shuffle_coulomb_columns(self):
        """
        shuffle_coulomb_columns is intended to go through every row of the output array and shuffle
        each coulomb column
        """
        print "Random Folding Initiated " + self.atom_name
        try: 
            print self.coulomb_column_array.shape[0]
            for i in range(0, self.coulomb_column_array.shape[0]):
                np.random.shuffle(self.coulomb_column_array[i])
        except:
            pass

    def scale_energies(self):
        # if args.scale_columns == True:
        #     preprocessing.scale(self.coulomb_column_array, args.column_norm)
        if type(self.energies) == list:
            self.energies = np.asarray(self.energies)
        self.energies = np.reshape(self.energies, (-1, 1))
        self.energies = preprocessing.scale(self.energies, axis=0)

    def calculate_energy_statistics(self):
        self.energies = np.asarray(self.energies)
        average_energy = np.average(self.energies)
        standard_deviation_energy = np.std(self.energies)
        median = np.median(self.energies)
        minimum = np.min(self.energies)
        maximum = np.max(self.energies)
        with open(self.atom_name+'_statistics', 'w') as statistics:
            statistics.write('Average Energy = ' + str(average_energy) + '\n')
            statistics.write('Standard Deviation = ' + str(standard_deviation_energy) + '\n')
            statistics.write('Median Energy = ' + str(median) + '\n')
            statistics.write('Minimum Energy = ' + str(minimum) + '\n')
            statistics.write('Maximum Energy = ' + str(maximum) + '\n')

    def save_out_data(self):
        """
        save_out_data is designed to save all of the processed data into numpy arrays
        """
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
        print self.atom_name
        print self.coulomb_column_array
        print self.energies


#TODO: roll this in with read_in_data?   
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
    """
    read_in_data sorts every input file into dictionarys based on the root wavefunction in the format
    wavefunction:[input1, input2, ..., inputX]
    """
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
    """
    generate_coulomb_matrix reads in the coordinates from the wavefunction
    file and converts them into a distance matrix, then inverts the distance matrix
    before multiplying it by the relevant charges, then returning the completed
    coulomb matrix and the corresponding atom labels 
    """
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
        #TODO can refer to class, but don't want to touch it right now
        charges.append(label_to_charge(element, AtomInputList))
    hcharges = np.asarray(charges)
    vcharges = hcharges.reshape(len(hcharges),-1)
    coulomb_in_progress = np.multiply(inverse_distance_matrix,hcharges)
    coulomb_matrix = np.multiply(coulomb_in_progress,vcharges)
    return labels, coulomb_matrix


def label_to_charge(label, AtomInputList):
    """
    label_to_charge takes in atom label and returns atom charge
    """
    if label in AtomInputList.keys():
        return AtomInputList[label][1]


def generate_coulomb_column(matrix, label):
    """
    generate_coulomb_column takes in one row of a coulomb matrix
    and removes all of the zero elements, before padding zeroes to the end.
    """
    coulomb_column = []
    for element in matrix:
        if element > 1e-7:
            coulomb_column.append(element)
    coulomb_column.sort()
    while (len(coulomb_column) != args.length_of_wavefunction):
        coulomb_column.append(0)
    return coulomb_column


def main(args):
    with open('Atom_Dict.json') as AtomsIn:
        AtomInputList = json.loads(AtomsIn.read())
    Atomic_Data_Dict = {}
    wavefunction_and_file_dict = read_in_data(AtomInputList)
    for atom_type in AtomInputList:
        New_Class = Atomic_Data(atom_type, AtomInputList[atom_type])
        Atomic_Data_Dict[atom_type] = New_Class
        New_Class.initialize_numpy_array(wavefunction_and_file_dict)
    for wavefunction in wavefunction_and_file_dict.keys():
        labels, coulomb_matrix = generate_coulomb_matrix(wavefunction+'.wfn', AtomInputList)
        if len(coulomb_matrix) > args.length_of_wavefunction:
            raise AttributeError("Length of wavefunction was set smaller than the largest wavefunction")
        for i in range(0, len(coulomb_matrix)):
            atom_label = ''.join([j for j in labels[i] if not j.isdigit()])
            Atomic_Data_Dict[atom_label].take_coulomb_column(
                generate_coulomb_column(coulomb_matrix[i], labels[i]))
        for intfile in wavefunction_and_file_dict[wavefunction]:
            for atom_type in Atomic_Data_Dict.keys():
                if args.WaveFunctionExtension + atom_type in intfile:
                    Atomic_Data_Dict[atom_type].add_energies(intfile)
    for atom_type in Atomic_Data_Dict.keys():
        Atomic_Data_Dict[atom_type].trim_zero_columns()
        if args.statistics == True:
            Atomic_Data_Dict[atom_type].calculate_energy_statistics()
        if args.scale_energies == True:
            Atomic_Data_Dict[atom_type].scale_energies()
        if args.shuffle == True:
            Atomic_Data_Dict[atom_type].shuffle_coulomb_columns()
        Atomic_Data_Dict[atom_type].save_out_data()

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
                        help='set this value equal to the largest number of atoms in a wavefunction, less one, default = 20',
                        type=int,
                        default=20)
    parser.add_argument('--shuffle',
                        dest='shuffle',
                        help='Shuffle coulomb columns before saving, default = True',
                        type=bool,
                        default=True)
    parser.add_argument('--scale_energies',
                        dest='scale_energies',
                        help='scale all energies using SKLearn scale function, default = False',
                        type=bool,
                        default=False)
    parser.add_argument('--calculate_energy_statistics',
                        dest='statistics',
                        help='calculate energy statistics such as mean, stdev, etc, default = True',
                        type=bool,
                        default=True)
    args = parser.parse_args()
    main(args)
