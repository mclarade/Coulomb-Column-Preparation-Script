from __future__ import division, print_function

import argparse
import json
import os

import numpy as np
from scipy.spatial.distance import pdist, squareform


class AtomicData:
    """
    This class contains all of the information for each atom type being processed

    Args:
        atom_type (str): the atom's symbol
        atom_info (tuple): a tuple (name, charge) containing the atom's info
    """
    def __init__(self, atom_type, atom_info, args):
        self.atom_type = atom_type
        self.atom_name = atom_info[0]
        self.charge = atom_info[1]
        self.position_in_array = 0
        self.energies = []
        self.coulomb_column_array = None
        self.args = args

    def add_energies(self, int_file):
        """
        Opens the appropriate input file, reads in the energy, and appends it
        to the energy list once for each time the coulomb column is randomized.

        Args:
            int_file (str): the name of the .int file
        """
        with open(int_file, 'r') as atomic_file:
            for line in atomic_file:
                if line.startswith('              K'):
                    energy = float(line.split()[3])
            for _ in range(0, self.args.n_randomizations):
                self.energies.append(energy)

    def take_coulomb_column(self, coulomb_column):
        """
        Places each coulomb column into the appropriate array once
        for each time it is to be randomized.

        Args:
            coulomb_column (np.ndarray): a Coulomb columb array to append
        """
        for _ in range(0, self.args.n_randomizations):
            self.coulomb_column_array[self.position_in_array] = coulomb_column
            self.position_in_array += 1

    def initialize_numpy_array(self, wavefunction_and_files):
        """
        Creates an empty array to store coulomb columns

        Args:
            wavefunction_and_files (dict): the dictionary containing
                .int files to process. The key is the basename.
        """
        counter = 0
        string_check = str(self.args.wavefunction_extension + self.atom_type)
        #TODO check functionality
        for intfile in wavefunction_and_files.values():
            for file_ in intfile:
                if string_check in str(file_):
                    counter += 1
        dimension0 = counter * self.args.n_randomizations
        if dimension0 > 0:
            self.coulomb_column_array = np.zeros([dimension0, self.args.wavefunction_length])

    def trim_zero_columns(self):
        """
        Removes all of the columns that contain only zeroes before
        randomization happens.
        """
        try:
            mask = (self.coulomb_column_array == 0).sum(axis=0)
            mask = mask != self.coulomb_column_array.shape[0]
            self.coulomb_column_array = self.coulomb_column_array[:, mask]
        except AttributeError:
            raise AttributeError('Please remove unused atom type from Atom_Dict.json, ' + self.atom_name)

    def shuffle_coulomb_columns(self):
        """
        Goes through every row of the output array and shuffles each coulomb column
        """
        print("Random Folding Initiated " + self.atom_name)
        print(self.coulomb_column_array.shape[0])
        for i in range(0, self.coulomb_column_array.shape[0]):
            np.random.shuffle(self.coulomb_column_array[i])

    def save_out_data(self):
        """
        Saves all of the processed data into numpy arrays
        """
        if self.args.cutoff:
            filename = (
                '{atom_name}_Cutoff{cutoff}Randomizations{n_randomizations}'.format(
                    atom_name=self.atom_name,
                    cutoff=self.args.cutoff,
                    n_randomizations=self.args.n_randomizations
                )
            )
            np.save(filename, self.coulomb_column_array)
            np.save(filename + 'Energies', self.energies)
        else:
            filename = (
                '{atom_name}Randomizations{n_randomizations}'.format(
                    atom_name=self.atom_name,
                    cutoff=self.args.cutoff,
                    n_randomizations=self.args.n_randomizations
                )
            )
            np.save(filename, self.coulomb_column_array)
            np.save(filename + 'Energies', self.energies)
        print(self.atom_name)
        print(self.coulomb_column_array)
        print(self.energies)


# TODO: roll this in with read_in_data?
def file_list_builder(input_extension):
    """
    Takes in the names of all of the files to be used and stores them in lists

    Args:
        input_extension (str): the extension to use when finding input files

    Returns:
        list: the list of files to be used
    """
    input_files = []
    for file_ in os.listdir(os.curdir):
        if file_.endswith(input_extension):
            input_files.append(file_)
    input_files.sort()
    return input_files


def read_in_data(input_extension, wavefunction_extension):
    """
    Sorts every input file into dictionarys based on the root wavefunction in the format
    wavefunction:[input1, input2, ..., inputX]

    Args:
        input_extension (str): the extension to use when finding input
            files
        wavefunction_extension (str): the extension to use when finding
            wavefunction files

    Returns:
        dict: a dictionary containing the wavefunction filenames, indexed by
            the basename of the wavefunction files
    """
    wavefunction_and_files = {}
    input_list = file_list_builder(input_extension)
    for atom_filename in input_list:
        wavefunction = atom_filename.split(wavefunction_extension)
        atom_files = wavefunction_and_files.setdefault(wavefunction[0], [])
        atom_files.append(atom_filename)
    return wavefunction_and_files


def generate_coulomb_matrix(filename, atom_input_list, cutoff=None):
    """
    Reads in the coordinates from the wavefunction file and builds
    the Coulomb matrix

    Args:
        filename (str): the wavefunction file
        atom_input_list (list): TODO: I'm not sure what this is yet
        cutoff (float, optional): the cutoff value

    Returns:
        *: the atom labels
        np.ndarray: the Coulomb matrix
    """
    with open(filename, 'r') as waveinput:
        x_list = []
        y_list = []
        z_list = []
        label_list = []
        for i, line in enumerate(waveinput):
            # TODO: A comment with a sample line of the input file would be
            #       helpful here
            linesplit = line.split()
            if i < 2:
                continue
            elif line.startswith("CENTRE ASSIGNMENTS"):
                break
            else:
                label_list.append((str(linesplit[0]) + str(linesplit[1])))
                x_list.append(float(linesplit[4]))
                y_list.append(float(linesplit[5]))
                z_list.append(float(linesplit[6]))

    position_matrix = np.stack((x_list, y_list, z_list), axis=-1)
    distance_matrix = squareform(pdist(position_matrix))
    distance_and_identity_matrix = distance_matrix + np.identity(len(distance_matrix))
    inverse_distance_and_identity_matrix = 1 / distance_and_identity_matrix
    inverse_distance_matrix = inverse_distance_and_identity_matrix - np.identity(len(distance_matrix))
    if cutoff:
        inverse_distance_matrix[distance_matrix >= cutoff] = 0
    labels = np.asarray(label_list)
    charges = []
    for element in labels:
        element = ''.join([char for char in element if not char.isdigit()])
        # TODO: can refer to class, but don't want to touch it right now
        charges.append(label_to_charge(element, atom_input_list))

    h_charges = np.asarray(charges)
    v_charges = h_charges.reshape(len(h_charges), -1)
    coulomb_matrix = inverse_distance_matrix * h_charges * v_charges
    return labels, coulomb_matrix


def label_to_charge(label, atom_input_list):
    """
    Returns an atom charge from its label

    Args:
        label (str): the atom symbol
        atom_input_list (dict): a dictionary mapping atom symbol to atomic number
    """
    if label in atom_input_list:
        return atom_input_list[label][1]


def generate_coulomb_column(row, max_length):
    """
    Takes in one row of a coulomb matrix and removes all of the zero elements and
    pads zeroes to the end.

    Args:
        row (np.ndarray): the row of the Coulomb matrix
        max_length (int): the desired length of all Coulomb columns
    """
    coulomb_column = []
    for element in row:
        if element > 1e-7:
            coulomb_column.append(element)
    coulomb_column.sort()
    while len(coulomb_column) != max_length:
        coulomb_column.append(0)
    return coulomb_column


def main(args):
    with open('Atom_Dict.json') as atoms_file:
        atom_input_list = json.load(atoms_file)

    atomic_data_dict = {}
    wavefunction_and_file_dict = read_in_data(args.input_extension, args.wavefunction_extension)
    for atom_type in atom_input_list:
        atomic_data = AtomicData(
            atom_type,
            atom_input_list[atom_type],
            args
        )
        atomic_data_dict[atom_type] = atomic_data
        atomic_data.initialize_numpy_array(wavefunction_and_file_dict)

    for wavefunction in wavefunction_and_file_dict.keys():
        labels, coulomb_matrix = generate_coulomb_matrix(wavefunction + '.wfn',
                                                         atom_input_list,
                                                         cutoff=args.cutoff)
        if len(coulomb_matrix) > args.wavefunction_length:
            raise AttributeError("Length of wavefunction was set smaller than "
                                 "the largest wavefunction")
        for i in range(0, len(coulomb_matrix)):
            atom_label = ''.join([j for j in labels[i] if not j.isdigit()])
            atomic_data_dict[atom_label].take_coulomb_column(
                generate_coulomb_column(coulomb_matrix[i],
                                        max_length=args.wavefunction_length)
            )

        for intfile in wavefunction_and_file_dict[wavefunction]:
            for atom_type in atomic_data_dict:
                if args.wavefunction_extension + atom_type in intfile:
                    atomic_data_dict[atom_type].add_energies(intfile)

    for atom_type in atomic_data_dict:
        atomic_data_dict[atom_type].trim_zero_columns()
        atomic_data_dict[atom_type].shuffle_coulomb_columns()
        atomic_data_dict[atom_type].save_out_data()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script automates the conversion of wfn and int '
                    'files into coulomb columns for use with neural network '
                    'inputs. Note that each wavefunction must have a unique '
                    'name, or behavior may become unpredictable'
    )
    parser.add_argument('-x', '--extension',
                        dest='input_extension',
                        help='Extension of aimpac output scripts',
                        default='.int')
    parser.add_argument('-w', '--wavefunction',
                        dest='wavefunction_extension',
                        help='Extension of wavefunction files',
                        default='.wfn')
    parser.add_argument('-c', '--cutoff',
                        dest='cutoff',
                        help='Set cutoff distance in Bohr, default = None, suggested values are between 2.0 and 11.0 Bohr',
                        type=float)
    parser.add_argument('-r', '--randomizations',
                        dest='n_randomizations',
                        help='Set number of randomizations to perform when preparing coulomb columns, default = 1',
                        type=int,
                        default=1)
    parser.add_argument('-a', '--atoms',
                        dest='atom_input_list',
                        help='Add to list of atoms to be inspected, takes input in the form Symbol:Name (eg, H:Hydrogen)',
                        type=dict,
                        default={'H': 'Hydrogen', 'C': 'Carbon', 'N': 'Nitgrogen', 'O': 'Oxygen'})
    parser.add_argument('-l', '--wavefunction-length',
                        dest='wavefunction_length',
                        help='set this value equal to the largest number of atoms in a wavefunction, less one, default = 18',
                        type=int,
                        default=18)

    main(parser.parse_args())
