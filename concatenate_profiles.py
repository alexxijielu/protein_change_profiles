'''Script for concatencating multiple protein localization profile files together by protein.
Requires a reference list of all proteins in the screens.
Author: Alex Lu
Email: alexlu@cs.toronto.edu
Last Updated: May 23th, 2018

Copyright (C) 2018 Alex Lu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<https://www.gnu.org/licenses/>.'''

from util import openGeneMatrix
from util import packageGeneMatrix
import numpy as np
import argparse

def sort_proteins (screen_list, reference_list):
    '''Given a list of screens and a list of all proteins analyzed in all screens, sort all screens by
    the list of all proteins (create a NAN vector if the protein isn't in a screen), and then concatencate
    the features of the screens together.'''

    # Open and read the reference list
    with open(reference_list, 'r') as f:
        reference = [line.rstrip('\n') for line in f]

    all_headers = ["PROTEIN"]
    all_matrix = []

    # Iterate through all screens
    for screen in screen_list:
        currmatrix = []
        headers, genelist, genematrix = openGeneMatrix(screen)
        genelist = np.array([x.strip(' ') for x in genelist])

        # Automatically get the screen name (to prefix features with for identification)
        screen_name = screen.split("/")[-1].split(".")[0]

        # Sort the screen by the reference list
        for protein in reference:
            index = np.where(genelist == protein)[0]
            if len(index) == 0:
                vector = np.empty((len(headers[:-1])))
                vector[:] = np.nan
                currmatrix.append(vector)
            else:
                currmatrix.append(genematrix[index][0])

        # Get a concatencated list of all feature names
        for header in headers[1:]:
            all_headers.append(screen_name + "_" + header)

        # Append the sorted matrices together
        currmatrix = np.array(currmatrix)
        if all_matrix == []:
            all_matrix = currmatrix
        else:
            all_matrix = np.hstack((all_matrix, currmatrix))

    return all_matrix, all_headers, reference

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='+', help="List of files to concatencate")
    parser.add_argument("--output", help="Output to write to.", type=str)
    parser.add_argument("--reference", help="Location of list containing all genes in screens.", type=str)
    args = parser.parse_args()

    all_matrix, all_headers, reference = sort_proteins(args.files, args.reference)
    packageGeneMatrix(args.output, all_headers, reference, all_matrix)
