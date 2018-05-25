'''Given a set of single cell features extracted by the image analysis software described in Handfield et al. 2013
(Code: http://www.moseslab.csb.utoronto.ca/louis-f/unsupervised/), select and calculate the truncated average for
the relevant features from the features file.
Input:
A directory containing all of the .txt files for the images you want to analyze. (Only the text files)

Calculates for each image the truncated average of features over 5 bins of cell cycle for each mother/bud cells:
N_Self_DST_Mean, N_Masc_DST_Mean, N_Edge_DST_Mean, N_Cent_DST_Mean, N_Budn_DST_Mean
(i.e. a total of 50 features)

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
<https://www.gnu.org/licenses/>.
'''

import numpy as np
import csv
import os
import argparse
from scipy import stats

def extract_columns(rawData):
    '''Given a full feature matrix, extracts only the relevant features. Extracts:
    CellID, Area, Cell_Type, Rel_Cell_ID, Cell_Prob,
    Intensity_Mean, N_Self_DST_Mean, N_Masc_DST_Mean,
    N_Edge_DST_Mean, N_Cent_DST_Mean, N_Budn_DST_Mean'''
    file = open(rawData)
    list = csv.reader(file, delimiter='\t')
    matrix = np.array([row for row in list])

    features = matrix[:, [1, 13, 21, 38, 42, 71, 75, 79, 83, 87, 91]]
    return features

def assign_bins(features):
    '''Given the 10 filtered features, assign each bud and mother a bin.
    Outputs an 2D array of bins, where index of array corresponds to CellID.'''
    bins = np.zeros((int(features[-1][0]) + 1))
    length = features.shape[0]

    # Bins are based upon bud size, and are calculated in Handfield et al. 2013
    # We merged adjacent bins compared to this paper
    for i in range(1, length):
        if (features[i][2] == 'b'):
            area = float(features[i][1])
            if (area >= 0 and area < 430.43):
                bins[int(features[i][0])] = 1
                bins[int(features[i][4])] = 1
            if (area >= 430.43 and area < 683.26):
                bins[int(features[i][0])] = 2
                bins[int(features[i][4])] = 2
            if (area >= 683.26 and area < 895.33):
                bins[int(features[i][0])] = 3
                bins[int(features[i][4])] = 3
            if (area >= 895.33 and area < 1084.61):
                bins[int(features[i][0])] = 4
                bins[int(features[i][4])] = 4
            if (area >= 1084.61):
                bins[int(features[i][0])] = 5
                bins[int(features[i][4])] = 5

    return bins

def tmean_bins(features, bins):
    '''Return all single-cell data in a bin as a truncated mean (0.05% trim off tails)
    N_Self_DST_Mean, N_Masc_DST_Mean, N_Edge_DST_Mean, N_Cent_DST_Mean, N_Budn_DST_Mean'''
    m0 = []
    m1 = []
    m2 = []
    m3 = []
    m4 = []
    b0 = []
    b1 = []
    b2 = []
    b3 = []
    b4 = []

    # Sort single cells into bins
    for i in range(1, bins.shape[0]):
        if (bins[i] != 0):
            index = np.where(features[:, 0] == str(i))[0][0]
            bin = bins[i]
            type = features[index][2]
            if (type == 'm' and bin == 1):
                m0.append(features[index])
            if (type == 'm' and bin == 2):
                m1.append(features[index])
            if (type == 'm' and bin == 3):
                m2.append(features[index])
            if (type == 'm' and bin == 4):
                m3.append(features[index])
            if (type == 'm' and bin == 5):
                m4.append(features[index])
            if (type == 'b' and bin == 1):
                b0.append(features[index])
            if (type == 'b' and bin == 2):
                b1.append(features[index])
            if (type == 'b' and bin == 3):
                b2.append(features[index])
            if (type == 'b' and bin == 4):
                b3.append(features[index])
            if (type == 'b' and bin == 5):
                b4.append(features[index])
    b0 = np.array(b0)
    b1 = np.array(b1)
    b2 = np.array(b2)
    b3 = np.array(b3)
    b4 = np.array(b4)
    m0 = np.array(m0)
    m1 = np.array(m1)
    m2 = np.array(m2)
    m3 = np.array(m3)
    m4 = np.array(m4)

    tmean_features = []
    try:
        # Go through each feature and append the tmean in order
        for i in range(6, 11):
            tmean_features.append(stats.trim_mean(np.array(b0[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(b1[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(b2[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(b3[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(b4[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(m0[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(m1[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(m2[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(m3[:, i].astype('float')), 0.05))
            tmean_features.append(stats.trim_mean(np.array(m4[:, i].astype('float')), 0.05))
    except:
        # If there are any empty bins (will cause exception), we'll return an empty vector
        tmean_features = []
    return tmean_features

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a directory of single-cell files into a truncated'
                                                 'mean summary file for each protein.')
    parser.add_argument("input", help="Input directory containing files", type=str)
    parser.add_argument("output", help="Output to write to.", type=str)
    args = parser.parse_args()

    inputdir = args.input
    if inputdir[-1] != "/":
        inputdir = inputdir + "/"

    # Open the output file
    output = open(args.output, "w")
    output.truncate()

    # Write the headers
    headers = ["ImageName", "SEF_B0", "SEF_B1", "SEF_B2", "SEF_B3", "SEF_B4", "SEF_M0", "SEF_M1",
               "SEF_M2", "SEF_M3", "SEF_M4", "MCT_B0", "MCT_B1", "MCT_B2", "MCT_B3", "MCT_B4", "MCT_M0",
               "MCT_M1", "MCT_M2", "MCT_M3", "MCT_M4", "EDG_B0", "EDG_B1", "EDG_B2", "EDG_B3", "EDG_B4",
               "EDG_M0", "EDG_M1", "EDG_M2", "EDG_M3", "EDG_M4", "CEN_B0", "CEN_B1", "CEN_B2", "CEN_B3",
               "CEN_B4", "CEN_M0", "CEN_M1", "CEN_M2", "CEN_M3", "CEN_M4", "NEC_B0", "NEC_B1", "NEC_B2",
               "NEC_B3", "NEC_B4", "NEC_M0", "NEC_M1", "NEC_M2", "NEC_M3", "NEC_M4"]
    for header in headers[:-1]:
        output.write(header)
        output.write("\t")
    output.write(headers[-1])
    output.write("\n")

    # Write the features
    for file in os.listdir(inputdir):
        filepath = inputdir + file
        try:
            features = extract_columns(filepath)
            bins = assign_bins(features)
            final_features = tmean_bins(features, bins)

            if final_features == []:
                print ("File skipped: Not enough valid cells in ", file)
            else:
                output.write(file)
                output.write("\t")
                for feature in final_features[:-1]:
                    output.write(str(feature))
                    output.write("\t")
                output.write(str(final_features[-1]))
                output.write("\n")
        except:
            print ("ERROR: Could not process file in directory - skipped ", file)
    output.close()

