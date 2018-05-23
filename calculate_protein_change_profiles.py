'''Given the profiled features for two screens, calculate the protein change profiles between the two screens.'''

from util import openGeneMatrix
from util import packageGeneMatrix
import argparse
import numpy as np
import sklearn.metrics.pairwise as skdist

def filter_matrices (reference, condition):
    '''Preprocessing operation - calculates the intersection of proteins between the reference and the
    condition, and sorts them so that they're in the same order.'''
    ref_headers, ref_genelist, ref_genematrix = openGeneMatrix(reference)
    cond_headers, cond_genelist, cond_genematrix = openGeneMatrix(condition)

    sorted_ref = []
    sorted_cond = []
    # Get the intersection of the list
    # Sometimes there can be duplicate proteins, so we'll just take the first occurrence if there is
    intersect = np.intersect1d(ref_genelist, cond_genelist)
    for protein in intersect:
        ref_index = np.where(ref_genelist == protein)[0][0]
        cond_index = np.where(cond_genelist == protein)[0][0]
        sorted_ref.append(ref_genematrix[ref_index])
        sorted_cond.append(cond_genematrix[cond_index])

    sorted_ref = np.array(sorted_ref)
    sorted_cond = np.array(sorted_cond)

    return intersect, sorted_ref, sorted_cond

def subtract_matrices (sorted_ref, sorted_cond):
    '''Produce a subtracted matrix given two sorted matrices'''
    return np.subtract(sorted_ref, sorted_cond)

def modWeights(k, geneMatrix, distMatrix, metric='euclidean'):
    '''Generates a mean and MAD vector for each protein in a protein feature matrix using the k closest genes using euclidean distance
    Inputs: k, wild-type protein feature matrix, change matrix
    Output: mean and MAD of k NN of proteins'''
    # Specify distance metric and get nearest neighbors
    dist = skdist.pairwise_distances(distMatrix, metric=metric)
    nearest = np.argsort(dist, axis=1)[:, 1:(k + 1)]

    # Calculate mean for each protein
    means = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        for index in nearest[gene]:
            for feature in range(0, len(geneMatrix[index])):
                means[gene][feature] += geneMatrix[index][feature]
        means[gene] = means[gene] / k
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated means for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))

    # Calculate medians and MAD for each protein
    medians = np.zeros(geneMatrix.shape)
    MAD = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        neighbors = []
        for index in nearest[gene]:
            neighbors.append(geneMatrix[index])

        neighbors = np.array(neighbors)
        for feature in range(0, len(neighbors[0])):
            medians[gene][feature] = np.median(neighbors[:, feature])
            MAD[gene][feature] = np.median(abs(medians[gene][feature] - neighbors[:, feature]))
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated medians for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))

    return medians, MAD

def calculateModZScores(geneMatrix, means, MAD):
    '''Calculates modified z-score vectors for each gene given mean and variance vectors of kNN neighbors
    Inputs: gene matrix and corresponding mean and variances of kNN neighbors
    Output: zscores'''
    zscores = np.zeros(geneMatrix.shape)
    for gene in range(0, geneMatrix.shape[0]):
        for feature in range(0, len(geneMatrix[gene])):
            zscores[gene][feature] = 0.6745 * (geneMatrix[gene][feature] - means[gene][feature]) / (MAD[gene][feature])
    return zscores

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create protein localization change profiles for a pair of files.'
                                                 'Calculates z-scores for each feature in each protein.')
    parser.add_argument("reference", help="Reference wild-type screen", type=str)
    parser.add_argument("condition", help="Condition screen", type=str)
    parser.add_argument("output", help="Output to write to.", type=str)
    args = parser.parse_args()

    genelist, sorted_ref, sorted_cond = filter_matrices(args.reference, args.condition)
    subtracted = subtract_matrices(sorted_ref, sorted_cond)
    means, variances = modWeights(50, subtracted, sorted_ref)
    zscores = calculateModZScores(subtracted, means, variances)

    headers, _, _ = openGeneMatrix(args.reference)
    packageGeneMatrix(args.output, headers, genelist, zscores)