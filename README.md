# Protein Change Profiles
Python code for unsupervised protein localization change detection for microscopy images, as described in Lu et al. 2018 (https://elifesciences.org/articles/31872)

No installation is required, just clone the repository and run the scripts as is.

batch_segmentation.py: Runs the Budding Yeast Morphologist segmentation and feature extraction pipeline on all of the images in a given input directory. Requires software here: https://github.com/lfhandfield/Budding-Yeast-morphologist
- To run: python batch_segmentation.py (directory of tif files to analyze) (location of the bin folder for the Budding Yeast Morphologist software)

average_single_cells.py: Averages single cell features extracted by the Budding Yeast Morphologist software
- To run: python average_single_cells.py (directory of feature files) (output file)

calculate_protein_change_profiles.py: Calculates the protein localization change profiles as described by Lu and Moses 2016 (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0158712)
- To run: python calculate_protein_change_profiles.py (file for untreated wild-type screen) (file for perturbation screen) (output file)

concatenate_profiles.py: Concatenate protein localization change profiles for different perturbations together:
- To run: python concatenate_profiles.py -files (list of protein localization change profile files) -output (output file) -reference (master list of all proteins)

# Dependencies #
- numpy 1.14.2
- sklearn 0.19.1
- scipy 1.0.0

