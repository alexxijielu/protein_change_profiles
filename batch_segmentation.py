import os
import subprocess
import shlex
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run batch segmentation and feature extraction on a directory of '
                                                 'tif files of yeast microscopy images.')
    parser.add_argument("directory", help="Input directory containing tif files", type=str)
    parser.add_argument("bin", help="Location of Budding Yeast Morphologist Programs.", type=str)
    args = parser.parse_args()

    # Standardize the extensions of the paths
    inputdir = args.directory
    if inputdir[-1] != "/":
        inputdir = inputdir + "/"

    PMbin = args.bin
    if PMbin[-1] != "/":
        PMbin = PMbin + "/"

    # Compile quality measures into a confidence matrix binary file
    # These quality measures are based upon manually curated cells, and used to estimate the probability of objects
    # being cells - the file can be found in the Budding Yeast Morphologist repository
    # (https://github.com/lfhandfield/Budding-Yeast-morphologist) under /example/Quality_measures.txt
    command = PMbin + "PMMakeConfidenceMatrices ./Conf_Matrices.scp " + \
              PMbin.rsplit("/", 2)[0] + "/example/Quality_measures.txt"
    subprocess.call(shlex.split(command))

    images = os.listdir(inputdir)
    for image in images:
        if image.endswith(".tif"):
            ID = image.split(".")[0]

            # Parse tiff files for red channels
            command = PMbin + "PMTiffManip -f 2,4,6,8 " + inputdir + image + " " + inputdir + ID + "_red.tiff"
            subprocess.call(shlex.split(command))

            # Seperate Background and Foreground in Images
            command = PMbin + "PMSegmentation -B 1.0f -b " + inputdir + ID + ".scp -d " + inputdir + ID + "_dist.tiff " + \
                      inputdir + ID + "_red.tiff " + inputdir + ID + "_seg.tiff"
            subprocess.call(shlex.split(command))

            # Identify cell centers
            command = PMbin + "PMFindMultiCover -o " + inputdir + ID + "_ellfit.mcv -R2 1.0 6.0 " + \
                      inputdir + ID + "_red.tiff " + inputdir + ID + "_seg.tiff " + inputdir + ID + "_dist.tiff "
            subprocess.call(shlex.split(command))

            # Find cell sub-segments
            command = PMbin + "PMWaterShed -M -b 1.0 " + inputdir + ID + "_red.tiff " + inputdir + ID + "_watershed.tiff"
            subprocess.call(shlex.split(command))

            # Cell identification from circle coordinated directed sub-segment agglomeration
            command = PMbin + "PMHiddenMapDirect -M -G " + inputdir + ID + "_watershed.tiff " + "-C " +\
                      inputdir + ID + "_cellseg.tiff " + inputdir + ID + "_seg.tiff " + inputdir + ID + \
                      "_dist.tiff " + inputdir + ID + "_ellfit.mcv "
            subprocess.call(shlex.split(command))

            # Parse tiff files for "green" channel
            command = PMbin + "PMTiffManip -f 1,3,5,7 " + inputdir + image + " " + inputdir + ID + "_gre.tiff"
            subprocess.call(shlex.split(command))

            # Copy temporary version of confidence measure
            command = "cp ./Conf_Matrices.scp " + inputdir + ID + "_conf_matrices.scp"
            subprocess.call(shlex.split(command))

            # Measure GFP Intensity and spatial spread
            command = PMbin + "PMExtractFeatures -c " + inputdir + ID + "_conf_matrices.scp -m " + \
                      inputdir + ID + "_ellfit.mcv -a 0.1 10 0.9 0.1 50 -t " + inputdir + ID + ".txt -S " + \
                      inputdir + ID + "_seg.tiff " + inputdir + ID + "_red.tiff " +\
                      inputdir + ID + "_gre.tiff " + inputdir + ID + "_cellseg.tiff"
            subprocess.call(shlex.split(command))

            # Produce a Displayable RGB Tiff file showing the cell boundary and bud necks:
            command = PMbin + "PMMakeDisplay " + inputdir + ID + ".txt " + inputdir + ID + "_cellseg.tiff " \
                      + inputdir + ID + "_red.tiff " + inputdir + ID + "_gre.tiff " + inputdir + ID + "_preview.tiff "
            subprocess.call(shlex.split(command))

            # Clean up intermediate files
            command = "rm " + inputdir + ID + "_red.tiff"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_gre.tiff"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_dist.tiff"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + ".scp"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_watershed.tiff"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_conf_matrices.scp"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_ellfit.mcv"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_cellseg.tiff"
            subprocess.call(shlex.split(command))
            command = "rm " + inputdir + ID + "_seg.tiff"
            subprocess.call(shlex.split(command))
