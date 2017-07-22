## Infer CNV profile from Hi-C data

These scripts were developed to segment 1D Hi-C data and thus estimate the CNV status.
They require a recent version of R with the following packages ; HiTC, rtracklayer, Matrix, RColorBrewer, ggplot2, gridExtra

To use them, you need ;  
1- A Hi-C dataset in the HiC-Pro format, i.e. triplet sparse format with a BED file that describes the bin coordinates.

Then, use the segment_hic_data.sh script to run the code.

>> bash segment_hic_data.sh  
usage : segment_hic_data.sh -i INPUT -b INPUT_BED -s BKP_STRINGENCY -o ODIR [-h]  
BKP_STRINGENCY must be normal, low, high  

The script will annotate the Hi-C data using GC content, fragment length and mappability information. According to the data resolution, this part can take some times.  
It then run the segmentation procedure on the annotated data.  

In case of question, please contact nicolas.servant@curie.fr
