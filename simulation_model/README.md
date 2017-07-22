## Simulation of cancer Hi-C data

These scripts were developed to simulate cancer Hi-C data from diploid real data.
They require a recent version of R with the following packages ; HiTC, rtracklayer, Matrix, ggplot2, reshape2


To use them, you need ;  
1- A diploid dataset in the HiC-Pro format, i.e. triplet sparse format with a BED file that describes the bin coordinates.  
2- A CNV file which simply contains the intervals of altered segments with their copy number to simulate (see the files folder for example).  

Then, use the simulateCNV.sh script to run the code.  

>> bash simulateCNV.sh  
usage : simulateCNV.sh -i INPUT -b INPUT_BED -c CNV_FILE -o OUTPUT [-h]  
Use option -h|--help for more information  

In case of question, please contact nicolas.servant@curie.fr
