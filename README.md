# cancer-hic-norm
Normalization of cancer Hi-C data

The project contains the scripts used in the manuscript : "Effective normalization for copy number variation in Hi-C data. Servant N, Varoquaux N, Heard E, Vert JP, Barillot E."

All results and data used in the publication are available at http://members.cbio.mines-paristech.fr/~nvaroquaux/normalization/

The sources are organized as follow :  
- simulation_model

Contains the script to simulation cancer Hi-C data from diploid real Hi-C data

- cnv_from_hic

Contains the script to infer the CNV profile from any Hi-C dataset

- CNV_norm

Scripts to apply the LOIC and CAIC normalization.  
Note that both methods are now included in the [iced](https://pypi.org/project/iced/#files) python package (>0.5.0).

In case of question, please contact ; nicolas.servant@curie.fr, nelle.varoquaux@gmail.com
