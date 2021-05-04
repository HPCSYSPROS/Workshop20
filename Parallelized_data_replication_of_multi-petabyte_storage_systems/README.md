# Parallelized Data Replication of Multi-Petabyte Storage Systems	

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737878.svg)](https://doi.org/10.5281/zenodo.4737878)
#[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4327599.svg)](https://doi.org/10.5281/zenodo.4327599)

**Authors**
* Honwai Leong, DDN
* Andrew Janke, The University of Sydney
* Daniel Richards, DDN
* Stephen Kolmann, The University of Sydney

**Abstract:**
This paper presents the architecture of a highly parallelized data replication workflow implemented at The University of Sydney that forms the disaster recovery strategy for two 8-petabyte research data storage systems at the University. The solution leverages DDN’s GRIDScaler appliances, the information lifecycle management feature of the IBM Spectrum Scale File System, rsync, GNU Parallel and the MPI dsync tool from mpiFileUtils. It achieves high performance asynchronous data replication between two storage systems at sites 40km’s apart. In this paper, the methodology, performance benchmarks, technical challenges encountered and fine-tuning improvements in the implementation are presented.
