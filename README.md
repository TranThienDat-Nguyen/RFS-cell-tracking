# Labeled RFS cell tracking algorithms
These codes are implementation of cell tracking and lineage inference algorithms based on labeled random finite set. The algorithms are described in the following paper:\
@article{NVVKC2021CellTracking, \
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; title={Tracking Cells and their Lineages via Labeled Random Finite Sets},\
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; author={Tran Thien Dat Nguyen and Ba-Ngu Vo and Ba-Tuong Vo and Du Yong Kim and Yu Suk Choi},\
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; year={2021},\
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; journal={IEEE Transactions on Signal Processing}\
}. \
Preprint version of this paper is avaialble at https://arxiv.org/pdf/2104.10964.pdf. \
Default settings are obtained from the paper, please follow the comments in demo files to adjust the parameters for your specific applications. The codes are developed for research purposes and not optimized for speed.
# Requirements
This implementation requires MATLAB with following toolboxes: Image Processing Toolbox, Statistics and Machine Learning Toolbox, Parallel Computing Toolbox, and Fuzzy Logic Toolbox. \
The codes were tested on Windows and Linux machines.
# Usages
cd to the folder contains the files 'demo_exp1.m', 'demo_exp2.m', and 'demo_exp3.m'. \
Run demo_exp1.m for experiment 1 in the paper.\
Run demo_exp2.m for experiment 2 in the paper.\
Run demo_exp3.m for experiment 3 in the paper.\
The filters are designed to run on multiple CPU cores, if you do not want this feature, please consider changing the 'parfor' loop to 'for' loop where appropriate. If you do not have enough memory, please consider reducing the number of cores in your parallel computation setting.
# Testing datasets
The synthetic datasets (in experiment 2) are generated using simcep simulation tool with certain modifications to account for cell division. \
The breast cancer cells dataset (MDA-MB-231) used in experiement 3 is provided in folder 'data/img_exp3'.
# Acknowledgments
This implementation is based on MATLAB RFS tracking toolbox provided by Prof. Ba-Tuong Vo at http://ba-tuong.vo-au.com/codes.html. \
The computation of TRA score is based on BaxterAlgorithms cell tracking package provided by Dr. Klas Magnusson at https://github.com/klasma/BaxterAlgorithms. \
The simcep simulation is provided by Dr. Antti Lehmussola based on the paper https://ieeexplore.ieee.org/document/4265752. \
The MPHD spot detector algorithm is provided by Dr. Hamid Rezatofighi at http://users.cecs.anu.edu.au/~hrezatofighi/publications.htm. \
For more details on licenses and information of softwares, functions used in this implementation, see 'third_party' folder.
# Contact
For any queries please contact me at tranthiendat.nguyen@gmail.com.\
Copyright (C) 2021, Tran Thien Dat Nguyen.
