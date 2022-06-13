# Labeled RFS cell tracking algorithms
These codes are implementation of cell tracking and lineage inference algorithms based on labeled random finite set. The algorithms are described in the following paper:
```
@ARTICLE{9536381,
  author={Nguyen, Tran Thien Dat and Vo, Ba-Ngu and Vo, Ba-Tuong and Kim, Du Yong and Choi, Yu Suk},
  journal={IEEE Transactions on Signal Processing},
  title={Tracking Cells and Their Lineages Via Labeled Random Finite Sets},
  year={2021},
  volume={69},
  number={},
  pages={5611-5626},
  doi={10.1109/TSP.2021.3111705}}
 ```
The paper is available at https://ieeexplore.ieee.org/document/9536381 (Open Access under CC BY 4.0 License).\
Default settings are obtained from the paper, please follow the comments in demo files to adjust the parameters for your specific applications. The codes are developed for research purposes and not optimized for speed.
# Requirements
This implementation requires MATLAB with following toolboxes: Image Processing Toolbox, Statistics and Machine Learning Toolbox, Parallel Computing Toolbox, and Fuzzy Logic Toolbox. \
The codes were tested on Windows computer with MATLAB 2019b and Linux computer with  MATLAB 2020b.
# Usages
cd to the folder contains the files 'demo_exp1.m', 'demo_exp2.m', and 'demo_exp3.m'. \
Run demo_exp1.m for tracking with simulated data.\
Run demo_exp2.m for tracking with synthetic cell migration datasets.\
Run demo_exp3.m for tracking with real migration sequence of breast cancer cells.\
The filters are designed to run on multiple CPU cores, if you do not want this feature, please consider changing the 'parfor' loop to 'for' loop where appropriate. If you do not have enough memory, please consider reducing the number of cores in your parallel computation setting. \
The filters will generate a folder 'output' while running to store temporary information for analysis. \
Note: the estimated results at each time step (as shown in MATLAB console) are only indicative since the recursive estimator corrects the past estimates given observed data upto current time step.
# Datasets
## Synthetic datasets
The synthetic datasets (in experiment 2) are generated using simcep simulation tool with certain modifications to account for cell division. When the file 'demo_exp2.m' is run, the codes will randomly generate synthetic sequences (using the prescribed ground truth and noise parameters but with random realization of noise). The generated sequences will be stored in 'output/Synthetic_X' folder. If you just want the datasets without running the codes, they are avaialble in the 'datasets' folder.  'synthetic_X' folders contain the synthetic cell migration sequences from the best quality 'synthetic_01' to the worst quality 'synthetic_05'. Folder 'synthetic_truth' contains the ground truth annotated using Cell Tracking Challenge format (https://public.celltrackingchallenge.net/documents/Naming%20and%20file%20content%20conventions.pdf).
## MDA-MB-231 dataset
The breast cancer cells dataset (MDA-MB-231) used in experiement 3 is provided in folder 'datasets/MDA_MB_231'. \
Datasets and labels (manual annotations) are also included for training purpose. 
# Acknowledgments
This implementation is based on MATLAB RFS tracking toolbox provided by Prof. Ba-Tuong Vo at http://ba-tuong.vo-au.com/codes.html. \
The computation of TRA score is based on BaxterAlgorithms cell tracking package provided by Dr. Klas Magnusson at https://github.com/klasma/BaxterAlgorithms. \
The simcep simulation is provided by Dr. Antti Lehmussola based on the paper https://ieeexplore.ieee.org/document/4265752. \
The MPHD spot detector algorithm is provided by Dr. Hamid Rezatofighi at http://users.cecs.anu.edu.au/~hrezatofighi/publications.htm. \
For more details on licenses and information of softwares, functions used in this implementation, see 'third_party' folder.
# Contact
For any queries please contact me at tranthiendat.nguyen@gmail.com.\
Copyright (C) 2021, Tran Thien Dat Nguyen.
