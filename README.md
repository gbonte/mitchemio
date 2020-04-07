# ULB MLG results related to the Pseudomonas aeruginosa dataset

This contains the prediction and feature assessment research work carried out in
[MLG ULB](http://mlg.ulb.ac.be) by [G. Bontempi](http://di.ulb.ac.be/map/gbonte/Welcome.html) and M. Mehrian. 

The dataset is released under the Public Domain license, i.e., unrestricted. We
ask that you cite the original, unaltered dataset as 

```
J. Stokes (March 2020), Pseudomonas aeruginosa dataset.
https://www.aicures.mit.edu/data
```

## About the data

The files `train.csv` and `test.csv` contain the SMILES of molecules tested for
activity against Pseudomonas aeruginosa.

To compare your results against our numbers, we have included the 10 splits used
for cross-validation, where indices index into the train set. 


## Feature extraction 
It is done in the script **mread.R**

3 types of feature extraction are performed thanks to the R rcdk package

* Molecular description of smiles: saved in file `descr.Rdata`
* Distance of smiles from known smiles (contained in the SMILES_COVID19-Sheet1.csv file)
* Distance of smiles from positive smiles in the training set

## Feature assessment 
It is done in the script **predfs.R**

## Prediction  testset
It is done in the script predts.R

