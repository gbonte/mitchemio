# ULB MLG results related to the [Pseudomonas aeruginosa dataset](https://www.aicures.mit.edu/data)

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

## Known smiles
They are contained in the file `SMILES_COVID19-Sheet1.csv` curated by M. Mehrian.

## Feature extraction 
It is done in the script `mread.R`

3 types of feature extraction are performed thanks to the R rcdk package

* Molecular description of smiles: saved in file `descr.Rdata`
* Distance of smiles from known smiles: saved in file `dmetrics1.Rdata`
* Distance of smiles from positive smiles in the training set: saved in file `dmetrics2.Rdata`

## Feature assessment 
It is done in the script `predfs.R`. in particular this script assess the significativity of the choice of known COVID molecules

## Prediction  testset
It is done in the script `predts.R`

