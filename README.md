# ULB MLG results related to the [Pseudomonas aeruginosa dataset](https://www.aicures.mit.edu/data)

This contains the prediction and feature assessment research work carried out in
[MLG ULB](http://mlg.ulb.ac.be) by [G. Bontempi](http://di.ulb.ac.be/map/gbonte/Welcome.html) and 
[M. Mehrian](https://mlg.ulb.ac.be/wordpress/members-2/mmehrian/). 

The dataset is released under the Public Domain license, i.e., unrestricted. We
ask that you cite the original, unaltered dataset as 

```
J. Stokes (March 2020), Pseudomonas aeruginosa dataset.
https://www.aicures.mit.edu/data
```

## About the data

The files `train.csv` and `test.csv` contain the SMILES of molecules tested for
activity against Pseudomonas aeruginosa and made available by [MIT AI cures](https://www.aicures.mit.edu/data).


## Known smiles
They are contained in the file `SMILES_COVID19-Sheet1.csv` curated by M. Mehrian.
The related Google sheets is [here](https://docs.google.com/spreadsheets/d/1Ll26liuImbjxnkfwunEBb9Hn9nH38lUvOsLJaYLayFQ/edit#gid=0).
It contains 130 molecules related to COVID and 70 negative controls.

If you are interested in update it please let us know.

## Feature extraction 
It is done in the script `mread.R`

3 types of feature extraction are performed thanks to the R rcdk package

* Molecular description of smiles: saved in file `descr.Rdata`
* Distance of smiles from known smiles (see above): saved in file `dmetrics1.Rdata`
* Distance of smiles from positive smiles in the training set: saved in file `dmetrics2.Rdata`

## Feature assessment 
The script `predfs.R`  selects the most informative features (among the ones in `dmetrics1.Rdata` )
and assesses the significativity of the choice of known COVID molecules.

The current version of feature selection script (combination of MRMR and Random Forest importance)
returns the following 16 molecules as the most relevant ones:

**Moxifloxacin (2 forms),       Anakinra ,          dexlansoprazole,     mesalamine ,   Ceftriaxone ,       Clomiphene citrate, Ivacaftor,          Gilteritinib,  Bazedoxifene,       Ribavirin,          Carvedilol, Niclosamide,  Oxyclozanide  ,     Droloxifene   ,     Lercanidipine**

Only 4 are negative controls (hypergeometric p-val 0.0008249591) were selected.
Associated cross-validated AUC is 0.853.

The preliminary conclusion is then that COVID related molecules are informative about the label of this dataset. Comments welcome.

## Boxplot distances
Script `visual.R`: it shows the boxplot of the distances of the training set molecules (distinguishing between class 0 and class 1) with respect to each of the molecules selected above.  This is a visualization of the discrimination power of the distance features.

## Prediction  testset
It is done in the script `predts.R`

## Comments/remarks
**Please consider that those are preliminary results made available for discussion.**
**We would be happy to receive comments and remarks.**
Please contact [Pr. G. Bontempi](mailto:gbonte@ulb.ac.be?subject=[mitchemio:Github]).

