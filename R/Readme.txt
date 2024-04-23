This folder is R code to achieve the results for the project Correspondence analysis: handling cell-wise outliers via the reconstitution algorithm

Section 3 How outliers originate:

smallmadeupexample.R

Section 5 Empirical studies/Results:

The CA for the original matrix: origbrands.R, origop.R. 
The CA for the reconstitution of order h: reconbrands.R, reconop.R. 
The CA for the MacroPCA: macropcabrands.R, macropcaop.R.
The CA for the supplementary points method: supplebrands.R, suppleop.R.

In addition, reconca.R is a function to implement reconstitution of order h; the code is written by adjusting imputeCA function in the package missMDA version 1.19.

Furthermore, to get the outlying rows in the second step of MacroPCA, we adjust the function MacroPCA in the R package cellWise. The difference is that we extract the information of outlying row in the second step "1. Projection pursuit"  res$Hstar = Hstar and "4. Reweighting" res$H0 = H0. Details for "1. Projection pursuit" and "4. Reweighting" are in Mia Hubert, Peter J. Rousseeuw & Wannes Van den Bossche (2019) MacroPCA: An All-in-One PCA Method Allowing for Missing Values as Well as Cellwise and Rowwise Outliers, Technometrics, 61:4, 459-473, DOI: 10.1080/00401706.2018.1562989.
