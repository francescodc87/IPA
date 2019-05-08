# IPA
Integrated Probabilistic Annotation (IPA) - A Bayesian annotation method for LC/MS data integrating biochemical relations, isotope patterns and adduct formation

## Getting started
Not sure if I need this section but I am keeping it for now

### Prerequisites
Prior installation, it is necessary to download and install
[R](https://www.r-project.org/) version 3.4.3 or later. R is a programming language and
free software environment. If you are unfamiliar with the R syntax and commands, it is
recommended to read the [introduction to R](https://cran.r-project.org/doc/manuals/R-intro.html)
document first.
<br />
<br />
The download an installation of [RStudio](https://www.rstudio.com/) is also necessary. RStudio
is a free and open-source integrated development environment for R.
It is also necessary to install the devtools package by simply using the following command in the R session:
```
install.packages(devtools)
```
<br />


## Installation
The installation of the IPA R package is extremely easy and it can be achieved by
simply using the following commands:
```
library(devtools)
install_github("francescodc87/IPA")
```



## Case study - Synthetic experiment

In order to show the functionalities of the IPA package we generated a syntheitc experiment.
This synthetic experiment contains 15 compounds involved in the mevalonate pathway and limonene synthesis.
Considering a positive ionisation, several adducts were simulated for each of the considered metabolites (M+H, M+Na,
M+2H and 2M+H for positive mode). Additionally, all the possible isotopes with a predicted relative abundance
higher than 5% were included in the dataset. A realistic experimental outcome was simulated by adding Gaussian noise to the m/z values, where appropriate the standard deviation was computed for each mass assuming an instrument accuracy of 10 ppm. For each detected mass, the measured intensities were chosen in order to be coherent with the theoretical abundance of the isotopes, and the same RT value (+-2 seconds) was assigned to all the simulated masses related to the same compound. After loading the package, the synthetic experiment just described can be imported in the R environment using the following commands together with the additonal datasets needed:
```
library(IPA)
data("SyntheticExperiment")
data("isotopes")
data("adducts")

```
The synthetic experiment was simulated considering positive ionisation and an instrument accuracy of 10 ppm:

```
ionisation="positive"
ppm=10

```

The first step consist in finding all the possible hits in the database. This is obtained with find.hits() function:


```
Hits <- find.hits(adducts.matrix= all_adducts_POS,
                  dataset=df.POS, ppm.thr= 5*ppm,
                  RTwin=60,relation.id = relation.id,
                  isotopes=isotopes, 
                  iso.threshold=1)

```

The output of this function can be used as input of the compute.Priors() function
in order to evaluate the prior probabilities:

```
Prior <- compute.Priors(Hits=Hits, dataset=df.POS,
                        pk=rep(1,nrow(Hits$all.formulas)),
                        ppm=ppm,unknown.ppm = 3*ppm, pr.lim = 1e-15)

```
In order to consider all the possible connections, it is necessary to build the connectivity matrices for adducts, isotopes and biochemical connections:

```
### building ADD matrix
ADD <- build.add.connenctivity.matrix(Prior=Prior,  DB=DB,
                                      ionisation,fully.connected=FALSE)

### building ISO matrix
ISO <- build.iso.connectivity.matrix(Prior=Prior, DB=DB, ratios=TRUE)

### building BIO matrix
BIO1 <- build.bio.connectivity.matrix(Prior=Prior, DB=DB,
                                      ionisation,
                                      connection.type="reactions")

```
Finally, it is possible to compute the posterior probabilities using the IPAposteriors() function.

```
PostISOADDBIO <- IPAposteriors(P=Prior,Iso = ISO, Add = ADD, Bio = BIO1,
                               Int = as.numeric(df.POS[,3]), ratio.toll = 0.8,
                               delta.iso = .1, delta.add = .1, delta.bio = 1,
                               allsamp = T, no.its = 5000, burn = 1000,
                               rel.id = relation.id)

```

It might be necessary to organize the results in a more readable fashion.
This can be obtained with the ParseIPAresults() function:

```
Final.res <- ParseIPAresults(Post, Prior,dataset =df.POS, DB = DB, IDs=IDs)
```

This creates a list with the same length of the number of feature considered. Each element contains a summary of the results for the considered features.

## Built With

* [R](https://www.r-project.org/)
* [RStudio](https://www.rstudio.com/)


## Authors

* **Francesco Del Carratore**
* **Simon Rogers**
* **Rainer Breitling**
