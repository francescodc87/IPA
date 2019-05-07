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
higher than 5% were included in the dataset. A realistic experimental outcome was simulated by adding Gaussian noise to the m/z values, where appropriate the standard deviation was computed for each mass assuming an instrument accuracy of 10 ppm. For each detected mass, the measured intensities were chosen in order to be coherent with the theoretical abundance of the isotopes, and the same RT value (+-2 seconds) was assigned to all the simulated masses related to the same compound. After loading the package, the synthetic experiment just described can be imported in the R environment using the following commands:
```
library(IPA)
install_github("francescodc87/IPA")
```




## Built With

* [R](https://www.r-project.org/)
* [RStudio](https://www.rstudio.com/)


## Authors

* **Francesco Del Carratore**
* **Simon Rogers**
* **Rainer Breitling**
