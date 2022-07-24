

This is the R package "Rstringmol", a tool for automating the analysis of runs of the stringmol automata chemistry. 



## Installing Rstringmol

The Rstringmol package was built using devtools, and is hosted on github. To install the package in R, you'll need to have devtools installed first:

```
install.packages("devtools")
```

Then install Rstringmol itself:

```
devtools::install_github("Rstringmol",build_vignettes = T,force = T,subdir = "Rstringmol")
```

### Vignettes

Long-form help documentation is available via R vignettes for this package. Each
vignette describes how to perform the analysis described in the academic papers listed below.

To access the vignettes from within R, enter: 

```
browseVignettes(package="Rstringmol")
```

..this should open a browser window with links to the following vignettes:

- **stringmol-analysis** which describes the analysis for paper [1] 
- **stringmol-analysis-II** which describes the analysis for paper [2]



### Linux-Ubuntu notes

If the devtools install fails in linux, you may need to install curl, which has to be done outside of R in a Terminal. Type the following: 

```
sudo apt-get install libcurl4-openssl-dev
```

## Python code and documentation

This package now has some python code which was used to generate the graphics in 
paper [2]. Documentation for this can be found in the `docs/` sub-directory. 


## Papers using Rstringmol


- [1] Hickinbotham, S., Stepney, S. and Hogeweg, P., "Nothing in evolution makes sense except in the light of parasitism: evolution of complex replication strategies", *Royal Society Open Science* **8 (8)** 210441
- [2] Stepney, S. and Hickinbotham, S. "On the Open-Endendness of Detecting Open-Endedness", *Artificial Life Journal, In press*


