

This is the R package "Rstringmol", a tool for automating the analysis of runs of the stringmol automata chemistry. 



# Installing Rstringmol

The Rstringmol package was built using devtools, and is hosted on github. To install the package in R, you'll need to have devtools installed first:

```
install.packages("devtools")
```

Then install Rstringmol itself:

```
devtools::install_github("Rstringmol",build_vignettes = T,force = T,subdir = "Rstringmol")
```

# Linux-Ubuntu notes

If the devtools install fails, you may need to install curl, which has to be done outside of R in a Terminal. Type the following: 

```
sudo apt-get install libcurl4-openssl-dev
```
