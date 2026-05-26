# AST5220-Cosmology-II

## About this repository
This repository contains the numerical work done in the course AST5220 at the University of Oslo Spring 2026.
The code follows the structure provided by [Hans Winther](https://github.com/HAWinther/AST5220-Cosmology), which again follows the general structure of [Petter Callin](https://arxiv.org/pdf/astro-ph/0606683) and the lecture notes developed by [Hans Winther](https://cmb.wintherscoming.no).


The repository contains a `C++` Einstein-Boltzmann solver (a CAMB like code), which culminates in reproducing the CMB power spectrum for a simplified model. The model developed in this repository is a simplification intended for master students leaving out helium, reionization, polarization and neutrinos. The plotting is implemented in `python`. The figures produced by the `m*_plots.py` can be found under `figures`. Simply comment out the desired plots at the bottom of the plotting functions. For the plotting you need the following packages:

```
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
```


## Compiling
You can compile the code by running 
```
make
```
in the terminal, but be aware that you need to have downloaded the `GSL` library. See below for a guide (copied from [Winther's GitHub](https://github.com/HAWinther/AST5220-Cosmology)). To run the simulation, simply run 
```
./cmb
```
and the simulation should start. Note that you need to make sure that relevant lines in the `Main.cpp` file are not commented out. For more information, inspect `Main.cpp` and comments therein.

## GSL library
Before you can run this code, you must ensure that you have installed the `GSL` library. See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory:
```
cd $HOME
```
- Make a local folder to hold libraries:
```
mkdir local
```
- Enter this directory:
```
cd local
```
- Download the code (if you don't have wget you need to get the file to this dir by other means):
```
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
```
- Untar the code:
```
tar -xvf gsl-2.6.tar.gz
```
- You should now have the gsl-2.6 folder. Enter it:
```
cd gsl-2.6
```
- Run the configure script:
```
./configure --prefix=$HOME/local
```
- Compile and install it:
```
make ; make install
```
- In the CMB code Makefile change the include and lib paths to point to the library:
```
INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas
```
- If this fails with "libgsl.so not found" then run the command:
```
export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"
```
and try to run `./cmb` again and it should work. To avoid having
to run this command every time you open a new terminal open
the `$HOME/.bashrc` file and add this line to the end of the file
and it will load everytime you open a new window.


