# SymmetrySTEM

Application to analyse SymmetrySTEM signals in 4D-STEM datasets. If you find the application useful in your reaserch please cite the following paper:

> Matus Krajnak, Joanne Etheridge  
> A symmetry-derived mechanism for atomic resolution imaging  
> Proceedings of the National Academy of Sciences Nov 2020, 117 (45) 27805-27810  
> https://doi.org/10.1073/pnas.2006975117


To build the application, you will require a Linux PC, standard build system libraries and ArrayFire library.
The executable can be compiled by:

``` 
mkdir build && cd build
cmake ..
make
```

The code then can be used as by executing ```./symmetry``` in terminal. Running the code without any parameters will give a help.
