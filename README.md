# The nHODLRdD and s-nHODLRdD accelerated iterative solver
This repository mainly contains two new fast iterative solvers (GMRES), accelerated using the following fast matrix-vector product algorithms.
1. The nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($nHODLRdD$)
2. The semi-nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($s-nHODLRdD$)
   
The $nHODLRdD$ and $s-nHODLRdD$ algorithms are the nested and sem-nested versions of the previously proposed $HODLRdD$ algorithm. The $HODLRdD$ fast algorithm can be found [here](https://github.com/SAFRAN-LAB/HODLRdD), which works for any user-given dimension $d$. Due to use of the nested / semi-nested bases the $nHODLRdD$ and $s-nHODLRdD$ algorithms are faster than the $HODLRdD$ algorithm.
Users can also run the $HODLRdD$ accelerated iterative solver using this repository by changing the flag in `Makefile` and the codes are self-explanatory. Currently, this repository works for $d=1,2,3$, i.e., in $1D$, $2D$ and $3D$.
______________________________________________________________

# How to run the code
User can run $3$ different types of fast iterative solver using this repository.
1. $nHODLRdD$ accelerated iterative solver.
2. $s-nHODLRdD$ accelerated iterative solver.
3. $HODLRdD$ accelerated iterative solver.
_____________________

### Dependencies
1. The [**Eigen**](https://eigen.tuxfamily.org) linear algebra libray.
2. The [**Boost**](https://www.boost.org/) library.
3. An OpenMP enabled modern compiler (optional).
______________________________________________________________________

### Download the Source code
The source code available in this repository can be cloned using:
```
git clone https://github.com/riteshkhan/nHODLRdD.git --recursive
```
______________________________________________________________________

### Testing
Make sure that you have set the path of the libraries correctly. The path can be modified in the `Makefile`. For testing purpose we have added two options (i) Integral Equation solver and (ii) RBF interpolation.

To test the code simply run the following commands in the terminal.
```
user@computer nHODLRdD$ mkdir build && cd build
```

Now the user can go to any directory to run the corresponding codes in $1D$ or $2D$ or $3D$. For example, to run the code in $1D$

```
user@computer 1D/codes$ make -f Makefile1D.mk clean && make -f Makefile1D.mk
```
In $1D$ if we set the following inputs in the `main` file, 
1. `atoi(argv[1])` $= N = 250000$ (Number of particles)
2. `atoi(argv[2])` $=n_{max} = 100$ (Number of maximum particles in leaf clusters)
3. `atoi(argv[3])` $=L=1$ (Semi-length of the cluster)
4. `atoi(argv[4])` $=12$ (NCA/ACA tolerance $= 10^{-12}$ and $\epsilon_{GMRES} = 10^{-12}$)
5. `atoi(argv[5])` $=1$ (Intergral equation solver=0 or RBF interpolation=1)

```
user@computer 1D/codes$ ./test 250000 100 1 12 1
```

If everything is fine one might get output like below.
```txt
GMRES Parameters 
Maximum Iterations : 500
GMRES Tolerance : 1e-12

Reached Solution before Max_Iterations 
Reached the desired tol after 11 iterations
Final resid 6.81945e-14
********** Summary of HSS / HBS / nHODLR1D accelerated GMRES to solve a system **********


The number of particles taken: 250000 and choice =  RBF interpolation

The maximum number of particles at leaf clusters: 100

Depth of the tree: 12

The final residual error is: 6.81945e-14

The total number of iteration is: 11

Total Assembly time: 47.0822s

GMRES time: 1.49966s

Storage (in GB): 0.658597 GB

Compression ratio: 0.00131719

The (norm-2) relative error in solution: 1.99987e-12
=========================================================================================
```

```txt
GMRES Parameters 
Maximum Iterations : 500
GMRES Tolerance : 1e-12

Reached Solution before Max_Iterations 
Reached the desired tol after 11 iterations
Final resid 6.63218e-14
********** Summary of HODLR / HODLR1D accelerated GMRES to solve a system **********


The number of particles taken: 250000 and choice =  RBF interpolation

The maximum number of particles at leaf clusters: 100

Depth of the tree: 12

The final residual error is: 6.63218e-14

The total number of iteration is: 11

Total Assembly time: 26.4311s

GMRES time: 2.37944s

Storage (in GB): 1.61 GB

Compression ratio: 0.00322

The (norm-2) relative error in solution: 8.68101e-13
====================================================================================


```
# References
[Paper](https://arxiv.org/pdf/2309.14085.pdf)
