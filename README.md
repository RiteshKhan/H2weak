# The $nHODLRdD$ and $s-nHODLRdD$ accelerated iterative solver
This repository mainly contains two new fast iterative solvers (GMRES), accelerated using the following fast matrix-vector product algorithms.
1. The nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($nHODLRdD$)
2. The semi-nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($s-nHODLRdD$)
   
The $nHODLRdD$ and $s-nHODLRdD$ algorithms are the nested and sem-nested versions of the previously proposed $HODLRdD$ algorithm. The previously developed $HODLRdD$ fast algorithm can be found [here](https://github.com/SAFRAN-LAB/HODLRdD), which works for any user-given dimension $d$. Due to use of the nested / semi-nested bases the $nHODLRdD$ and $s-nHODLRdD$ algorithms are faster than the $HODLRdD$ algorithm.
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

### Downloading the Source code
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
user@computer 1D$ make -f Makefile1D.mk clean && make -f Makefile1D.mk
```
By default it will work for $N = 100000$, $n_{max} = 100$ and $NCA/ACA tol. = 10^{-12}$
