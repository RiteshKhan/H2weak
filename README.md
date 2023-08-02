# The $nHODLRdD$ and $s-nHODLRdD$ accelerated iterative solver 
This repository mainly contains two different fast iterative solvers (GMRES), accelerated using the following fast matrix-vector product algorithms.
1. The nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($nHODLRdD$)
2. The semi-nested hierarchically off-diagonal low-rank matrix in $d$ dimensions ($s-nHODLRdD$)
   
The $nHODLRdD$ and $s-nHODLRdD$ algorithms are the nested and sem-nested versions of the previously proposed $HODLRdD$ algorithm. The previously developed $HODLRdD$ fast algorithm can be found [here](https://github.com/SAFRAN-LAB/HODLRdD), which works for any user-given dimension $d$. Due to use of the nested / semi-nested bases the $nHODLRdD$ and $s-nHODLRdD$ algorithms are faster than the $HODLRdD$ algorithm.
Users can also run the $HODLRdD$ accelerated iterative solver using this repository by changing the flag in `Makefile`. Currently, this repository works for $d=1,2,3$, i.e., in $1D$, $2D$ and $3D$.

# How to run the code
User can run $3$ different types of fast iterative solver using this repository.
1. $nHODLRdD$ accelerated iterative solver.
2. $s-nHODLRdD$ accelerated iterative solver.
3. $HODLRdD$ accelerated iterative solver.


<u>sometext</u>
