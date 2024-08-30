This is the code for self-consistent mean field solution of the Kitaev+Zeeman+Kappa Hamiltonian. In this code along with checking iterative self-onsistency I use lagrange multipliers for satisfying the constraint for the solution to be in the physical subspace. I use the hybrd1 solver of minpack.

The hybrd solver solves for a set of n equations in n variables. This set of equations is given in the subroutine fcn. In this code the set of equations is the six constraint equations 3 for each sublattic. hybrd then solves for lagrange mulipliers that ensure that the constraints are satisfied.

The Hamiltonian I solve for is

$H=K\sum_{\langle ij\rangle_{\alpha}}\sigma_i^a\sigma_j^a+\sum+{i}(h^x,h^y,h^z).(\sigma_i^x,\sigma_i^y,\sigma_i^z)+\kappa\sum_{<jlk>_{R_2,R_I}}\sigma_j^a\sigma_k^b\sigma_l^c$


The sign of \Kappa depends on the sign we take for the zeeman term. If H_zeeman=h.S, \kappa is positive. If H_zeeman=-h.s, Kappa is negative. In this code I have chosen this sign to be positive by default. We can switch the zeeman/kappa term on or off by choosing zeeman_info/kappa_info = 1 or 0 in the input file para.in

Note that we can also switch the signs of the zeeman and Kappa term to be negative by choosing zeeman_info, kappa_info to be -1.

This is because in the code I construct my hamiltonian as

$H=KH_K+zeeman_infoH_zeeman+kappa_infoH_{\kappa}$
