use the finite difference method to sove the 3d pnp-ns equation for proteins in nanopore   
i made following attempt
1. use the oct tree to mesh protein adaptively.
2. build the matrix
3. solve it with amgx, gmres + amg
## trouble
1. the trouble is the fine mesh of the electrical double layer (million equation), i could not hadnle it on my own pc
2. the ns equation needs schur complement approach
3. i do not know why the pnp equation can not converge nicely (during sovling the jacobian matrix), i guess the because of the EDL, i could not make them finer. 

i give up this, i convert to fenicsx, much much easier. 
