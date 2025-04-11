use the finite difference method to sove the 3d pnp-ns equation for proteins in nanopore   
i made following attempt
1. use the oct tree to mesh protein adaptively.
2. build the matrix
3. solve it with amgx, gmres + amg
## trouble
1. the trouble is the fine mesh of the electrical double layer (million equation), i could not hadnle it on my own pc
2. the ns equation needs schur complement approach
3. i do not know why the pnp equation can not converge nicely (during sovling the jacobian matrix), i guess the because of the EDL, i could not make them finer. 

i give up this, i convert to fenicsx, much much easier, but still have trouble.  


## some results 

### mesh 
![image](https://github.com/user-attachments/assets/ee81ae21-7973-4381-b12f-fac949abe541)
### fluids  (from fenics)
![image](https://github.com/user-attachments/assets/b60db630-e927-4900-935c-7898752810d2)
### potential   
![image](https://github.com/user-attachments/assets/22a5792a-5f98-4f57-bb88-3dd4e1d8b386)
