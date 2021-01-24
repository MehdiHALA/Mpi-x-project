### Mpi-x-project
TP de Programmation hybride MPI + X du module IT389

## Load Compilers

```sh
module load compiler/gcc/9.1.0 compiler/intel/2019_update4 mpi/openmpi/4.0.3
```

## Compilation 

```sh
make
```
PS : J'ai encore un problème que je n'arrive pas à résoudre. Quand je push le Makefile, il détecte les tabs de celui-ci comme étant 8 espaces ainsi le Makefile ne fonctionne pas donc il faut refaire toutes les tabulations

## Execution
Pour la version basique et la version omp:
```sh
./build/stencil-%
```
Pour la version MPI et MPI+X:
```sh
mpirun -np 4 ./build/stencil%
```

## Script de génération 
```sh
./scripts/generation.sh
```
Il faudra changer le nom de l’éxecutable à l'intérieur du fichier bash si on souhaite tester un autre que la version basique.
