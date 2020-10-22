# ptatin-timeseries-svd-data
Collection of PETSc vecs to be stacked into a snapshot matrix.

The collection is a time series of temperature fields generated from a 3D time dependent model of continental rifting. The forward model used was pTatin3D. The time steps saved in the diretory `data/` include:

* 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400

File names reflect the timestep each file is associated with.



See file `data/ptatin.options-2020.09.18_23:27:48` for rift model specific information.



### Dependencies

* PETSc (I used v3.13 but 3.14 will also be fine )



### Compiling 

* Make sure the environment variables `PETSC_DIR` and `PETSC_ARCH` are defined.

* Then execute the following command

* ```
  make -f $PETSC_DIR/share/petsc/Makefile.user v2m
  ```

  



### Notes

Executing `./v2m` will collect a set of snapshots (temperature solutions), assemble them into a single PETSc Mat object and write the matrix out into a binary file called `snapshot.pbmat`.

The size of the default snapshot matrix is 1094049 x 13.

The number of snapshots can be reduced via the command line argument 

```
-truncate <value>
```

where `<value>` is an integer and must be less than or equal to 13.

A vts file (loadable by ParaView) called `snapshot0.vts` is also created. It contains the temperature field associated with the first vector in snapshot series. The coordinates associated with the mesh plotted are not the true coordinates, but rather a hard coded set of values for the domain [0,12] x [-1.5,0] x [0,6].

The solution and true coordinates of the domain for a given time step can be visualized with the option

```
-solution_view <step>
```

where `<step>` is an integer. The VTS file generated will be called `step<step>_temperature.vts`.

The function `PetscErrorCode SnapshotMatCreate(Mat *snapshots)` can be re-used in a SLEPc application to directly create and load the snapshot matrix. This might be desired as the binary ouput of tha matrix is signficantly larger than simply the sum of the individual `*.pbvec` files.



### Reference (reduced) SVD result

The singular values for the full data set are

```
singular values:
 [2.22402501e+05 7.79461315e+04 3.24114516e+04 1.74890097e+04
 9.56670518e+03 5.36755762e+03 2.99723043e+03 2.06755101e+03
 1.31076808e+03 6.93460366e+02 4.35445057e+02 2.90975623e+02
 4.17900032e-10]
```

These can be computed using `python numpy-svd.py`. On Dave's laptop, computing the reduced SVD using `numpy.linalg.svd()` required < 2 seconds.

```
np.linalg.svd: 1.3746e+00 (sec)
```

Note that the snapshot matrix was constructed such that it had zero mean. See lines 45-47 of `numpy-svd.py`. If lines 45-47 are commented out, the singular values obtained are

```
singular values:
 [3.90129922e+06 2.15015904e+05 7.78146522e+04 3.16355182e+04
 1.72906096e+04 9.11472586e+03 5.25599082e+03 2.89014839e+03
 2.05205149e+03 1.30710758e+03 6.93455391e+02 4.34965945e+02
 2.90971420e+02]
```



