# ptatin-timeseries-svd-data
Collection of PETSc vecs to be stacked into a snapshot matrix



### Dependencies

* PETSc (I used v3.13 but 3.14 will also be fine )



### Notes

Executing `v2m.app` will collect a set of snapshots (temperature solutions), assemble them into a single PETSc Mat object and write the matrix out into a binary file called `snapshot.pbmat`.

The size of the default snapshot matrix is 1094049 x 13.

The number of snapshots can be reduced via the command line argument 

```
-truncate <value>
```

where `<value>` must be less than or equatl to 13.

A vts file (loadable by ParaView) called `snapshot0.vts` is also created. It contains the temperature field associated with the first vector in snapshot series.

The function `PetscErrorCode SnapshotMatCreate(Mat *snapshots)` can be re-used in a SLEPc application to directly create and load the snapshot matrix. This might be desired as the binary ouput of tha matrix is signficantly larger than simply the sum of the individual `*.pbvec` files.