This inference is carried out using only T3A and T3C1 cases.
Use the augmentation function parameters in iteration 81 to obtain optimal results.

To make executables:
```
cd JOEfiles
make -j
cd ..
ln -s JOEfiles/joe_direct .
ln -s JOEfiles/joe_adjoint .
```

To recreate inference results, create the following folders in the appropriate directory (Baseline
or one of the iteration directories) and run appropriate simulations using the restart files provided:
```
DIRECT/T3A/
DIRECT/T3C1/
ADJOINT/T3A/
ADJOINT/T3C1/
```
