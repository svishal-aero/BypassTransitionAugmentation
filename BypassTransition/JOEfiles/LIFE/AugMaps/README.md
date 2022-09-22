# Class: GenAugMap (Generic Augmentation Map)

This is a virtual class with purely virtual methods. Only fully defined child classes can be instantiated.

## Protected Variables

- `mpi_rank` is the MPI rank of the process
- `mpi_size` is the number of MPI processes
- `name` is the name given to the augmentation function - files containing parameters and sensitivities are named accordingly.
- `params` is a 1D vector of doubles containing parameters.
- `dJ_dParams` is the sensitivity of the objective w.r.t. the parameter vector.

## Purely Virtual Public Methods

- `init` allocates the vectors `params` and `dJ_dParams`.
- `calc_beta` takes in features (`std::array<double, [num_features]>`), and calculates the augmentation (`double`) and its sensitivity w.r.t. features (`std::array<double, [num_features]>`).
- `update_dJ_dParams` takes in the sensitivity of the objective w.r.t. the augmentation (`double`) along with the features (`std::array<double, [num_features]>`), and updates `dJ_dParams`.
