#define SERIAL_BEG if(mpi_rank>0) MPI_Recv(NULL, 0, MPI_INT, mpi_rank-1, 123, mpi_comm, MPI_STATUS_IGNORE);
#define SERIAL_END if(mpi_rank<mpi_size-1) MPI_Send(NULL, 0, MPI_INT, mpi_rank+1, 123, mpi_comm); MPI_Barrier(mpi_comm);
#define PRINT(...) MPI_Barrier(mpi_comm); if(mpi_rank==0) printf(__VA_ARGS__); MPI_Barrier(mpi_comm);
#define PRINT_SERIAL(...) if(mpi_rank==0) printf(__VA_ARGS__);
#define FOPEN_WRITE_SERIAL(fp, name) if(mpi_rank==0) fp = fopen(name,"w"); else fp = fopen(name,"a");
