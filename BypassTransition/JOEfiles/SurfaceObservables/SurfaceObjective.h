template<typename DBLE>
class SurfaceObjective
{
  public:
    
    int nTargetLocations, nTargets;
    int nFaces, nFacesTotal;
    std::vector<int> indexForFaces;
    std::vector<std::array<double,3>> x_targets;
    std::vector<std::array<double,3>> x_faces;
    std::vector<std::vector<double>> targets, targets_interp;
    std::vector<std::vector<DBLE>> values;
    std::vector<double> lambdas;
    DBLE local_obj;
    std::vector<double> total_obj;
  
  public:
    
    inline void readFile(std::string filename)
    {
      PRINT("\nStarted reading surface targets from file \"%s\"...\n", filename.c_str());
      SERIAL_BEG
        FILE *fp = fopen(filename.c_str(),"r");
        if(fp==NULL) printf("Target file \"%s\" for objective evaluation"
                            " not found\n", filename.c_str());
        else
        {
          fscanf(fp, "%d", &nTargetLocations);
          fscanf(fp, "%d", &nTargets);
          PRINT_SERIAL("nLocations: %d nTargets: %d\n", nTargetLocations, nTargets);
          x_targets.resize(nTargetLocations);
          targets.resize(nTargetLocations);
          lambdas.resize(nTargets);
          total_obj.resize(nTargets);
          for(int j=0; j<nTargets; j++)
          {
            fscanf(fp, "%le", &(lambdas[j]));
          }
          PRINT_SERIAL("lambdas read...\n");
          for(int i=0; i<nTargetLocations; i++)
          {
            targets[i].resize(nTargets);
            fscanf(fp, "%le", &(x_targets[i][0]));
            fscanf(fp, "%le", &(x_targets[i][1]));
            fscanf(fp, "%le", &(x_targets[i][2]));
            for(int j=0; j<nTargets; j++)
            {
              fscanf(fp, "%le", &(targets[i][j]));
            }
          }
          fclose(fp);
        }
        nFaces = 0;
      SERIAL_END
      PRINT("Finished reading surface targets...\n\n");
    }

    inline void evaluate()
    {
      DBLE local_obj_j;
      double local_obj_val;
      for(int j=0; j<nTargets; j++)
      {
        local_obj_j = 0.0;
        for(int i=0; i<nFaces; i++)
        {
          // __CHANGE__ here if objective function needs to be changed
          local_obj_j = local_obj_j + pow(values[i][j]-targets_interp[i][j],2);
        }
        local_obj = local_obj + lambdas[j]*local_obj_j;
        #ifndef ADJOINT_MODE
        double local_obj_val = local_obj_j;
        #else
        double local_obj_val = local_obj_j.value();
        #endif
        MPI_Allreduce(&local_obj_val, &(total_obj[j]), 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
      }
    }

    inline void write(std::string filename)
    {
      FILE *fp;
      PRINT("Writing marker file\n");
      PRINT("nFacesTotal: %d\n", nFacesTotal);
      SERIAL_BEG
        FOPEN_WRITE_SERIAL(fp,filename.c_str());
        #ifndef ADJOINT_MODE
        for(int i=0; i<nFaces; i++)
        {
          fprintf(fp, "%+.10le %+.10le %+.10le", x_faces[i][0], x_faces[i][1], x_faces[i][2]);
          for(int j=0; j<nTargets; j++)
          {
            fprintf(fp, " %+.10le %+.10le", targets_interp[i][j], values[i][j]);
          }
          fprintf(fp, "\n");
        }
        #else
        for(int i=0; i<nFaces; i++)
        {
          fprintf(fp, "%+.10le %+.10le %+.10le", x_faces[i][0], x_faces[i][1], x_faces[i][2]);
          for(int j=0; j<nTargets; j++)
          {
            fprintf(fp, " %+.10le %+.10le", targets_interp[i][j], values[i][j].value());
          }
          fprintf(fp, "\n");
        }
        #endif
        fclose(fp);
      SERIAL_END
    }
};
