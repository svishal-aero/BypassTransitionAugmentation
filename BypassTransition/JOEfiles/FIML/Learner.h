#ifndef LEARNER_H
#define LEARNER_H

#include "adolc.h"
#include <vector>
#include <string>

class Learner
{
  public:

    std::vector<double> params;
    std::vector<double> derivs;

  public:

    Learner(void) { }
    
    virtual ~Learner() { }

    virtual    void initParamsAndDerivs(void) { return; }
    virtual  double calculate(std::vector<double> features) { return 0.0; }
    virtual adouble calculate_AD(std::vector<adouble> features_AD) { return 0.0; }
    virtual    void calcDerivs(double dJdbeta, std::vector<double> features) { return; }
    
    virtual void addGridSearchVar(double min, double max, int nBins) { return; }

    void zerofyDerivs(void) { for(int i=0; i<derivs.size(); i++) derivs[i] = 0.0; }

    void readParams(std::string filename)
    {
      int rank, size;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if(rank>0) MPI_Recv(NULL, 0, MPI_CHAR, rank-1, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      FILE *fp = fopen(filename.c_str(), "r");
      
      for(int i=0; i<params.size(); i++)
      {
        int rtn = fscanf(fp, "%le", &(params[i]));
      }
      
      fclose(fp);
      
      if(rank<size-1) MPI_Send(NULL, 0, MPI_CHAR, rank+1, 123, MPI_COMM_WORLD);
    }

    void writeDerivs(std::string filename)
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      std::vector<double> derivs_temp = derivs;
      MPI_Allreduce(&(derivs_temp[0]), &(derivs[0]), derivs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(rank==0)
      {
        FILE *fp = fopen(filename.c_str(), "w");

        for(int i=0; i<derivs.size(); i++)
        {
          fprintf(fp, "%+.12le\n", derivs[i]);
        }

        fclose(fp);
      }
    }
};

#endif
