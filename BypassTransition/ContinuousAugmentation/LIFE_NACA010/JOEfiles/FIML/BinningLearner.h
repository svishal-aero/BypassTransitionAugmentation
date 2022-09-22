#ifndef BinningLearner_h
#define BinningLearner_h

#include <cstdio>
#include <iostream>
#include <vector>
#include <cassert>

class ParamIndexer
{
  public:

    int     nFtrs;
    int    *nBins;
    double *params;
    double *derivs;

  public:

    void init
    (
      int nFtrs,
      std::vector<int>& nBins,
      std::vector<double>& params,
      std::vector<double>& derivs
    )
    {
      this->nFtrs = nFtrs;
      this->nBins = &(nBins[0]);
      this->params = &(params[0]);
      this->derivs = &(derivs[0]);
    }

    double getParam(std::vector<int>& indices)
    {
      int ind = 0;
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        ind *= nBins[iFtr];
        ind += indices[iFtr];
      }
      return params[ind];
    }

    void setDeriv(std::vector<int>& indices, double deriv_val)
    {
      int ind = 0;
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        ind *= nBins[iFtr];
        ind += indices[iFtr];
      }
      derivs[ind] = deriv_val;
    }
};

class BinningLearner
{
  public:
    
    double beta;
    std::vector<double> params;
    std::vector<double> derivs;

  private:

    int nFtrs;
    ParamIndexer indxr;
    std::vector<int> nBins;
    std::vector<int> indices;
    std::vector<double> max_ftrs;
    std::vector<double> min_ftrs;

  public:
    
    void init
    (
      std::vector<int> nBins,
      std::vector<double> max_ftrs,
      std::vector<double> min_ftrs
    )
    {
      this->nBins = nBins;
      this->max_ftrs = max_ftrs;
      this->min_ftrs = min_ftrs;
      this->nFtrs = this->nBins.size();
      this->indices.resize(nFtrs, 0);
      this->params.resize(get_nParams(), 1.0);
      this->derivs.resize(get_nParams(), 0.0);
      this->indxr.init(nFtrs, nBins, params, derivs);
    }

    int get_nParams()
    {
      int nParams = 1;
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        nParams *= nBins[iFtr];
      }
      return nParams;
    }

    void readParams(std::string filename)
    {
      SERIAL_BEG;
        FILE *fp = fopen(filename.c_str(), "r");
        for(int iParam=0; iParam<params.size(); iParam++)
        {
          int dummy = fscanf(fp, "%le", &(params[iParam]));
        }
        fclose(fp);
      SERIAL_END;
    }

    void writeDerivs(std::string filename)
    {
      std::vector<double> derivs_total(derivs.size(), 0.0);
      MPI_Allreduce(&(derivs[0]), &(derivs_total[0]), derivs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      if(mpi_rank==0)
      {
        FILE *fp = fopen(filename.c_str(), "w");
        for(int iParam=0; iParam<params.size(); iParam++)
        {
          fprintf(fp, "%+.15le\n", derivs_total[iParam]);
        }
        fclose(fp);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    double calculate
    (
      std::vector<double>& ftrs,
      std::vector<double>& dbeta_dftrs
    )
    {
      int ind; double lval, cval, rval, dx, d_ftr;
      
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        d_ftr = (max_ftrs[iFtr]-min_ftrs[iFtr]) / nBins[iFtr];
        indices[iFtr] = std::lround(floor((ftrs[iFtr]-min_ftrs[iFtr])/d_ftr));
      }

      ind = 0;
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        ind *= nBins[iFtr];
        ind += indices[iFtr];
      }
      cval = params[ind];
      beta = cval;
      
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        d_ftr = (max_ftrs[iFtr]-min_ftrs[iFtr]) / nBins[iFtr];
        
        if(indices[iFtr]==0)
        {
          lval = cval;
        }
        else
        {
          indices[iFtr] -= 1;
          lval = 0.5 * cval;
          ind = 0;
          for(int iFtr=0; iFtr<nFtrs; iFtr++)
          {
            ind *= nBins[iFtr];
            ind += indices[iFtr];
          }
          lval += 0.5 * params[ind];
          indices[iFtr] += 1;
        }
        
        if(indices[iFtr]==nBins[iFtr]-1)
        {
          rval = cval;
        }
        else 
        {
          indices[iFtr] += 1;
          rval = 0.5 * cval;
          ind = 0;
          for(int iFtr=0; iFtr<nFtrs; iFtr++)
          {
            ind *= nBins[iFtr];
            ind += indices[iFtr];
          }
          rval += 0.5 * params[ind];
          indices[iFtr] -= 1;
        }
        
        dx = (ftrs[iFtr] - (indices[iFtr]+0.5)*d_ftr);
        dbeta_dftrs[iFtr] = (rval - lval) / d_ftr;
        beta += dbeta_dftrs[iFtr] * dx;
      }

      return beta;
    }

    void calcDerivs(double dJdbeta, std::vector<double>& ftrs)
    {
      int ind; double lval, cval, rval, dx, d_ftr, dJ_dlval, dJ_drval, dJ_dcval;
      
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        d_ftr = (max_ftrs[iFtr]-min_ftrs[iFtr]) / nBins[iFtr];
        indices[iFtr] = std::lround(floor((ftrs[iFtr]-min_ftrs[iFtr])/d_ftr));
      }

      dJ_dcval = 0.0;
      
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        d_ftr = (max_ftrs[iFtr]-min_ftrs[iFtr]) / nBins[iFtr];
        dx = (ftrs[iFtr] - (indices[iFtr]+0.5)*d_ftr);
        dJ_drval =  dx / d_ftr;
        dJ_dlval = -dx / d_ftr;

        if(indices[iFtr]==0)
        {
          dJ_dcval += dJ_dlval;
        }
        else
        {
          indices[iFtr] -= 1;
          dJ_dcval += 0.5 * dJ_dlval;
          ind = 0;
          for(int iFtr=0; iFtr<nFtrs; iFtr++)
          {
            ind *= nBins[iFtr];
            ind += indices[iFtr];
          }
          derivs[ind] += 0.5 * dJ_dlval * dJdbeta;
          indices[iFtr] += 1;
        }
        
        if(indices[iFtr]==nBins[iFtr]-1)
        {
          dJ_dcval += dJ_drval;
        }
        else 
        {
          indices[iFtr] += 1;
          dJ_dcval += 0.5 * dJ_drval;
          ind = 0;
          for(int iFtr=0; iFtr<nFtrs; iFtr++)
          {
            ind *= nBins[iFtr];
            ind += indices[iFtr];
          }
          derivs[ind] += 0.5 * dJ_drval * dJdbeta;
          indices[iFtr] -= 1;
        }
      }

      ind = 0;
      for(int iFtr=0; iFtr<nFtrs; iFtr++)
      {
        ind *= nBins[iFtr];
        ind += indices[iFtr];
      }
      derivs[ind] += dJ_dcval * dJdbeta;
    }

    void zerofyDerivs(void)
    {
      for(int iParam=0; iParam<params.size(); iParam++) derivs[iParam] = 0.0;
    }
};

#endif
