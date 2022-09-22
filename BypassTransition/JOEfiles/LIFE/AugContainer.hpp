#ifndef AugContainer_hpp
#define AugContainer_hpp

#include "AugMaps/GenAugMap.hpp"

typedef  double  vec3_t[3];
typedef adouble avec3_t[3];

struct AugVarContainer
{
  int *ncv;

  double **rho;
  vec3_t **rhou;
  double **rhoE;
  double **strMag;
  double **vortMag;
  double **mu;
  double **mu_T_FS;
  double **wallDist;
  double **ReThetaT;
  std::vector<double**> scal;
  std::vector<vec3_t**> grad_scal;

  adouble **rho_AD;
  avec3_t **rhou_AD;
  adouble **rhoE_AD;
  adouble **strMag_AD;
  adouble **vortMag_AD;
  adouble **mu_AD;
  std::vector<adouble**> scal_AD;
  std::vector<avec3_t**> grad_scal_AD;
};

template <class AugMap, const int nFtrs>
struct AugContainer
{
  AugMap augMap;
  double *psi;
  double *dR_dBeta;
  AugVarContainer* vars;

  double beta;
  adouble beta_AD;
  std::array<double, nFtrs> ftrs;
  std::array<adouble, nFtrs> ftrs_AD;
  std::array<double, nFtrs> dBeta_dFtrs;

  AugContainer(
    void
  )
  {
    psi      = NULL;
    dR_dBeta = NULL;
    vars     = NULL;
  }

  virtual void init(
    AugVarContainer* vars_
  ) = 0;

  virtual void evalFtrs(
    int icv
  ) = 0;

  virtual void evalFtrs_AD(
    int icv
  ) = 0;

  inline void calculate(
    int icv
  )
  { 
    evalFtrs(icv);
    beta = augMap.calculate_value(ftrs);
  }

  inline void calculate_AD(
    int icv
  )
  {
    evalFtrs_AD(icv);
    beta = augMap.calculate_value(ftrs, dBeta_dFtrs.data());
    beta_AD = beta;
    for(int iFtr=0; iFtr<nFtrs; iFtr++)
    {
      beta_AD = beta_AD +
        (ftrs_AD[iFtr] - ftrs[iFtr]) * dBeta_dFtrs[iFtr];
    }
  }

  inline void writeSens(
    void
  )
  {
    augMap.zerofySens();
    
    for(int icv=0; icv<vars->ncv[0]; icv++)
    {
      evalFtrs(icv);
      augMap.calculate_sens(ftrs, psi[icv]*dR_dBeta[icv]);
    }

    augMap.writeSens();
  }
};

#endif
