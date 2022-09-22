#ifndef FIML_AUG_GAMMA_H
#define FIML_AUG_GAMMA_H

#include "BinningLearner.h"

typedef  double  vec3_t[3];
typedef adouble avec3_t[3];

class FIML_Aug_Gamma
{
  public:

    int *ncv;
    int rank, size;

    double *psi, *dRdbeta;

    double beta;
    adouble beta_AD;
    
    std::vector<double>  features;
    std::vector<double>  dbeta_dftrs;
    std::vector<adouble> features_AD;
    
    std::string    name;
    BinningLearner learner;

    double **rho, **rhoE, **strMag, **vortMag, **wallDist, **mu, **ReThetaT, **mu_T_FS;
    vec3_t **rhou;
    std::vector<double**> scal;
    std::vector<vec3_t**> grad_scal;

    adouble **rho_AD, **rhoE_AD, **strMag_AD, **vortMag_AD, **mu_AD;
    avec3_t **rhou_AD;
    std::vector<adouble**> scal_AD;
    std::vector<avec3_t**> grad_scal_AD;

  public:

    FIML_Aug_Gamma(void)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if(rank==0) printf("FIML_Aug_Gamma()\n");
    }

    void init(std::string name_val)
    {
      name    = name_val;

      features.resize(3);
      features_AD.resize(3);

      dbeta_dftrs.resize(3);
      
      if(rank==0) printf("Initializing Learner\n");
      learner.init({30, 10, 10}, {3.0001, 1.0001, 1.0001}, {0.0, 0.0, 0.0});
      learner.readParams(name+"_params.dat");
      if(rank==0) printf("Learner initialized\n");
    }

    void evalFeatures(int icv)
    {
      double kine_cv   = (*(scal[0]))[icv];
      double omega_cv  = (*(scal[1]))[icv];
      double gammai_cv = (*(scal[2]))[icv];

			double ReThetaT_cv = (*ReThetaT)[icv];
      
      double U_cv  = 0.0;
      
      for(int i=0; i<3; i++)
        U_cv += pow((*rhou)[icv][i],2);
      
      U_cv = sqrt(U_cv)/(*rho)[icv];

      double d_cv    = (*wallDist)[icv];
      double vort_cv = (*vortMag)[icv];
      double nu_cv   = (*mu)[icv] / (*rho)[icv];
      
      features[0] = fmin(d_cv * d_cv * vort_cv / (2.188 * nu_cv * ReThetaT_cv), 3.0);
      features[1] = d_cv / (d_cv + sqrt(kine_cv)/omega_cv);
			features[2] = nu_cv / (nu_cv + kine_cv/omega_cv);
    }

    void evalFeatures_AD(int icv)
    {
      adouble kine_cv   = (*(scal_AD[0]))[icv];
      adouble omega_cv  = (*(scal_AD[1]))[icv];
      adouble gammai_cv = (*(scal_AD[2]))[icv];

			double ReThetaT_cv = (*ReThetaT)[icv];
      
      adouble U_cv  = 0.0;
      
      for(int i=0; i<3; i++)
        U_cv += pow((*rhou_AD)[icv][i],2);
      
      U_cv = sqrt(U_cv)/(*rho_AD)[icv];

      adouble d_cv    = (*wallDist)[icv];
      adouble vort_cv = (*vortMag_AD)[icv];
      adouble nu_cv   = (*mu_AD)[icv] / (*rho_AD)[icv];
      
      features_AD[0] = fmin(d_cv * d_cv * vort_cv / (2.188 * nu_cv * ReThetaT_cv), 3.0);
      features_AD[1] = d_cv / (d_cv + sqrt(kine_cv)/omega_cv);
			features_AD[2] = nu_cv / (nu_cv + kine_cv/omega_cv);
      
      features[0] = features_AD[0].value();
      features[1] = features_AD[1].value();
      features[2] = features_AD[2].value();
    }

    void calculate(int icv)
    {
      evalFeatures(icv);
      beta = learner.calculate(features, dbeta_dftrs);
    }

    void calculate_AD(int icv)
    {
      evalFeatures_AD(icv);
      beta = learner.calculate(features, dbeta_dftrs);
      beta_AD = beta;
      for(int iFtr=0; iFtr<features.size(); iFtr++)
      {
        beta_AD = beta_AD + (features_AD[iFtr]-features[iFtr]) * dbeta_dftrs[iFtr];
      }
    }

    void writeDerivs(void)
    {
      learner.zerofyDerivs();
      for(int icv=0; icv<(*ncv); icv++)
      {
        evalFeatures(icv);
        learner.calcDerivs(psi[icv]*dRdbeta[icv], features);
      }
      learner.writeDerivs(name+"_derivs.dat");
    }
};

#endif
