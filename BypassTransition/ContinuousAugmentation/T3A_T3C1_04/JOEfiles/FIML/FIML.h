#ifndef FIML_AUG_H
#define FIML_AUG_H

#include "Learner.h"
#include "GridSearchLearner.h"

typedef  double  vec3_t[3];
typedef adouble avec3_t[3];

class FIML_Aug
{
  public:

    int *ncv;
    int rank, size;

    double *psi, *dRdbeta;

    double beta;
    adouble beta_AD;
    
    std::vector<double>  features;
    std::vector<adouble> features_AD;
    
    std::string  name;
    Learner     *learner;

    double **rho, **rhoE, **strMag, **vortMag, **wallDist, **mu;
    vec3_t **rhou;
    std::vector<double**> scal;
    std::vector<vec3_t**> grad_scal;

    adouble **rho_AD, **rhoE_AD, **strMag_AD, **vortMag_AD, **mu_AD;
    avec3_t **rhou_AD;
    std::vector<adouble**> scal_AD;
    std::vector<avec3_t**> grad_scal_AD;

  public:

    FIML_Aug(void)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if(rank==0) printf("FIML_Aug()\n");
      learner = nullptr;
    }

    void init(std::string name_val)
    {
      name    = name_val;

      features.resize(3);
      features_AD.resize(3);
      
      if(rank==0) printf("Initializing Learner\n");
      learner = new GridSearchLearner();
      learner->addGridSearchVar(0.0, 1.0, 30);
      learner->addGridSearchVar(0.0, 1.0, 30);
      learner->addGridSearchVar(0.0, 1.0, 30);
      //learner->addGridSearchVar(0.0, 1.0, 30);
      learner->initParamsAndDerivs();
      learner->readParams(name+"_params.dat");
      if(rank==0) printf("Learner initialized\n");
    }

    ~FIML_Aug(void) { if(learner!=nullptr) delete learner; }

    void evalFeatures(int icv)
    {
      double kine_cv   = (*(scal[0]))[icv];
      double omega_cv  = (*(scal[1]))[icv];
      double gammai_cv = (*(scal[2]))[icv];
      double thetaT_cv = (*(scal[3]))[icv];
      
      double U_cv  = 0.0;
      
      for(int i=0; i<3; i++)
        U_cv += pow((*rhou)[icv][i],2);
      
      U_cv = sqrt(U_cv)/(*rho)[icv];

      double d_cv    = (*wallDist)[icv];
      double vort_cv = (*vortMag)[icv];
      double nu_cv   = (*mu)[icv] / (*rho)[icv];
      
      features[0] = U_cv / (U_cv + d_cv * vort_cv);
      features[1] = thetaT_cv / (d_cv + thetaT_cv);
      features[2] = nu_cv / (kine_cv / omega_cv + nu_cv);
      //features[3] = gammai_cv;
    }

    void evalFeatures_AD(int icv)
    {
      adouble kine_cv   = (*(scal_AD[0]))[icv];
      adouble omega_cv  = (*(scal_AD[1]))[icv];
      adouble gammai_cv = (*(scal_AD[2]))[icv];
      adouble thetaT_cv = (*(scal_AD[3]))[icv];
      
      adouble U_cv  = 0.0;
      
      for(int i=0; i<3; i++)
        U_cv += pow((*rhou_AD)[icv][i],2);
      
      U_cv = sqrt(U_cv)/(*rho_AD)[icv];

      adouble d_cv    = (*wallDist)[icv];
      adouble vort_cv = (*vortMag_AD)[icv];
      adouble nu_cv   = (*mu_AD)[icv] / (*rho_AD)[icv];
      
      features_AD[0] = U_cv / (U_cv + d_cv * vort_cv);
      features_AD[1] = thetaT_cv / (d_cv + thetaT_cv);
      features_AD[2] = nu_cv / (kine_cv / omega_cv + nu_cv);
      //features_AD[3] = gammai_cv;
    }

    void calculate(int icv)
    {
      evalFeatures(icv);
      beta = learner->calculate(features);
    }

    void calculate_AD(int icv)
    {
      evalFeatures_AD(icv);
      beta_AD = learner->calculate_AD(features_AD);
    }

    void writeDerivs(void)
    {
      learner->zerofyDerivs();
      for(int icv=0; icv<(*ncv); icv++)
      {
        evalFeatures(icv);
        learner->calcDerivs(psi[icv]*dRdbeta[icv], features);
      }
      learner->writeDerivs(name+"_derivs.dat");

      if(rank>0) MPI_Recv(NULL, 0, MPI_CHAR, rank-1, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      FILE *fp;
      if(rank==0) fp = fopen("verify_adjoints.dat", "w");
      if(rank!=0) fp = fopen("verify_adjoints.dat", "a");
      
      for(int icv=0; icv<(*ncv); icv++)
      {
        int rtn = fprintf(fp, "%+.12le    %+.12le\n", psi[icv], dRdbeta[icv]);
      }
      
      fclose(fp);
      
      if(rank<size-1) MPI_Send(NULL, 0, MPI_CHAR, rank+1, 123, MPI_COMM_WORLD);
    }
};

#endif
