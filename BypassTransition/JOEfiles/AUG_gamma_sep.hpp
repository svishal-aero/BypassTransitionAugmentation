#ifndef AUG_gamma_sep_hpp
#define AUG_gamma_sep_hpp

#include "LIFE/AugContainer.hpp"
#include "LIFE/AugMaps/PwMultiLinAugMap/PwMultiLinAugMap.hpp"

class GammaSepAug : public AugContainer<PwMultiLinAugMap<3>, 3>
{
  public:

    GammaSepAug() : AugContainer<PwMultiLinAugMap<3>, 3>() { }

    inline void init(
      AugVarContainer* vars_
    )
    {
      this->vars = vars_;
      this->augMap.init("gamma_sep", {51, 11, 11}, {0, 0, 0},
                        {5+1e-5, 1+1e-5, 1+1e-5});
      this->augMap.readParams();
    }

    inline void evalFtrs(
      int icv
    )
    {
      double     kine_cv = this->vars->scal[0][0][icv];
      double    omega_cv = this->vars->scal[1][0][icv];
      double    gamma_cv = this->vars->scal[2][0][icv];
      double ReThetaT_cv = this->vars->ReThetaT[0][icv];
      
      double U_cv = 0.0;
      for(int i=0; i<3; i++)
        U_cv += pow(this->vars->rhou[0][icv][i], 2);
      U_cv = sqrt(U_cv) / this->vars->rho[0][icv];

      double    d_cv = this->vars->wallDist[0][icv];
      double vort_cv = this->vars->vortMag[0][icv];
      double   nu_cv = this->vars->mu[0][icv] / this->vars->rho[0][icv];

      this->ftrs[0] = fmin(d_cv*d_cv*vort_cv/(2.188*nu_cv*ReThetaT_cv), 3);
      this->ftrs[1] = d_cv*omega_cv/(d_cv*omega_cv+sqrt(kine_cv));
      this->ftrs[2] = nu_cv*omega_cv/(nu_cv*omega_cv+kine_cv);
      //this->ftrs[3] = d_cv*vort_cv/(d_cv*vort_cv+U_cv);
    }

    inline void evalFtrs_AD(
      int icv
    )
    {
      adouble     kine_cv = this->vars->scal_AD[0][0][icv];
      adouble    omega_cv = this->vars->scal_AD[1][0][icv];
      adouble    gamma_cv = this->vars->scal_AD[2][0][icv];
      double  ReThetaT_cv = this->vars->ReThetaT[0][icv];
      
      adouble U_cv = 0.0;
      for(int i=0; i<3; i++)
        U_cv += pow(this->vars->rhou_AD[0][icv][i], 2);
      U_cv = sqrt(U_cv) / this->vars->rho_AD[0][icv];

      double     d_cv = this->vars->wallDist[0][icv];
      adouble vort_cv = this->vars->vortMag_AD[0][icv];
      adouble   nu_cv = this->vars->mu_AD[0][icv] / this->vars->rho_AD[0][icv];

      this->ftrs_AD[0] = fmin(d_cv*d_cv*vort_cv/(2.188*nu_cv*ReThetaT_cv), 3);
      this->ftrs_AD[1] = d_cv*omega_cv/(d_cv*omega_cv+sqrt(kine_cv));
      this->ftrs_AD[2] = nu_cv*omega_cv/(nu_cv*omega_cv+kine_cv);
      //this->ftrs_AD[3] = d_cv*vort_cv/(d_cv*vort_cv+U_cv);

      this->ftrs[0] = this->ftrs_AD[0].value();
      this->ftrs[1] = this->ftrs_AD[1].value();
      this->ftrs[2] = this->ftrs_AD[2].value();
      //this->ftrs[3] = this->ftrs_AD[3].value();
    }
};

#endif
