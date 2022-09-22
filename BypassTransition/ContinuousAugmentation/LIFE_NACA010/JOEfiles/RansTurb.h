#ifndef RANSTURB_H
#define RANSTURB_H

#include "UgpWithCvCompFlow.h"
#include "LIFE/AugContainer.hpp"
#include "AUG_gamma_max.hpp"
#include "AUG_gamma_sep.hpp"
#include "WallLinker/WallLinker.h"

class RansTurb : virtual public UgpWithCvCompFlow
{
  public:

    double *muT, *muLam, *wallDist, *mu_T_FS, *ReThetaT_local, *ReThetaT;
    double (*grad_kine)[3];
    double *omega, (*grad_omega)[3];
    double *gammai, (*grad_gammai)[3];
    double *testVar;
    AugVarContainer augVars;
    GammaMaxAug gamma_max;
    GammaSepAug gamma_sep;
	WallLinker wLink;
	#include "WallLinker/WallLinkerSubroutines.h"

  public:

		virtual ~RansTurb(void) { deleteWallLinker(&wLink); }

    RansTurb(void)
    {
      if(mpi_rank==0) std::cout<<"RansTurb()"<<std::endl;
      turbModel = KOMGAM;
      registerScalar(gamma_max.dR_dBeta, "gamma_max_dR_dBeta", CV_DATA);
      registerScalar(gamma_sep.dR_dBeta, "gamma_sep_dR_dBeta", CV_DATA);

      registerDefaultVariables();
      registerModelEquations();

      augVars.scal.resize(3);
      augVars.grad_scal.resize(3);
      augVars.ncv          = &ncv;
      augVars.rho          = &rho;
      augVars.rhou         = &rhou;
      augVars.rhoE         = &rhoE;
      augVars.strMag       = &strMag;
      augVars.vortMag      = &vortMag;
      augVars.wallDist     = &wallDist;
      augVars.mu           = &muLam;
			augVars.ReThetaT     = &ReThetaT;
      augVars.scal[0]      = &(getScalarTransportData("kine")->phi);
      augVars.scal[1]      = &(getScalarTransportData("omega")->phi);
      augVars.scal[2]      = &(getScalarTransportData("gammai")->phi);
      augVars.grad_scal[0] = &(getScalarTransportData("kine")->grad_phi);
      augVars.grad_scal[1] = &(getScalarTransportData("omega")->grad_phi);
      augVars.grad_scal[2] = &(getScalarTransportData("gammai")->grad_phi);

      gamma_max.init(&augVars);
      gamma_sep.init(&augVars);
    }

    void registerDefaultVariables(void)
    {
      strMag         = NULL; registerScalar(strMag,         "strMag",           CV_DATA);
      vortMag        = NULL; registerScalar(vortMag,        "vortMag",          CV_DATA);
      diverg         = NULL; registerScalar(diverg,         "diverg",           CV_DATA);
      muT            = NULL; registerScalar(muT,            "muT",              CV_DATA);
      muLam          = NULL; registerScalar(muLam,          "muLam",            CV_DATA);
      wallDist       = NULL; registerScalar(wallDist,       "wallDist",         CV_DATA);
      testVar        = NULL; registerScalar(testVar,        "testVar",          CV_DATA);
      mu_T_FS        = NULL; registerScalar(mu_T_FS,        "mu_T_FS",          CV_DATA);
      ReThetaT       = NULL; registerScalar(ReThetaT,       "ReThetaT",         CV_DATA);
      ReThetaT_local = NULL; registerScalar(ReThetaT_local, "ReThetaT_local",   CV_DATA);
    }

    void registerModelEquations(void)
    {
      registerModelEquation("kine", 1e-10, 1e10);
      registerModelEquation("omega", 1e-4, 1e15);
      registerModelEquation("gammai", 1e-12, 2.0);
    }

    void registerModelEquation
    (
      const std::string &name,
      const double &lb,
      const double &ub
    )
    {
      ScalarTranspEq *eq = registerScalarTransport(name.c_str(), CV_DATA);
      eq->relax      = getDoubleParam(("RELAX_"+name).c_str(),"0.4");
      eq->phiZero    = 1e-8;
      eq->phiZeroRel = 1e-2;
      eq->phiMaxiter = 500;
      eq->lowerBound = lb;
      eq->upperBound = ub;
    }

    void initialHookScalarRansTurbModel(void)
    {
      if(mpi_rank==0) printf("initialHookScalarRansTurbModel()\n");
    
      initializeWallDistance();
      updateCvDataG1G2(wallDist, REPLACE_DATA);
    
      if(mpi_rank==0) printf("Initializing scalar variables and gradients\n");
      initializeModelEquations();
      if(mpi_rank==0) printf("initialHookScalarRansTurbModel() done\n");
    }

    void initializeWallDistance(void)
    {
      if(checkParam("NO_WALLDIST"))
      {
        for(int icv=0; icv<ncv_gg; icv++) wallDist[icv] = 1e5;
      }
      else
      {
        for(int icv=0; icv<ncv; icv++) wallDist[icv] = 0.0;
				if(checkParam("USE_FS"))
				{
					initWallLinker(&wLink);
					calcWallDistanceModified(&wLink, wallDist);
				}
				else
					calcWallDistance(NULL, wallDist);
      }
    }

    void initializeModelEquations(void)
    {
      initializeModelEquation("kine",   &kine,   &grad_kine  );
      initializeModelEquation("omega",  &omega,  &grad_omega );
      initializeModelEquation("gammai", &gammai, &grad_gammai);
    }
    
    void initializeModelEquation
    (
      std::string name,
      double **var,
      double (**grad_var)[3]
    )
    {
      Param *pmy;
      bool chk = getParam(pmy, name+"_INITIAL");
      if(!chk && mpi_rank==0){ printf(("Could not find "+name+"_INITIAL").c_str()); throw(-1); }
      if(mpi_rank==0) printf(("Initializing scalar "+name).c_str());
      double value = pmy->getDouble(1);
      ScalarTranspEq *eq = getScalarTransportData(name.c_str());
      (*var) = eq->phi;
      (*grad_var) = eq->grad_phi;
      if(!checkScalarFlag((char*)(name.c_str())))
        for(int icv=0; icv<ncv; icv++) (*var)[icv] = value;
    }

    void getInterpolationFactors
    (
      int ifa,
      int &icv0,
      int &icv1,
      double &w0,
      double &w1
    )
    {
      icv0 = cvofa[ifa][0];
      icv1 = cvofa[ifa][1];
      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      w0 = sqrt(vecDotVec3d(dx0, dx0));
      w1 = sqrt(vecDotVec3d(dx1, dx1));
      double ws = w0 + w1; w0 /= ws; w1 /= ws;
    }

    void calcRansTurbViscMuet(void)
    {
      calcGradVel(); calcStrainRateAndDivergence(); calcVorticity();

      for(int ifa=nfa_b; ifa<nfa; ifa++)
      {
        int icv0, icv1; double w0, w1;
        getInterpolationFactors(ifa, icv0, icv1, w0, w1);
        
        double rho_fa  = w1*rho[icv0]   + w0*rho[icv1];
        double kine_fa = w1*kine[icv0]  + w0*kine[icv1];
        double om_fa   = w1*omega[icv0] + w0*omega[icv1];
        
        mut_fa[ifa] = rho_fa*kine_fa/om_fa;
      }

      for(list<FaZone>::iterator ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
      if(ZZ->getKind() == FA_ZONE_BOUNDARY)
      {
        if(zoneIsWall(ZZ->getName()))
        for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
        {
          mut_fa[ifa] = 0.0;        
        }
        
        if(!zoneIsWall(ZZ->getName()))
        for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
        {
          int icv1 = cvofa[ifa][1];
          mut_fa[ifa] = rho[icv1]*kine[icv1]/omega[icv1];
        }
      }
    
      for (int icv=0; icv<ncv; icv++)
      {
        muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
      }
    }

    void diffusivityHookScalarRansTurb(const std::string &name)
    {
      ScalarTranspEq *eq; double sigma, sigmal;
      
      if((name == "kine") || (name == "omega") || (name =="gammai"))
      {
        eq = getScalarTransportData(name);
        
        if(name == "kine")   { sigma = 0.5; sigmal = 1.0; }
        if(name == "omega")  { sigma = 0.5; sigmal = 1.0; }
        if(name == "gammai") { sigma = 5.0; sigmal = 5.0; }
    
        for (int ifa=nfa_b; ifa<nfa; ifa++)
        {
          double w0, w1; int icv0, icv1;
          getInterpolationFactors(ifa, icv0, icv1, w0, w1);
          
          double rho_fa  = w1*rho[icv0]   + w0*rho[icv1];
          double kine_fa = w1*kine[icv0]  + w0*kine[icv1];
          double om_fa   = w1*omega[icv0] + w0*omega[icv1];
          
          eq->diff[ifa] = mul_fa[ifa]/sigmal + sigma*rho_fa*kine_fa/om_fa;
        }
    
        for(list<FaZone>::iterator ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
        if(ZZ->getKind() == FA_ZONE_BOUNDARY)
        {
          if(zoneIsWall(ZZ->getName()))
          for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa]/sigmal;
          }
          
          if(!zoneIsWall(ZZ->getName()))
          for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
          {
            int icv1 = cvofa[ifa][1];
            eq->diff[ifa] = mul_fa[ifa]/sigmal + sigma*rho[icv1]*kine[icv1]/omega[icv1];
          }
        }
      }
    }

    void boundaryHookScalarRansTurb
    (
      double *phi,
      FaZone *zone,
      const std::string &name
    )
    {
      if(name=="omega")
      if(zone->getKind()==FA_ZONE_BOUNDARY)
      if(zoneIsWall(zone->getName()))
      for(int index=0; index<zone->faVec.size(); index++)
      {
        int ifa = zone->faVec[index];
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];
        double mu_cv = calcMuLam(icv0);
        phi[icv1] = 6.0*mu_cv/(rho[icv0]*0.075*pow(wallDist[icv0],2));
      }
    }
    
    void sourceHookRansTurbCoupled
    (
      double **rhs,
      double ***A,
      int nScal,
      int flagImplicit
    )
    {
      int noc00, index; double rhoVol, src, d_src, dRdbeta, nuT, nuL, Prod;
    
      for(int icv=0; icv<ncv; icv++)
      {
        rhoVol = rho[icv] * cv_volume[icv];

        noc00 = nbocv_i[icv];
        nuT = kine[icv]/omega[icv];
        nuL = calcMuLam(icv)/rho[icv];
        muLam[icv] = nuL * rho[icv];
        Prod = fmax(nuT*vortMag[icv]*vortMag[icv] - (2./3.)*kine[icv]*diverg[icv], 0.0);
    
        // kine :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        index = getScalarTransportIndex("kine");
        src   = fmin(2.0, gammai[icv]) * Prod - 0.09 * kine[icv] * omega[icv];
        d_src = -0.09 * omega[icv];
    
        rhs[icv][5+index] += rhoVol * src;
        if(flagImplicit) A[noc00][5+index][5+index] -= d_src * cv_volume[icv];
    
        // omega ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        index = getScalarTransportIndex("omega");
        src   = (5./9.) * omega[icv] / kine[icv] * Prod - 0.075 * omega[icv] * omega[icv];
        d_src = -0.15 * omega[icv];
    
        rhs[icv][5+index] += rhoVol * src; 
        if(flagImplicit) A[noc00][5+index][5+index] -= d_src * cv_volume[icv];

        // gammai :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        index = getScalarTransportIndex("gammai");
        gamma_sep.calculate(icv); double bfSep = getBlendingFactor(icv);
        gamma_max.calculate(icv); double bfMax = 1; //1-bfSep;

        src   = (pow(gamma_max.beta, bfMax) * pow(gamma_sep.beta, bfSep) - gammai[icv])
              * vortMag[icv] * sqrt(gammai[icv]);
        d_src = -1.5 * vortMag[icv] * sqrt(gammai[icv]);
    
        rhs[icv][5+index] += rhoVol * src;
        if(flagImplicit) A[noc00][5+index][5+index] -= d_src * cv_volume[icv];

        gamma_max.dR_dBeta[icv] = rhoVol * vortMag[icv] * sqrt(gammai[icv])
                                * bfMax * pow(gamma_max.beta, bfMax-1) * pow(gamma_sep.beta, bfSep);

        gamma_sep.dR_dBeta[icv] = rhoVol * vortMag[icv] * sqrt(gammai[icv])
                                * bfSep * pow(gamma_max.beta, bfMax) * pow(gamma_sep.beta, bfSep-1);

        // ReThetaT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double Umag = sqrt(pow(rhou[icv][0],2)+pow(rhou[icv][1],2)+pow(rhou[icv][2],2))/rho[icv];
				ReThetaT_local[icv] = Umag * sqrt((7.0*nuL)/(9.0*omega[icv])) / nuL;
      }
    }

    double sigmoid(double x, double s)
    {
      return 1.0 / (1.0 + exp(-x/s));
    }

    void writeDerivsFIML(void) { gamma_max.writeSens(); gamma_sep.writeSens(); }

    void writeFeatures(void)
    {
      SERIAL_BEG
        FILE* fp;
        if(mpi_rank==0) fp = fopen("features.out","w");
        else fp = fopen("features.out","a");
        for(int icv=0; icv<ncv; icv++)
        {
          double ftr0 = fmin(pow(wallDist[icv],2) * vortMag[icv] / (2.188 * muLam[icv]/rho[icv] * ReThetaT[icv]), 3.0);
          double ftr1 = wallDist[icv] / (wallDist[icv] + sqrt(kine[icv])/omega[icv]);
          double ftr2 = muLam[icv] / (muLam[icv] + rho[icv]*kine[icv]/omega[icv]);
          fprintf(fp, "%9d    %+.10le    %+.10le    %+.10le    %+.10le    %+.10le\n", (int)(global_id[icv]), x_cv[icv][0], x_cv[icv][1], ftr0, ftr1, ftr2);
        }
        fclose(fp);
      SERIAL_END
    }

    void calcWallNormal(int icv, double *wn)
    {
      //--------------------------------------------------------------------------------------
      // Calculate gradient of wallDist (Replace fake cells with internal values)
      //--------------------------------------------------------------------------------------
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;

      for(int i=0; i<3; i++) wn[i] = 0.0;

      for(int noc = noc_f; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_all_v[noc]; if(icv_nbr>=ncv_gg) icv_nbr = icv;
        for(int i=0; i<3; i++)
          wn[i] = wn[i] + nbocv_all_grad_coeff[noc][i] * wallDist[icv_nbr];
      }

      //--------------------------------------------------------------------------------------
      // Normalize gradient vector to get approximate wall normal
      //--------------------------------------------------------------------------------------
      double mag = 0.0;
      for(int i=0; i<3; i++) mag += wn[i] * wn[i];
      mag = sqrt(mag);
      
      for(int i=0; i<3; i++) wn[i] /= mag;
    }

    //double getBlendingFactor(int icv)
    //{
    //  double wn[3]; calcWallNormal(icv, wn);
    //  // Directional derivative of vorticity magnitude in the normal direction
    //  double dvn = 0.0; 
    //  int noc_f = nbocv_all_i[icv];
    //  FOR_I3 dvn = nbocv_all_grad_coeff[noc_f][i] * vortMag[icv] * wn[i];
    //  int noc_l = nbocv_all_i[icv+1]-1;
    //  for(int noc = noc_f + 1; noc <= noc_l; noc++)
    //  {
    //    int icv_nbr = nbocv_all_v[noc];
    //    double nbrVort = (icv_nbr >= ncv_gg) ? vortMag[icv] : vortMag[icv_nbr];
    //    FOR_I3 dvn = dvn + nbocv_all_grad_coeff[noc][i] * nbrVort * wn[i];
    //  }
    //  double vort_by_dist = vortMag[icv] / sqrt(muLam[icv] / rho[icv] / vortMag[icv]);
    //  Param *swtch = getParam("SWITCH_AUG_AT");
    //  double result = dvn / (fabs(dvn) + vort_by_dist) * omega[icv]/(sqrt(2)*vortMag[icv]+omega[icv]);
    //  return sigmoid(  result - swtch->getDouble(1), 0.003) *
    //         sigmoid(- result + swtch->getDouble(2), 0.003);
    //}

    double getBlendingFactor(int icv)
    {
      //--------------------------------------------------------------------------------------
      // Calculate wall normal
      //--------------------------------------------------------------------------------------
      double wn[3];
      calcWallNormal(icv, wn);
      
      //--------------------------------------------------------------------------------------
      // Calculate wall-normal derivative of vortMag (Replace fake cells with internal values)
      //--------------------------------------------------------------------------------------
      int noc_f = nbocv_all_i[icv];
      int noc_l = nbocv_all_i[icv+1]-1;

      double dvn = 0.0; 

      for(int noc = noc_f; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_all_v[noc]; if(icv_nbr>=ncv_gg) icv_nbr = icv;
        for(int i=0; i<3; i++)
          dvn = dvn + nbocv_all_grad_coeff[noc][i] * vortMag[icv_nbr] * wn[i];
      }

      //--------------------------------------------------------------------------------------
      // Calculate f_sigma
      //--------------------------------------------------------------------------------------
      double nu      = muLam[icv] / rho[icv];
      double invlen  = sqrt(vortMag[icv] / nu);
      double fun1    = dvn / (fabs(dvn) + vortMag[icv]*invlen);
      double fun2    = omega[icv] / (sqrt(2)*vortMag[icv] + omega[icv]);
      double f_sigma = fun1 * fun2;

      //--------------------------------------------------------------------------------------
      // Calculate sigma
      //--------------------------------------------------------------------------------------
      Param *swtch = getParam("SWITCH_AUG_AT");
      double lb    = swtch->getDouble(1);
      double ub    = swtch->getDouble(2);
      double s     = 0.003;
      double sigma = sigmoid(f_sigma-lb, s) * sigmoid(ub-f_sigma, s);
      return sigma;
    }
};

#endif
