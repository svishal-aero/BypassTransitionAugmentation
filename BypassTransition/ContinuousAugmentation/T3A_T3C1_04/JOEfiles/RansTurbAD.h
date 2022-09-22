#ifndef RANSTURB_AD_H
#define RANSTURB_AD_H

#include "RansTurb.h"
#include "UgpWithCvCompFlowAD.h"

class RansTurb_AD : virtual public UgpWithCvCompFlow_AD, public RansTurb
{
  public:
    
    adouble (*grad_kine)[3], *kine_diff;
    adouble *omega, (*grad_omega)[3], *omega_diff;
    adouble *gammai, (*grad_gammai)[3], *gammai_diff;

  public:

    RansTurb_AD(void)
    {
      if(mpi_rank==0) std::cout<<"RansTurb_AD()"<<std::endl;
      registerScalar(gamma_max.psi, "gamma_max_psi", CV_DATA);
      
      augVars.scal_AD.resize(3);
      augVars.grad_scal_AD.resize(3);
      augVars.rho_AD          = &rho_AD;
      augVars.rhou_AD         = &rhou_AD;
      augVars.rhoE_AD         = &rhoE_AD;
      augVars.strMag_AD       = &strMag;
      augVars.vortMag_AD      = &vortMag;
      augVars.wallDist        = &wallDist;
      augVars.mu_AD           = &muLam_AD;
      augVars.scal_AD[0]      = &kine;
      augVars.scal_AD[1]      = &omega;
      augVars.scal_AD[2]      = &gammai;
      augVars.grad_scal_AD[0] = &grad_kine;
      augVars.grad_scal_AD[1] = &grad_omega;
      augVars.grad_scal_AD[2] = &grad_gammai;
    }

    void copy_turb_adjoint()
    {
      int nScal = scalarTranspEqVector.size();
      for(int i=0; i<nScal; i++)
      {
        if(strcmp(scalarTranspEqVector_AD[i].name,"gammai")==0)
        {
          for(int icv=0; icv<ncv_gg; icv++)
          {
            gamma_max.psi[icv] = scalarTranspEqVector_psi[i].phi[icv];
          }
        }
      }
    }

    void initialHookScalarRansTurbModel_AD(void)
    {
      if(mpi_rank==0) printf("initialHookScalarRansTurbModel()\n");

      connectPointers("kine",     &kine, &grad_kine,   &kine_diff);    
      connectPointers("omega",   &omega, &grad_omega,  &omega_diff);    
      connectPointers("gammai", &gammai, &grad_gammai, &gammai_diff);    

      updateCvDataG1G2(wallDist, REPLACE_DATA);
    }

    void initialize_turb_adjoint(void) { readAdjointsFromFile(); }

    void writeAdjointsToFile()
    {
      int ncv_global;
      MPI_Allreduce(&ncv, &ncv_global, 1, MPI_INT, MPI_SUM, mpi_comm);
      PRINT("Writing adjoints to file \n");
      PRINT("\tTotal number of cv's across all processors = %d\n", ncv_global);
      int nScal = scalarTranspEqVector.size();
      double *psi__kine, *psi__omega, *psi__gammai;
      for(int i=0; i<nScal; i++)
      {
        if(strcmp(scalarTranspEqVector_AD[i].name, "kine")==0)
        { psi__kine   = scalarTranspEqVector_psi[i].phi; PRINT("\tkine pointer connected\n"); }
        if(strcmp(scalarTranspEqVector_AD[i].name, "omega")==0)
        { psi__omega  = scalarTranspEqVector_psi[i].phi; PRINT("\tomega pointer connected\n"); }
        if(strcmp(scalarTranspEqVector_AD[i].name, "gammai")==0)
        { psi__gammai = scalarTranspEqVector_psi[i].phi; PRINT("\tgammai pointer connected\n"); }
      }
      SERIAL_BEG;
      FILE *fp;
      if(mpi_rank==0) fp = fopen("adjoint_values.txt", "w");
      else            fp = fopen("adjoint_values.txt", "a");
      for(int icv=0; icv<ncv; icv++)
      {
        fprintf(fp, "%09ld  %+.15le  %+.15le  %+.15le  %+.15le  %+.15le  %+.15le  %+.15le  %+.15le\n",
            std::lround(global_id[icv]),
            psi_rho[icv], psi_rhou[icv][0], psi_rhou[icv][1], psi_rhou[icv][2], psi_rhoE[icv],
            psi__kine[icv], psi__omega[icv], psi__gammai[icv]);
      }
      fclose(fp);
      SERIAL_END;
      MPI_Barrier(mpi_comm);
    }

    void readAdjointsFromFile()
    {
      PRINT("\n\nReading adjoints values from text file, if provided\nCreating local_id vector\n");

      int ncv_global;
      PRINT("Number of cv's on this processor = %d\n", ncv);
      MPI_Allreduce(&ncv, &ncv_global, 1, MPI_INT, MPI_SUM, mpi_comm);
      PRINT("Total number of cv's across processors = %d\n", ncv_global);
      std::vector<int> local_id(ncv_global, -1);
      for(int icv=0; icv<ncv; icv++)
        local_id[std::lround(global_id[icv])] = icv;

      int nScal = scalarTranspEqVector.size();
      double *psi__kine, *psi__omega, *psi__gammai;
      for(int i=0; i<nScal; i++)
      {
        if(strcmp(scalarTranspEqVector_AD[i].name, "kine")==0)
        { psi__kine   = scalarTranspEqVector_psi[i].phi; PRINT("kine pointer connected\n"); }
        if(strcmp(scalarTranspEqVector_AD[i].name, "omega")==0)
        { psi__omega  = scalarTranspEqVector_psi[i].phi; PRINT("omega pointer connected\n"); }
        if(strcmp(scalarTranspEqVector_AD[i].name, "gammai")==0)
        { psi__gammai = scalarTranspEqVector_psi[i].phi; PRINT("gammai pointer connected\n"); }
      }

      PRINT("Checking if file exists and reading...\n");
      
      SERIAL_BEG;
      printf("Reading on processor %d\n", mpi_rank);
      FILE *fp = fopen("adjoint_values.txt", "r");
      if(fp!=NULL)
      {
        int rtnval, gloId, locId;
        double temp;
        printf("\tReading line 000000000");
        for(int iLine=0; iLine<ncv_global; iLine++)
        {
          printf("\b\b\b\b\b\b\b\b\b%09d", iLine+1);
          rtnval = fscanf(fp, "%d", &gloId); locId = local_id[gloId];
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi_rho[locId]     = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi_rhou[locId][0] = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi_rhou[locId][1] = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi_rhou[locId][2] = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi_rhoE[locId]    = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi__kine[locId]   = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi__omega[locId]  = temp;
          rtnval = fscanf(fp, "%le", &temp); if(locId>=0) psi__gammai[locId] = temp;
        }
        printf("\n");
        fclose(fp);
      }
      SERIAL_END;
      
      updateCvDataG1G2(psi_rho, REPLACE_DATA);
      updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
      updateCvDataG1G2(psi_rhoE, REPLACE_DATA);
      updateCvDataG1G2(psi__kine, REPLACE_DATA);
      updateCvDataG1G2(psi__omega, REPLACE_DATA);
      updateCvDataG1G2(psi__gammai, REPLACE_DATA);
      PRINT("\n\n\n");
    }

    void connectPointers(const std::string &name, adouble **q, adouble (**grad_q)[3], adouble **q_diff)
    {
      int nScal = scalarTranspEqVector.size();
      for(int i=0; i<nScal; i++)
      {
        if(strcmp(scalarTranspEqVector_AD[i].name, name.c_str())==0)
        {
          (*q)      = scalarTranspEqVector_AD[i].phi;
          (*grad_q) = scalarTranspEqVector_AD[i].grad_phi;
          (*q_diff) = scalarTranspEqVector_AD[i].diff;
        }
      }
    }

    void calcRansTurbViscMuet_AD
    (
      adouble *rho,
      adouble (*rhou)[3]
    )
    {
      calcStrainRateAndDivergence_AD(); calcVorticity_AD();
    
      for(int ifa=nfa_b; ifa<nfa_b2gg; ifa++)
      if(ifa<nfa || (ifa>=nfa_b2 && ifa<nfa_b2gg))
      {
        double w0, w1; int icv0, icv1;
        getInterpolationFactors(ifa, icv0, icv1, w0, w1);

        adouble rho_fa  = w1*rho[icv0]   + w0*rho[icv1];
        adouble kine_fa = w1*kine[icv0]  + w0*kine[icv1];
        adouble om_fa   = w1*omega[icv0] + w0*omega[icv1];
        
        mut_fa[ifa] = rho_fa*kine_fa/om_fa;
      }

      for(list<FaZone>::iterator ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
      if(ZZ->getKind() == FA_ZONE_BOUNDARY)
      {
        if(zoneIsWall(ZZ->getName()))
        for(int index=0; index<ZZ->faVec.size(); index++)
        {
          int ifa = ZZ->faVec[index];
          mut_fa[ifa] = 0.0;        
        }
        if(!zoneIsWall(ZZ->getName()))
        for(int index=0; index<ZZ->faVec.size(); index++)
        {
          int ifa = ZZ->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          
          mut_fa[ifa] = rho[icv1]*kine[icv1]/omega[icv1];
        }
      }
    }

    void diffusivityHookScalarRansTurb_AD(const std::string &name)
    {
      ScalarTranspEq *eq; double sigma, sigmal; adouble *diff;
      
      if((name == "kine") || (name == "omega") || (name =="gammai"))
      {
        if(name=="kine")   { diff = kine_diff;   sigma=0.5, sigmal=1.0; }
        if(name=="omega")  { diff = omega_diff;  sigma=0.5, sigmal=1.0; }
        if(name=="gammai") { diff = gammai_diff; sigma=5.0; sigmal=5.0; }
        
        for (int ifa=nfa_b; ifa<nfa_b2gg; ifa++)
        if(ifa<nfa || (ifa>=nfa_b2 && ifa<nfa_b2gg))
        {
          double w0, w1; int icv0, icv1; getInterpolationFactors(ifa, icv0, icv1, w0, w1);

          adouble rho_fa  = (w1*rho_AD[icv0] + w0*rho_AD[icv1])/(w0+w1);
          adouble kine_fa = (w1*kine[icv0]   + w0*kine[icv1])/(w0+w1);
          adouble om_fa   = (w1*omega[icv0]  + w0*omega[icv1])/(w0+w1);

          diff[ifa] = mul_fa[ifa]/sigmal + sigma*rho_fa*kine_fa/om_fa;
        }

        for(list<FaZone>::iterator ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
        if(ZZ->getKind() == FA_ZONE_BOUNDARY)
        {
          for(int index=0; index<ZZ->faVec.size(); index++)
          {
            int ifa = ZZ->faVec[index];
            int icv0 = cvofa[ifa][0];
            
            diff[ifa] = mul_fa[ifa]/sigmal;
          }
          if(!zoneIsWall(ZZ->getName()))
          for(int index=0; index<ZZ->faVec.size(); index++)
          {
            int ifa = ZZ->faVec[index];
            int icv1 = cvofa[ifa][1];
            
            diff[ifa] = mul_fa[ifa]/sigmal + sigma*rho_AD[icv1]*kine[icv1]/omega[icv1];
          }
        }
      }
    }

    void boundaryHookScalarRansTurb_AD
    (
      adouble *phi,
      FaZone *zone,
      const std::string &name
    )
    {
      if(name=="omega")
      if(zone->getKind()==FA_ZONE_BOUNDARY)
      if(zoneIsWall(zone->getName()))
      {
        for(int index=0; index<zone->faVec.size(); index++)
        {
          int ifa = zone->faVec[index];
          int icv0 = cvofa[ifa][0];
          int icv1 = cvofa[ifa][1];
          adouble mu_cv = calcMuLam_AD(icv0);
          phi[icv1] = 6.0*mu_cv/(rho_AD[icv0]*0.075*pow(wallDist[icv0],2));
        }
      }
    }
    
    virtual void sourceHookRansTurbCoupled_AD
    (
      adouble **rhs,
      double ***A,
      int flagImplicit
    )
    {
      adouble src, Prod, nuT, nuL, d_src; int noc00, index;
    
      for (int icv=0; icv<ncv_g; icv++)
      {
        //PRINT("Checkpoint1 icv=%d\n", icv);
        noc00 = nbocv_i[icv];
        nuT = kine[icv]/omega[icv];
        muLam_AD[icv] = calcMuLam_AD(icv);
        nuL = muLam_AD[icv]/rho_AD[icv];
        Prod = fmax(nuT*vortMag[icv]*vortMag[icv] - (2./3.)*kine[icv]*diverg[icv], 0.0);
        
        // kine :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //PRINT("Checkpoint2 icv=%d\n", icv);
        index = getScalarTransportIndex("kine");
        src   = fmin(2.0, gammai[icv]) * Prod - 0.09 * kine[icv] * omega[icv];
        d_src = -0.09 * omega[icv];
    
        rhs[icv][5+index] += rho_AD[icv] * src * cv_volume[icv];
        if(flagImplicit && icv<ncv) A[noc00][5+index][5+index] -= d_src.value() * cv_volume[icv];
    
        // amega ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //PRINT("Checkpoint3 icv=%d\n", icv);
        index = getScalarTransportIndex("omega");
        src   = (5./9.) * omega[icv] / kine[icv] * Prod - 0.075 * omega[icv] * omega[icv];
        d_src = -0.15 * omega[icv];
    
        rhs[icv][5+index] += rho_AD[icv] * src * cv_volume[icv];
        if(flagImplicit && icv<ncv) A[noc00][5+index][5+index] -= d_src.value() * cv_volume[icv];

        // gammai :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //PRINT("Checkpoint4 icv=%d\n", icv);
        gamma_max.calculate_AD(icv);
        index = getScalarTransportIndex("gammai");
        src   = (gamma_max.beta_AD - gammai[icv]) * vortMag[icv] * sqrt(gammai[icv]);
        d_src = -1.5 * vortMag[icv] * sqrt(gammai[icv]);
    
        rhs[icv][5+index] += rho_AD[icv] * src * cv_volume[icv];
        if(flagImplicit && icv<ncv) A[noc00][5+index][5+index] -= d_src.value() * cv_volume[icv];
      }
    }

    adouble sigmoid_AD(adouble x, adouble x0, adouble s)
    {
      return 1.0 / (1.0 + exp((x0-x)/s));
    }
};

#endif
