#include <unistd.h>
#include "macros.h"
#include "JoeWithModels.h"
#include "RansTurb.h"
#include "SurfaceObservables/SurfaceObjective.h"

#ifdef ADJOINT_MODE
  #include "JoeWithModelsAD.h"
  #include "RansTurbAD.h"
#endif

#ifndef ADJOINT_MODE
class MyJoe : virtual public JoeWithModels, virtual public RansTurb
#else
class MyJoe : virtual public JoeWithModels_AD, virtual public RansTurb_AD
#endif
{
  public:
    
    std::vector<std::string> markerNames;
    #ifdef ADJOINT_MODE
    SurfaceObjective<adouble> surfObj;
    #else
    SurfaceObjective<double> surfObj;
    #endif
  
  public:
    
    #include "SurfaceObservables/extractFacesForData.h"
    #include "SurfaceObservables/extractMarkerData.h"
    
    #ifndef ADJOINT_MODE
    MyJoe(char *name) : JoeWithModels(name), UgpWithCvCompFlow(name)
    #else
    MyJoe(char *name) : JoeWithModels(name), JoeWithModels_AD(name), UgpWithCvCompFlow(name), UgpWithCvCompFlow_AD(name)
    #endif
    {
      PRINT("MyJoe::MyJoe()\n");
      Param *param = getParam("MARKER_NAMES");
      for(int i=1; i<param->getSize(); i++)
        markerNames.push_back(param->getString(i));
    }

    ~MyJoe() { }
    
    void initialHook()
    {
      PRINT("MyJoe::initialHook()\n");
      JoeWithModels::initialHook();
      #ifndef ADJOINT_MODE
      if(checkParam("CALCULATE_OBJECTIVE"))
      {
        surfObj.readFile(getStringParam("REF_FILENAME", "_.dat"));
        surfObj.indexForFaces.resize(nfa_b, -1);
        setSurfaceObjectiveFaces();
        surfObj.write("verification.dat");
      }

      Param *param = getParam("U_INITIAL");
      double U2 = param->getDouble(1)*param->getDouble(1)
                + param->getDouble(2)*param->getDouble(2)
                + param->getDouble(3)*param->getDouble(3);

			double ReThetaT_val = sqrt(U2 * (7.0 * getDoubleParam("RHO_REF"))/(9.0 * getDoubleParam("omega_INITIAL") * getDoubleParam("MU_REF")));
			double mu_T_FS_val = getDoubleParam("RHO_REF") * getDoubleParam("kine_INITIAL") / getDoubleParam("omega_INITIAL");
      
			for(int icv=0; icv<ncv_gg; icv++)
      {
        ReThetaT[icv] = ReThetaT_val;
        mu_T_FS[icv] = mu_T_FS_val;
      }

      PRINT("MyJoe::initialHook() completed\n");
      #endif
    }
    
    void initialHookScalarRansTurbModel(void)
    {
      RansTurb::initialHookScalarRansTurbModel();
    }
    
    void boundaryHook(double *temp, double (*vel)[3], double *press, FaZone *zone)
    {
      boundaryFn<double>(temp, vel, press, zone);
    }
    
    template<typename DBLE> void boundaryFn(DBLE *temp, DBLE (*vel)[3], DBLE *press, FaZone *zone) { }

    #ifndef ADJOINT_MODE
    void temporalHook()
    {
      if(checkParam("USE_FS"))
      {
				bool checkStep = (step%100==0);
        if(checkStep) { PRINT("Setting FS data for step\n"); }
        setFS(&wLink, wallDist, ReThetaT_local, ReThetaT, "ReThetaT.dat", checkStep);
        setFS(&wLink, wallDist, muT, mu_T_FS, "mu_T_FS.dat", checkStep);
      }

      if(step%5==0)
      {
        writeFeatures();
      }

      if(step%100==0)
      {
        extractMarkerData<double>(rho, rhou, rhoE);
        if(checkParam("CALCULATE_OBJECTIVE"))
        {
          FILE *fp;
          surfObj.write("markerObjective.out");
          surfObj.evaluate();
          PRINT("------------------------------------------------------------------------\n");
          for(int i=0; i<surfObj.nTargets; i++)
          {
            PRINT(" OBJECTIVE      %d: %le\n", i, surfObj.total_obj[i] / surfObj.nFacesTotal);
          }
          PRINT("------------------------------------------------------------------------\n");

          if(mpi_rank==0)
          {
            fp = fopen("obj.dat", "w");
            for(int i=0; i<surfObj.nTargets; i++)
              fprintf(fp, "%+.10le\n", surfObj.total_obj[i] / surfObj.nFacesTotal);
            fclose(fp);
          }
        }
        SERIAL_BEG;
        FILE *fp;
        if(mpi_rank==0) fp = fopen("features.dat", "w");
        if(mpi_rank!=0) fp = fopen("features.dat", "a");
        for(int icv=0; icv<ncv; icv++)
        {
          gamma_max.evalFtrs(icv);
          fprintf(fp, "%9d", (int)(global_id[icv]));
          for(int i=0; i<gamma_max.ftrs.size(); i++)
          {
            fprintf(fp, "  %+.12le", gamma_max.ftrs[i]);
          }
          fprintf(fp, "\n");
        }
        fclose(fp);
        SERIAL_END;
      }
    }
    #endif

    void finalHook() { if(mpi_rank==0) system("touch completed"); }

    #ifdef ADJOINT_MODE
    void initialHook_AD()
    {
      PRINT("MyJoe::initialHook_AD()\n");
      surfObj.readFile(getStringParam("REF_FILENAME", "_.dat"));
      surfObj.indexForFaces.resize(nfa_b, -1);
      setSurfaceObjectiveFaces();
      surfObj.write("verification.dat");
    }

    void boundaryHook_AD(adouble *temp, adouble (*vel)[3], adouble *press, FaZone *zone)
    {
      boundaryFn<adouble>(temp, vel, press, zone);
    }

    void calcFunctional_AD(adouble *rho, adouble (*rhou)[3], adouble *rhoE)
    {
      extractMarkerData<adouble>(rho, rhou, rhoE);
      surfObj.evaluate();
      functional = surfObj.local_obj;
    }

    void temporalHook_AD(int step_val) { if(step_val%100==0 && step_val>0) { copy_turb_adjoint(); writeDerivsFIML(); writeAdjointsToFile(); } }
    
    void finalHook_AD() { if(mpi_rank==0) system("touch completed"); }
    #endif
};


int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv); MPI_Comm mpi_comm_val = MPI_COMM_WORLD;
	
	initMpiStuff(mpi_comm_val); int run = 1; char inputfilename[50]; sprintf(inputfilename, "Joe.in");
	
  for(int iArg=1; iArg<argc; iArg+=2)
  {
		if(!strcmp((char*)("--file"), argv[iArg])) strcpy(inputfilename, argv[iArg+1]);
		if(!strcmp((char*)("--run"), argv[iArg]))  run = atoi(argv[iArg+1]);
	}

	if (mpi_rank==0) printf("Specified Input name = %s\n", inputfilename);
	if (mpi_rank==0) printf("Specified Run        = %d\n", run);
	
  try
  {
		MyJoe *joe = new MyJoe(inputfilename);
		double wtime, wtime0; MPI_Barrier(mpi_comm); if (mpi_rank==0) wtime = MPI_Wtime();
		
    switch(run)
    {
			case 0    : joe->setGlobalId();       break;
      #ifndef ADJOINT_MODE
			case 1    : joe->run(1);              break;
			#else
			case 1    : joe->runAdjoint(1);       break;
			#endif
		}
		
    MPI_Barrier(mpi_comm); wtime0 = wtime; wtime = MPI_Wtime();
		if (mpi_rank==0) printf(" > Total runtime [s]: %le\n", wtime - wtime0);
		delete joe;
	}
	
  catch(int e) { cerr << "Exception: " << e << endl; return -1; }
	
  catch(...) { cerr << "Unhandled Exception\n" << endl; return -1; }
	
  MPI_Barrier(mpi_comm_val); MPI_Finalize(); system("rm -f ADOLC-*");
	return 0;
}
