template<typename DBLE>
void extractMarkerData(DBLE *rho, DBLE (*rhou)[3], DBLE *rhoE)
{
  #ifndef ADJOINT_MODE
  FILE *fp;
  #endif
  //-----------------
  // Variables
  //-----------------
  int icv0 = 0, icv1 = 0;
  double area = 0.0, walld = 0.0;
  double nVec[3] = {0.0, 0.0, 0.0}, sVec[3] = {0.0, 0.0, 0.0};
  DBLE u[3], un, u_face[3];
  DBLE press_face, Cp;
  DBLE temp_cell, temp_face, q_face, lambda, Ch;
  DBLE muTotal, tau_wall[3], Cf;
  //---------------------------------
  // Calculate dynamic pressure
  //---------------------------------
  Param *param   = getParam("U_INITIAL");
  double rho_ref = getDoubleParam("RHO_REF", "1.225");
  double p_ref   = getDoubleParam("P_REF",   "101325.0");
  double vel_ref = 0.0;
  for(int i=1; i<=3; i++)
  {
    vel_ref += pow(param->getDouble(i),2);
  }
  vel_ref = sqrt(vel_ref);
  double Q_ref   = 0.5 * rho_ref * pow(vel_ref, 2);
  //----------------------------------
  // Loop over marker faces
  //----------------------------------
  for(auto &markerName : markerNames)
  {
    for(auto Z=faZoneList.begin();Z!=faZoneList.end();Z++)
    {
      if(std::string(Z->getName())==markerName)
      {
        #ifndef ADJOINT_MODE
    		fp = fopen((markerName+"___"+std::to_string(mpi_rank)+".dat").c_str(), "w");
        #endif
    	for(int ifa=Z->ifa_f;ifa<=Z->ifa_l;ifa++)
        {
          //---------------------------------------------------------------
          // Calculate area, face normal and wall distance to nearest cell
          //---------------------------------------------------------------
    		  icv0=cvofa[ifa][0]; icv1=cvofa[ifa][1];
          area=normVec3d(nVec,fa_normal[ifa]);
          vecMinVec3d(sVec,x_fa[ifa],x_cv[icv0]);
          walld=fabs(vecDotVec3d(sVec,nVec));
          //---------------------------------------------------------------
          // Calculate normal and tangential velocities
          //---------------------------------------------------------------
          un = 0.0;
          for(int i=0; i<3; i++)
          {
            u[i]=rhou[icv0][i]/rho[icv0];
          }
          for(int i=0; i<3; i++)
          {
            un += u[i]*nVec[i];
          }
    		  for(int i=0; i<3; i++)
          {
            u_face[i]=u[i]-un*nVec[i];
          }
    		  //---------------------------------------------------------------
          // Calculate pressure and Cp
          //---------------------------------------------------------------
          press_face = (GAMMA-1.0) * 
                       (rhoE[icv0]-0.5*(rho[icv0]*u_face[0]*u_face[0] +
                                        rho[icv0]*u_face[1]*u_face[1] +
                                        rho[icv0]*u_face[2]*u_face[2]));
          Cp = (press_face - p_ref) / Q_ref;
          //---------------------------------------------------------------
          // Calculate temperature and heat transfer coefficient
          //---------------------------------------------------------------
          temp_face = press_face / (R_gas * rho[icv0]);
          temp_cell = (GAMMA-1.0) / R_gas / rho[icv0] *
                      (rhoE[icv0]-0.5*(rhou[icv0][0]*u[0] +
                                       rhou[icv0][1]*u[1] +
                                       rhou[icv0][2]*u[2]));
          lambda    = GAMMA*R_gas/(GAMMA-1.0)*muTotal/Pr;
          q_face    = lambda*(temp[icv1]-temp_cell)/walld;
          if(temp_face!=temp[icv1])
          {
            Ch = q_face/(temp[icv1]-temp_face);
          }
          else
          {
            Ch = 0.0;
          }
          //---------------------------------------------------------------
          // Calculate viscosity and skin friction
          //---------------------------------------------------------------
          if (viscMode == "SUTHERLAND")
          {
            muTotal = mu_ref * pow(temp_face/SL_Tref,1.5) * (SL_Tref + SL_Sref) / (temp_face + SL_Sref);
          }
    		  if (viscMode == "POWERLAW")
          {
            muTotal = mu_ref * pow(temp_face/T_ref, mu_power_law);
          }
          for (int i=0; i<3; i++)
          {
            tau_wall[i] = muTotal*u_face[i]/walld;
          }
    		  Cf = sqrt(pow(tau_wall[0],2) + pow(tau_wall[1],2) + pow(tau_wall[2],2)) / Q_ref;
          #ifndef ADJOINT_MODE
          //--------------------------------------------------------------
          // Write data to file
          //--------------------------------------------------------------
          fprintf(fp, "%+.10le\t", x_fa[ifa][0]);
          fprintf(fp, "%+.10le\t", x_fa[ifa][1]);
          fprintf(fp, "%+.10le\t", x_fa[ifa][2]);
          fprintf(fp, "%+.10le\t", area);
          fprintf(fp, "%+.10le\t", nVec[0]);
          fprintf(fp, "%+.10le\t", nVec[1]);
          fprintf(fp, "%+.10le\t", nVec[2]);
          fprintf(fp, "%+.10le\t", u_face[0]);
          fprintf(fp, "%+.10le\t", press_face);
          fprintf(fp, "%+.10le\t", tau_wall[0]);
          fprintf(fp, "%+.10le\t", tau_wall[1]);
          fprintf(fp, "%+.10le\t", tau_wall[2]);
          fprintf(fp, "%+.10le\t", temp_face);
          fprintf(fp, "%+.10le\t", q_face);
          fprintf(fp, "%+.10le\t", Cp);
          fprintf(fp, "%+.10le\t", Cf);
          fprintf(fp, "%+.10le\n", Ch);
          #endif
          //--------------------------------------------------------------
          // Calculate observables
          //--------------------------------------------------------------
          if(checkParam("CALCULATE_OBJECTIVE"))
          {
            if(surfObj.indexForFaces[ifa]>=0)
            {
              int iLoc = surfObj.indexForFaces[ifa];
              //------------------------------------------
              // __CHANGE__ here to change observables
              //------------------------------------------

              if(getStringParam("CASE_NAME")=="NACA")
              {
                if(x_fa[ifa][0]>0.08 && x_fa[ifa][0]<0.116)
                {
                  surfObj.values[iLoc][0] = sqrt(tau_wall[0]*tau_wall[0] + tau_wall[1]*tau_wall[1]);
                  
                  if(tau_wall[0]<0.0)
                  {
                    surfObj.values[iLoc][0] *= -1;
                  }
                }
                else
                {
                  surfObj.values[iLoc][0] = surfObj.targets_interp[iLoc][0];
    		        }
              }

              if(getStringParam("CASE_NAME")=="T3")
              {
                if(x_fa[ifa][0]>0.025 && x_fa[ifa][0]<1.7)
                {
                  surfObj.values[iLoc][0] = sqrt(tau_wall[0]*tau_wall[0] + tau_wall[1]*tau_wall[1]) / Q_ref;
                  
                  if(tau_wall[0]<0.0)
                  {
                    surfObj.values[iLoc][0] *= -1;
                  }
                }
                else
                {
                  surfObj.values[iLoc][0] = surfObj.targets_interp[iLoc][0];
    		        }
              }
            }
          }
        }
        #ifndef ADJOINT_MODE
    		fclose(fp);
        #endif
    	}
    }
  }
  #ifndef ADJOINT_MODE
  for(auto markerName:markerNames)
  {
  	MPI_Barrier(mpi_comm);
  	if(mpi_rank==0)
    {
      system(("rm -f "+markerName+".dat; cat "+markerName+"___*.dat >> "+markerName+".dat; rm "+markerName+"___*.dat;").c_str());
  	}
    MPI_Barrier(mpi_comm);
  }
  #endif
}
