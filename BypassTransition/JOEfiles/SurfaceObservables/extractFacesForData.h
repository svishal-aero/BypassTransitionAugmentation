virtual void setSurfaceObjectiveFaces()
{
  FILE *fp;
  //--------------------------------------------
  // Initialize variables
  //--------------------------------------------
  double area, nVec[3];
  double x1[2], x2[2];  // Distance from successive datapoints
  double nDist1, nDist2;
  double fac1, fac2;
  double tol = getDoubleParam("MARKER_DATA_TOL");
  //--------------------------------------------
  // Loop over faces
  //--------------------------------------------
  for(auto &markerName : markerNames)
  {
    for(auto ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
    {
      if(std::string(ZZ->getName())==markerName)
      {
    		fp = fopen((markerName+"___"+std::to_string(mpi_rank)+".dat").c_str(), "w");
        PRINT("Assigning faces for marker \"%s\"...\n", markerName.c_str());
        for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
        {
          fprintf(fp, "%+.10le\t", x_fa[ifa][0]);
          fprintf(fp, "%+.10le\t", x_fa[ifa][1]);
          fprintf(fp, "%+.10le\n", x_fa[ifa][2]);
          area=normVec3d(nVec,fa_normal[ifa]);
          //------------------------------------------------------
          // Calculate x2, nDist2 and fac2 for first datapoint
          //------------------------------------------------------
          for(int i=0; i<2; i++)
          {
            x2[i] = surfObj.x_targets[0][i] - x_fa[ifa][i];
          }
          nDist2 = x2[0]*nVec[0] + x2[1]*nVec[1];
          fac2 = 0.0;
          for(int i=0; i<2; i++)
          {
            x2[i] -= nDist2*nVec[i];
            fac2 += x2[i]*x2[i];
          }
          fac2 = sqrt(fac2);
          //------------------------------------------------------
          // Loop over datapoints
          //------------------------------------------------------
          for(int iLoc=1; iLoc<surfObj.nTargetLocations; iLoc++)
          {
            for(int i=0; i<2; i++)
            {
              x1[i] = x2[i];
            }
            nDist1 = nDist2;
            fac1 = fac2;
            //----------------------------------------------------
            // Calculate x2, nDist2 and fac2 for current datapoint
            //----------------------------------------------------
            for(int i=0; i<2; i++)
            {
              x2[i] = surfObj.x_targets[iLoc][i] - x_fa[ifa][i];
            }
            nDist2 = x2[0]*nVec[0] + x2[1]*nVec[1];
            fac2 = 0.0;
            for(int i=0; i<2; i++)
            {
              x2[i] -= nDist2*nVec[i];
              fac2 += x2[i]*x2[i];
            }
            fac2 = sqrt(fac2);
            //----------------------------------------------------
            // Check if the datapoint iLoc-1 is exactly
            // on the face center
            //----------------------------------------------------
            if(fac1==0)
            {
              surfObj.x_faces.push_back(std::array<double,3>{x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]});
              surfObj.indexForFaces[ifa] = surfObj.nFaces;
              surfObj.targets_interp.push_back(std::vector<double>(surfObj.nTargets));
              for(int iTgt=0; iTgt<surfObj.nTargets; iTgt++)
              {
                surfObj.targets_interp[surfObj.nFaces][iTgt] = surfObj.targets[iLoc-1][iTgt];
              }
              #ifndef ADJOINT_MODE
              surfObj.values.push_back(std::vector<double>(surfObj.nTargets));
              #else
              surfObj.values.push_back(std::vector<adouble>(surfObj.nTargets));
              #endif
              surfObj.nFaces++;
              break;
            }
            //----------------------------------------------------
            // Check if the datapoint iLoc is exactly
            // on the face center
            //----------------------------------------------------
            else if(fac2==0)
            {
              surfObj.x_faces.push_back(std::array<double,3>{x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]});
              surfObj.indexForFaces[ifa] = surfObj.nFaces;
              surfObj.targets_interp.push_back(std::vector<double>(surfObj.nTargets));
              for(int iTgt=0; iTgt<surfObj.nTargets; iTgt++)
              {
                surfObj.targets_interp[surfObj.nFaces][iTgt] = surfObj.targets[iLoc][iTgt];
              }
              #ifndef ADJOINT_MODE
              surfObj.values.push_back(std::vector<double>(surfObj.nTargets));
              #else
              surfObj.values.push_back(std::vector<adouble>(surfObj.nTargets));
              #endif
              surfObj.nFaces++;
              break;
            }
            //----------------------------------------------------
            // Check if the face center is between datapoints
            //----------------------------------------------------
            else if(((x1[0]*x2[0] + x1[1]*x2[1])<=0.0) && (nDist1<tol) && (nDist2<tol))
            {
              surfObj.x_faces.push_back(std::array<double,3>{x_fa[ifa][0], x_fa[ifa][1], x_fa[ifa][2]});
              surfObj.indexForFaces[ifa] = surfObj.nFaces;
              surfObj.targets_interp.push_back(std::vector<double>(surfObj.nTargets));
              for(int iTgt=0; iTgt<surfObj.nTargets; iTgt++)
              {
                surfObj.targets_interp[surfObj.nFaces][iTgt] = (surfObj.targets[iLoc-1][iTgt]*fac2 + surfObj.targets[iLoc][iTgt]*fac1)/(fac1+fac2);
              }
              #ifndef ADJOINT_MODE
              surfObj.values.push_back(std::vector<double>(surfObj.nTargets));
              #else
              surfObj.values.push_back(std::vector<adouble>(surfObj.nTargets));
              #endif
              surfObj.nFaces++;
              break;
            }
          }
        }
        fclose(fp);
      }
    }
  }
  PRINT("\n");
  MPI_Allreduce(&surfObj.nFaces, &surfObj.nFacesTotal, 1, MPI_INT, MPI_SUM, mpi_comm);
  PRINT("Number of faces recorded for markers: %d\n\n", surfObj.nFacesTotal);
  for(auto markerName:markerNames)
  {
  	MPI_Barrier(mpi_comm);
  	if(mpi_rank==0)
    {
      system(("rm -f "+markerName+".dat; cat "+markerName+"___*.dat >> "+markerName+".dat; rm "+markerName+"___*.dat;").c_str());
  	}
    MPI_Barrier(mpi_comm);
  }
}
