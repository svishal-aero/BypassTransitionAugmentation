void initWallLinker(WallLinker *wLink)
{
  if(mpi_rank==0) printf("Initializing wall linker\n");

  // Initialize the number of wall boundary faces to zero
  wLink->nMyWallFaces = 0;
  
  // Add number of wall boundary faces from all wall boundary zones on this processor
  for(auto ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
  if(ZZ->getKind()==FA_ZONE_BOUNDARY)
  {
    Param *param; if(getParam(param, ZZ->getName())) if(param->getString()=="WALL")
    for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++) wLink->nMyWallFaces++;
  }

  // Sum up number of wall boundary faces from all processors
  MPI_Allreduce(&(wLink->nMyWallFaces), &(wLink->nWallFaces), 1, MPI_INT, MPI_SUM, mpi_comm);

  // Set up the sending part of the communication
  wLink->send_counts = new    int[mpi_size];
  wLink->send_displs = new    int[mpi_size];
  wLink->send_buffer = new double[wLink->nMyWallFaces];

  // Set up the receiving part of the communication
  wLink->recv_counts = new    int[mpi_size];
  wLink->recv_displs = new    int[mpi_size];
  wLink->recv_buffer = new double[wLink->nWallFaces];

  // Set up the array containing global indices of wall boundary faces closest to all cv's
  wLink->wallFace_i  = NULL;
  registerScalar(wLink->wallFace_i, "wallFace_i", CV_DATA);

  // Set up counters (to check repetition) and buffers to temporarily store variables
  wLink->globalVar   = new double[wLink->nWallFaces];
  wLink->globalCtr   = new double[wLink->nWallFaces];
  wLink->localVar    = new double[wLink->nWallFaces];
  wLink->localCtr    = new double[wLink->nWallFaces];

  // Set the send_counts as nMyWallFaces and send_displs to zero for all processors
  for(int i=0; i<mpi_size; i++){ wLink->send_counts[i] = wLink->nMyWallFaces; wLink->send_displs[i] = 0; }

  // Gather all the send_counts into recv_counts (must be the same mpi_size-sized array on all processors)
  MPI_Allgather(&(wLink->nMyWallFaces), 1, MPI_INT, wLink->recv_counts, 1, MPI_INT, mpi_comm);

  // Calculate recv_displs from recv_counts
  wLink->recv_displs[0] = 0; for(int i=1; i<mpi_size; i++){ wLink->recv_displs[i] = wLink->recv_displs[i-1] + wLink->recv_counts[i-1]; }

  if(mpi_rank==0) printf("Finalizing wall linker\n");
}

void deleteWallLinker(WallLinker *wLink)
{
  delete[] wLink->send_counts;
  delete[] wLink->send_displs;
  delete[] wLink->send_buffer;
  delete[] wLink->recv_counts;
  delete[] wLink->recv_displs;
  delete[] wLink->recv_buffer;
  delete[] wLink->globalVar;
  delete[] wLink->globalCtr;
  delete[] wLink->localVar;
  delete[] wLink->localCtr;
}

void calcWallDistanceModified(WallLinker *wLink, double *wd)
{
  PRINT("calcWallDistanceModified: Checkpoint 1\n");

  // Create send and receive buffers
  double *myFaceVars = new double[wLink->nMyWallFaces*6];
  double *faceVars = new double[wLink->nWallFaces*6];
  
  // Create a counter for send buffer entry
  int count = 0;

  // Fill the send buffer in with coordinates and normals of all wall boundary faces on this processor
  for(auto ZZ=faZoneList.begin(); ZZ!=faZoneList.end(); ZZ++)
  if(ZZ->getKind()==FA_ZONE_BOUNDARY)
  {
    Param *param; if(getParam(param, ZZ->getName())) if(param->getString()=="WALL")
    for(int ifa=ZZ->ifa_f; ifa<=ZZ->ifa_l; ifa++)
    {
      myFaceVars[wLink->nMyWallFaces*0+count] = x_fa[ifa][0];
      myFaceVars[wLink->nMyWallFaces*1+count] = x_fa[ifa][1];
      myFaceVars[wLink->nMyWallFaces*2+count] = x_fa[ifa][2];
      myFaceVars[wLink->nMyWallFaces*3+count] = fa_normal[ifa][0];
      myFaceVars[wLink->nMyWallFaces*4+count] = fa_normal[ifa][1];
      myFaceVars[wLink->nMyWallFaces*5+count] = fa_normal[ifa][2];
      count++;
    }
  }

  // Check that the counting is correct
  assert(count==wLink->nMyWallFaces);

  // Broadcast all quantities
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*0]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*0]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*1]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*1]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*2]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*2]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*3]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*3]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*4]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*4]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  MPI_Alltoallv(&(myFaceVars[wLink->nMyWallFaces*5]), wLink->send_counts, wLink->send_displs, MPI_DOUBLE,
                    &(faceVars[wLink->nWallFaces*5]), wLink->recv_counts, wLink->recv_displs, MPI_DOUBLE, mpi_comm);
  
  delete[] myFaceVars;
  
  // Find minimum wall distance =========================================================

  PRINT("calcWallDistanceModified: Starting loop over cells (ncv=%d)\n", ncv);
  for (int icv=0; icv<ncv; icv++)
  {
    double minDist = 1.0e20; int ifaSave;

    // Find the wall boundary face (global) with the minimum distance from the cv
    for (int ifa=0; ifa<wLink->nWallFaces; ifa++)
    {
      double dist = sqrt(pow(x_cv[icv][0]-faceVars[wLink->nWallFaces*0+ifa], 2.0)
                       + pow(x_cv[icv][1]-faceVars[wLink->nWallFaces*1+ifa], 2.0)
                       + pow(x_cv[icv][2]-faceVars[wLink->nWallFaces*2+ifa], 2.0));
      
      if (dist < minDist) { ifaSave = ifa; minDist = dist; }
    }

    // Evaluate the face area (nmag) and displacement vector (s_half)
    double fa_x[3] = {faceVars[wLink->nWallFaces*0+ifaSave], faceVars[wLink->nWallFaces*1+ifaSave], faceVars[wLink->nWallFaces*2+ifaSave]};
    double fa_n[3] = {faceVars[wLink->nWallFaces*3+ifaSave], faceVars[wLink->nWallFaces*4+ifaSave], faceVars[wLink->nWallFaces*5+ifaSave]};
    double n[3], s_half[3];
    double nmag = normVec3d(n, fa_n);
    vecMinVec3d(s_half, fa_x, x_cv[icv]);
    
    // Find normal and tangential distance^2 of the cv center from face center
    double sMag2 = vecDotVec3d(s_half, s_half);
    double dNorm2 = vecDotVec3d(s_half, n)*vecDotVec3d(s_half, n);
    double dTang2 = sMag2-dNorm2;

    // Find the effective radius^2 of the face
    double specRadius2 = 0.5*nmag;
    
    // If the tangential distance^2 is within the effective radius^2, then set the distance as normal distance
    // Otherwise, set it to the minimum distance calculated in the loop above
    if(dTang2<=specRadius2)
    {
      wd[icv] = sqrt(dNorm2);
      wLink->wallFace_i[icv] = ifaSave;
    }
    else
    {
      wd[icv] = minDist;
      wLink->wallFace_i[icv] = -1;
    }
  }
  
  PRINT("calcWallDistanceModified: Checkpoint 3\n");
  MPI_Barrier(mpi_comm);
  
  delete[] faceVars;
}

void setFS(WallLinker *wLink, double *wd, double *q, double *qFS, std::string name, bool write)
{
  double distFS = getDoubleParam("DIST_FS");
  double tolFS = getDoubleParam("TOL_FS");

  // Zero out all the counters and buffers
  for(int ifa=0; ifa<wLink->nWallFaces; ifa++) { wLink->localCtr[ifa] = 0.0; wLink->localVar[ifa] = 0.0; } 

  // Check which cv's satisfy the freestream bounds, and if they do, assign the freestream variables
  for(int icv=0; icv<ncv; icv++)
  {
    int face_id = (int)(wLink->wallFace_i[icv]+0.5);
    if(fabs(wd[icv]-distFS)<tolFS && face_id>=0)
    {
      wLink->localCtr[face_id] += 1.0;
      wLink->localVar[face_id] += q[icv];
    }
  }

  // Sum up the contributions from all processors on all processors
  MPI_Allreduce(wLink->localCtr, wLink->globalCtr, wLink->nWallFaces, MPI_DOUBLE, MPI_SUM, mpi_comm);
  MPI_Allreduce(wLink->localVar, wLink->globalVar, wLink->nWallFaces, MPI_DOUBLE, MPI_SUM, mpi_comm);
 
  FILE *fp; if(mpi_rank==0 && write) fp = fopen(name.c_str(), "w");

  // Normalize with counters to avoid redundant summing above
  for(int ifa=0; ifa<wLink->nWallFaces; ifa++)
  {
    if(wLink->globalCtr[ifa]==0.0) wLink->globalCtr[ifa]=1.0;
    wLink->globalVar[ifa] /= wLink->globalCtr[ifa];
    if(mpi_rank==0 && write) fprintf(fp, "%6d\t%+.10le\n", ifa, wLink->globalVar[ifa]);
  }

  if(mpi_rank==0 && write) fclose(fp);

  // Finally, set all the cv's to freestream values (within bounds)
  for(int icv=0; icv<ncv; icv++)
  {
    int face_id = (int)(wLink->wallFace_i[icv]+0.5);
    qFS[icv] = wLink->globalVar[face_id];
  }

  updateCvDataG1G2(qFS, REPLACE_DATA);
}
