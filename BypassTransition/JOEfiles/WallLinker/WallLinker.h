struct WallLinker
{
  int    nWallFaces;
  int    nMyWallFaces;
  int    *send_counts;
  int    *send_displs;
  double *send_buffer;
  int    *recv_counts;
  int    *recv_displs;
  double *recv_buffer;
  double *wallFace_i;
  double *globalVar;
  double *globalCtr;
  double *localVar;
  double *localCtr;
};
