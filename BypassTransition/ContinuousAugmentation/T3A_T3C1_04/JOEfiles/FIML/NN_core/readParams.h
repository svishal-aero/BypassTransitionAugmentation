void readParams()
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int temp;

  if(rank>0) MPI_Recv(NULL, NULL, 0, MPI_CHAR, rank-1, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
  FILE *fp = fopen((name+"_param_values.dat").c_str(), "r");
    
  fscanf(fp, "%d", &temp); assert(temp==nodes.size());
  
  for(int i=0; i<nodes.size(); i++)
  {
    fscanf(fp, "%d", &temp); assert(temp==nodes[i].size());
    fscanf(fp, "%d", &temp); assert(temp==act[i]);
  }

  for(int i=1; i<nodes.size(); i++)
  {
    for(int j=0; j<nodes[i].size; j++)
    {
      fscanf(fp, "%le", &(biases[i][j]));

      for(int k=0; k<nodes[i-1].size(); k++)
      {
        fscanf(fp, "%le", &(weights[i][j][k]));
      }
    }
  }

  fclose(fp);

  if(rank<size-1) MPI_Send(NULL, NULL, 0, MPI_CHAR, rank+1, 123, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);
}
