void writeParams()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0)
  {
    FILE *fp = fopen((name+"_param_values.dat").c_str(), "w");
    
    fprintf(fp, "%d", nodes.size());
    
    for(int i=0; i<nodes.size(); i++)
    {
      fprintf(fp, " %d %d", nodes[i].size(), act[i]);
    }

    for(int i=1; i<nodes.size(); i++)
    {
      for(int j=0; j<nodes[i].size; j++)
      {
        fprintf(fp, "\n%+.12le", biases[i][j]);

        for(int k=0; k<nodes[i-1].size(); k++)
        {
          fprintf(fp, "\n%+.12le", weights[i][j][k]);
        }
      }
    }

    fclose(fp);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
