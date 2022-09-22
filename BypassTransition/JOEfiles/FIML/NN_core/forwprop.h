void forwprop(void)
{
  for(int i=1; i<nodes.size(); i++)
  {
    for(int j=0; j<nodes[i].size(); j++)
    {
      nodes[i][j] = biases[i][j];

      for(int k=0; k<nodes[i-1].size(); k++)
      {
        nodes_deact[i][j] += weights[i][j][k] * nodes[i-1][k];
      }

      nodes[i][j] = activate(nodes_deact[i][j], act[i]);
    }
  }
}

void forwprop_AD(void)
{
  for(int i=1; i<nodes_AD.size(); i++)
  {
    for(int j=0; j<nodes_AD[i].size(); j++)
    {
      nodes_AD[i][j] = biases[i][j];

      for(int k=0; k<nodes_AD[i-1].size(); k++)
      {
        nodes_AD[i][j] += weights[i][j][k] * nodes_AD[i-1][k];
      }

      nodes_AD[i][j] = activate(nodes_AD[i][j], act[i]);
    }
  }
}
