void backprop(void)
{
  for(int i=nodes.size()-1; i>0; i--)
  {
    for(int k=0; k<nodes[i-1].size(); k++) d_nodes[i-1][k] = 0.0;
    
    for(int j=0; j<nodes[i].size(); j++)
    {
      d_nodes[i][j]  = d_nodes[i][j] * d_activate(nodes_deact[i][j], act[i]);
      d_biases[i][j] = d_nodes[i][j];

      for(int k; k<nodes[i-1].size(); k++)
      {
        d_nodes[i-1][k]   += d_nodes[i][j] * weights[i][j][k];
        d_weights[i][j][k] = d_nodes[i][j] * nodes[i-1][k];
      }
    }
  }
}
