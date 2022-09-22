void init(std::string name_val, std::vector<int> nNodes, std::vector<int> actFun, RansTurb *rans_val, RansTurb_AD *rans_AD_val)
{
  name = name_val;
  act = actFun;
  
  features.resize(nNodes[0], 0.0);
  beta.resize(nNodes[nNodes.size()-1], 0.0);
  nodes.resize(nNodes.size()); nodes_deact.resize(nNodes.size()); d_nodes.resize(nNodes.size());
  biases.resize(nNodes.size()); d_biases.resize(nNodes.size());
  weights.resize(nNodes.size()); d_weights.resize(nNodes.size());
  
  for(int i=0; i<nNodes.size(); i++)
  {
    nodes[i].resize(nNodes[i], 0.0); nodes_deact.resize(nNodes[i], 0.0); d_nodes[i].resize(nNodes[i], 0.0);
    biases[i].resize(nNodes[i], 0.0); d_biases[i].resize(nNodes[i], 0.0);
    weights[i].resize(nNodes[i]); d_weights[i].resize(nNodes[i]);
    
    if(i>0)
    {
      for(int j=0; j<nNodes[i]; j++)
      {
        weights[i][j].resize(nNodes[i-1], 0.0); d_weights[i][j].resize(nNodes[i-1], 0.0);
      }
    }
  }
  
  features_AD.resize(nNodes[0]);
  beta_AD.resize(nNodes[nNodes.size()-1]);
  nodes_AD.resize(nNodes.size());
  for(int i=0; i<nNodes.size(); i++)
  {
    nodes_AD[i].resize(nNodes[i], 0.0);
  }
  
  rans    = rans_val;
  rans_AD = rans_AD_val;
}
