#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include "adolc.h"

#define NN_LINEAR  0
#define NN_SIGMOID 1
#define NN_TANH    2
#define NN_SWISH   3

class NN
{
  private:

    std::string      name;
    std::vector<int> act;

    std::vector<double> features;
    std::vector<double> beta;
    std::vector<std::vector<double>> nodes, d_nodes;
    std::vector<std::vector<double>> biases, d_biases;
    std::vector<std::vector<std::vector<double>>> weights, d_weights;

    std::vector<adouble> features_AD;
    std::vector<adouble> beta_AD;
    std::vector<std::vector<adouble>> nodes_AD;

    RansTurb    *rans;
    RansTurb_AD *rans_AD;

  public:

    #include "NN_core/init.h"
    #include "NN_core/actFun.h"
    #include "NN_core/forwprop.h"
    #include "NN_core/backprop.h"
    #include "NN_core/readParams.h"
    #include "NN_core/writeParams.h"

    #include "calcFeatures.h"
    #include "calcFilter.h"

    void calculate(int icv)
    {
      calcFeatures(icv); for(int j=0; j<nodes[0].size(); j++) nodes[0][j] = features[j];
      forwprop();
      for(int j=0; j<nodes[nodes.size()-1].size(); j++) beta[j] = nodes[nodes.size()-1][j] * calcFilter(icv);
    }

    void calculate_AD(int icv)
    {
      calcFeatures_AD(icv); for(int j=0; j<nodes_AD[0].size(); j++) nodes_AD[0][j] = features_AD[j];
      forwprop_AD();
      for(int j=0; j<nodes_AD[nodes_AD.size()-1].size(); j++) beta_AD[j] = nodes_AD[nodes_AD.size()-1][j] * calcFilter_AD(icv);
    }

    void evalSens(int icv, int iOutput, double dJdbeta)
    {
      calcFeatures(icv); for(int j=0; j<nodes[0].size(); j++) nodes[0][j] = features[j];
      forwprop();
      for(int j=0; j<nodes[nodes.size()-1].size(); j++) d_nodes[nodes.size()-1][j] = 0.0;
      d_nodes[nodes.size()-1][iOutput] = dJdbeta * calcFilter(icv);
      backprop();
    }
};
