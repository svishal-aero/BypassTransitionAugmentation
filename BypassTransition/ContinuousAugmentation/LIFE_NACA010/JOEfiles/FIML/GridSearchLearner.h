#ifndef GRIDSEARCHLEARNER_H
#define GRIDSEARCHLEARNER_H

#include <cstdio>
#include "Learner.h"

class GridSearchVar
{
  public:

    double                    min;
    double                    max;
    int                     nBins;
    std::vector<double> binBounds;

  public:

    GridSearchVar(void) { }

    ~GridSearchVar() {}

    void init(double min_val, double max_val, int nBins_val)
    {
      min = min_val; max = max_val; nBins = nBins_val;
      binBounds.resize(nBins+1);
      binBounds[0] = min;
      double binSize = (max-min) / nBins;
      for(int i=0; i<nBins; i++) binBounds[i] = binBounds[i-1] + binSize;
      binBounds[nBins] = max;
    }

    void init(std::vector<double> binBounds_val)
    {
      binBounds = binBounds_val;
      nBins = binBounds.size()-1;
      min = binBounds[0]; max = binBounds[nBins];
    }

    int operator()(double input)
    {
      for(int i=0; i<nBins; i++)
      {
        if(input>=binBounds[i] && input<=binBounds[i+1])
          return i;
      }
      return -1;
    }
};

class GridSearchLearner : public Learner
{
  public:

    std::vector<GridSearchVar> vars;

  public:

    void addGridSearchVar(double min, double max, int nBins)
    {
      vars.push_back(GridSearchVar());
      vars[vars.size()-1].init(min, max, nBins);
    }

    void initParamsAndDerivs(void)
    {
      int arrSize = 1;
      for(int i=0; i<vars.size(); i++)
        arrSize *= vars[i].nBins;
      params.resize(arrSize);
      derivs.resize(arrSize);
    }

    double getLinearizedFunction(std::vector<double> features)
    {
      int temp, index;
      temp = vars[0](features[0]);
      if(temp>=0) index = temp; else return 1.0;

      for(int i=1; i<vars.size(); i++)
      {
        index *= vars[i].nBins;
        temp = vars[i](features[i]);
        if(temp>=0) index += temp; else return 1.0;
      }
      return params[index];
    }
};

#endif
