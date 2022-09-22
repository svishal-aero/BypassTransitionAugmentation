#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../PwMultiLinAugMap.hpp"

int main()
{
  MPI_Init(NULL, NULL);
  
  PwMultiLinAugMap<2> augMap;
  augMap.init("testMap", {11, 11}, {0, 0}, {1.000001, 1.000001});

  for(int i=0; i<11; i++) for(int j=0; j<11; j++)
  {
    augMap.params[i*11+j] = exp(-0.25*((i-3)*(i-3) + (j-8)*(j-8)));
  }
  
  FILE *fp = fopen("plot.dat", "w");
  FILE *gp = fopen("plot_check.dat", "w");
  
    for(int i=0; i<101; i++) for(int j=0; j<101; j++)
    {
      double val = augMap.calculate_value({0.01*i, 0.01*j});
      fprintf(fp, "%le %le %le\n", 0.01*i, 0.01*j, val);
      fprintf(gp, "%le %le %le\n", 0.01*i, 0.01*j, exp(-0.01*0.25*((i-30)*(i-30) + (j-80)*(j-80))));
    }
  
  fclose(gp);
  fclose(fp);
  
  MPI_Finalize();
}
