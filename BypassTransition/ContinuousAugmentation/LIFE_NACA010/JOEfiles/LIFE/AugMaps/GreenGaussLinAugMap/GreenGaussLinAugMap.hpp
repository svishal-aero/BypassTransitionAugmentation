#ifndef GreenGaussLinAugMap_hpp
#define GreenGaussLinAugMap_hpp

#include <cstdio>
#include "../GenAugMap.hpp"

template<const int nFtrs>
class GreenGaussLinAugMap : public GenAugMap<nFtrs>
{
    protected:

        std::array<int,    nFtrs> nCells;
        std::array<int,    nFtrs> offsets;
        std::array<double, nFtrs> minFtrs;
        std::array<double, nFtrs> maxFtrs;
        std::array<double, nFtrs> h;

    public:

        GreenGaussLinAugMap() : GenAugMap<nFtrs>() { }

        inline void init(
            const std::string& name_,
            const std::array<int, nFtrs>& nCells_,
            const std::array<double, nFtrs>& minFtrs_,
            const std::array<double, nFtrs>& maxFtrs_
        )
        {
            this->name = name_;
            nCells = nCells_;
            minFtrs = minFtrs_;
            maxFtrs = maxFtrs_;

            int nParams = 1;
            for(auto& nCellsFtr : nCells)
                nParams *= nCellsFtr;
            this->params.resize(nParams);
            this->dJ_dParams.resize(nParams);

            offsets[nFtrs-1] = 1;
            for(int iFtr=nFtrs-2; iFtr>=0; iFtr--)
                offsets[iFtr] = nCells[iFtr+1] * offsets[iFtr+1];

            for(int iFtr=0; iFtr<nFtrs; iFtr++)
                h[iFtr] = (maxFtrs[iFtr]-minFtrs[iFtr]) / nCells[iFtr];
        }

        inline double calculate_value(
            const std::array<double, nFtrs>& ftrs,
            double* dBeta_dFtrs_export = nullptr
        )
        {
            static std::array<int, nFtrs> iCellFtr;
            static std::array<double, nFtrs> dBeta_dFtrs;
            static int iCell;

            iCell = 0;
            for(int iFtr=0; iFtr<nFtrs; iFtr++)
            {
                iCellFtr[iFtr] = (int)((ftrs[iFtr]-minFtrs[iFtr])/h[iFtr]);
                iCell += offsets[iFtr] * iCellFtr[iFtr];
            }

            double value = this->params[iCell];

            for(int iFtr=0; iFtr<nFtrs; iFtr++)
            {
                int iPrev = iCellFtr[iFtr]==0 ? iCell : iCell-offsets[iFtr];
                int iNext = iCellFtr[iFtr]==nCells[iFtr]-1 ? iCell : iCell+offsets[iFtr];
                dBeta_dFtrs[iFtr] = 0.5 * (this->params[iNext]-this->params[iPrev]) / h[iFtr];
                if(dBeta_dFtrs_export!=nullptr)
                    dBeta_dFtrs_export[iFtr] = dBeta_dFtrs[iFtr];
                value += dBeta_dFtrs[iFtr] * (ftrs[iFtr] -(iCellFtr[iFtr]+0.5)*h[iFtr]);
            }

            return value;
        }

        inline void calculate_sens(
            const std::array<double, nFtrs>& ftrs,
            const double& dJ_dBeta
        )
        {
            static std::array<int, nFtrs> iCellFtr;
            static std::array<double, nFtrs> dBeta_dFtrs;
            static int iCell;

            iCell = 0;
            for(int iFtr=0; iFtr<nFtrs; iFtr++)
            {
                iCellFtr[iFtr] = (int)((ftrs[iFtr]-minFtrs[iFtr])/h[iFtr]);
                iCell += offsets[iFtr] * iCellFtr[iFtr];
            }

            for(int iFtr=0; iFtr<nFtrs; iFtr++)
            {
                int iPrev = iCellFtr[iFtr]==0 ? iCell : iCell-offsets[iFtr];
                int iNext = iCellFtr[iFtr]==nCells[iFtr]-1 ? iCell : iCell+offsets[iFtr];
                this->dJ_dParams[iCell] += dJ_dBeta;
                this->dJ_dParams[iNext] += dJ_dBeta*0.5/h[iFtr];
                this->dJ_dParams[iPrev] -= dJ_dBeta*0.5/h[iFtr];
            }
        }
};

#endif
