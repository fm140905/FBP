#if !defined(CPIMAGING_H)
#define CPIMAGING_H

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TEllipse.h"

#include "omp.h"

#include "getInput.h"

class neutronPair
{
private:
    /* data */
public:
    float_t delT=0;
    // Energy deposited in the first detector, MeV 
    float_t Energy0=0;
    UShort_t CH0=0;
    UShort_t CH1=0;
    // constructor
    neutronPair();
    neutronPair(float_t delt, Float_t erg0, UShort_t ch0, UShort_t ch1): delT(delt), Energy0(erg0), CH0(ch0), CH1(ch1) {};
    // ~neutronPair();

};

class LO2EnergyTable
{
private:
    std::map<Float_t, Float_t> table;
public:
    LO2EnergyTable(/* args */);
    LO2EnergyTable(const std::string filepath);
    // ~LO2Energy();
    Float_t getEnergy(const Float_t LO);
};

// LO2Energy::LO2Energy(/* args */)
// {
// }

// LO2Energy::~LO2Energy()
// {
// }

// neutronPair::neutronPair(/* args */)
// {
// }

// neutronPair::~neutronPair()
// {
// }


// void mergeSort(const Parameters& setting, const std::vector< std::vector< Event > >& multiCHEvents, std::vector<neutronPair>& tuples);
void loadNeutronPairs(const Parameters settings, std::vector<neutronPair>& tuples);

void cpFBP(const Parameters settings, const std::vector<neutronPair>& tuples, std::vector<std::vector<Float_t>>& ImageFBP);

void plotImage(const Parameters settings, TH2D *histo, TCanvas *canvas, const std::vector<std::vector<Float_t>>& imageFBP);

bool saveImage2Txt(const std::string fpath, const std::vector<std::vector<Float_t>>& imageFBP);

#endif // CPIMAGING_H