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

/**
 * @brief A pair of neutron coincidences.
 * 
 */
class NeutronPair
{
private:
    /* data */
public:
    // time of flight between two detections, ns.
    float_t delT=0;
    // Energy deposited in the first detector, MeV.
    float_t Energy0=0;
    // detector number of the first detector.
    UShort_t CH0=0;
    // detector number of the second detector.
    UShort_t CH1=0;
    // constructor
    NeutronPair();
    NeutronPair(float_t delt, Float_t erg0, UShort_t ch0, UShort_t ch1): delT(delt), Energy0(erg0), CH0(ch0), CH1(ch1) {};
    // ~NeutronPair();

};

/**
 * @brief Convert light output (MeVee) to energy deposited (MeV).
 * 
 */
class LO2EnergyTable
{
private:
    // map storing the LO conversion function. 
    // Key is the LO, MeVee.
    // Value is the corresponding energy deposited, MeV.
    std::map<Float_t, Float_t> table;
public:
    LO2EnergyTable(/* args */);
    LO2EnergyTable(const std::string filepath);
    // ~LO2Energy();
    /**
     * @brief Get the Energy deposited correponding to the LO value.
     * 
     * @param LO Light output, MeVee.
     * @return Float_t Energy in MeV.
     */
    Float_t getEnergy(const Float_t LO);
};

// void mergeSort(const Parameters& setting, const std::vector< std::vector< Event > >& multiCHEvents, std::vector<neutronPair>& tuples);

/**
 * @brief Load the neutron coincidences from file.
 * 
 * @param settings Input parameters from the input file.
 * @param tuples List of neutron coincidences.
 */
void loadNeutronPairs(const Parameters settings, std::vector<NeutronPair>& tuples);

/**
 * @brief Perform back projection image reconstruction.
 * 
 * @param settings Input parameters.
 * @param tuples List of neutron coincidences.
 * @param ImageFBP Reconstructed image.
 */
void cpFBP(const Parameters settings, const std::vector<NeutronPair>& tuples, std::vector<std::vector<Float_t>>& ImageFBP);

/**
 * @brief Plot the reconstructed image using ROOT.
 * 
 * @param settings Input parameters.
 * @param histo ROOT 2D histogram storing the image.
 * @param canvas Plot canvas.
 * @param imageFBP Reconstructed image.
 */
void plotImage(const Parameters settings, TH2D *histo, TCanvas *canvas, const std::vector<std::vector<Float_t>>& imageFBP);

/**
 * @brief Save the reconstructed image to a text file.
 * 
 * @param fpath Path to the output file.
 * @param imageFBP Reconstructed image.
 * @return true if data saving was successful.
 * @return false if data saving failed.
 */
bool saveImage2Txt(const std::string fpath, const std::vector<std::vector<Float_t>>& imageFBP);

#endif // CPIMAGING_H