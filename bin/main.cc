#include <vector>
#include <ctime>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <tuple>

/* ROOT */
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
// #include "TF1.h"
// #include "TGraph.h"
// #include "TMultiGraph.h"

#include "getInput.h"
#include "cpimaging.h"

int main(int argc, char *argv[])
{
    struct stat st = {0};
    // create directories to save the output files and pics, works in Linux
    if (stat("output", &st) == -1) {
        mkdir("output", 0700);
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    
    TApplication *myapp = new TApplication("myApp", 0, 0);

    // read input file
    const std::string inputFilePath("input/input.txt");
    Parameters setting;
    getInput(inputFilePath, setting);

    // load coincident events
    std::cout << "Reading neutron coincidences from " << setting.CoincidenceFile << std::endl;
    std::vector<neutronPair> neutronPairs;
    loadNeutronPairs(setting, neutronPairs);

    // create image
    const UInt_t thetaBins = 179;
    const UInt_t phiBins = 360;
    std::vector<std::vector<Float_t>> imageFBP(thetaBins, std::vector<Float_t>(phiBins, 0));
    cpFBP(setting, neutronPairs, imageFBP);

    // plot
    TH2D *sino = new TH2D("ROI", " ; Azimuth; Elevation", phiBins, -180, 180, thetaBins, -89.5, 89.5);
    TCanvas *canvas4 = new TCanvas("THCanvas","Sinocanvas", 1000, 500);
    getImage(setting, sino, canvas4, imageFBP);
    canvas4->SaveAs(setting.OutImageFile.c_str());

    auto endTime = std::chrono::high_resolution_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    std::cout << "Ctrl + C to exit..." << std::endl;
    myapp->Run();
    return 0;
}