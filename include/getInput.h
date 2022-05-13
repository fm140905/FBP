#ifndef GETINPUT_H
#define GETINPUT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <map>

/**
 * @brief Struct for stroring the detector positions.
 * 
 */
struct DetectorInfo
{
    // // channel number
    // u_int ChannelNum;
    // x coordinate, cm
    double pos_x;
    // y coordinate, cm
    double pos_y;
    // z coordinate, cm
    double pos_z;
    // uncertainty associated with x coordinate, cm
    double sigma_x;
    // uncertainty associated with y coordinate, cm
    double sigma_y;
    // uncertainty associated with z coordinate, cm
    double sigma_z;
};

/**
 * @brief Struct storing parameters read from the input file.
 * 
 */
struct Parameters
{
    // Input file 1, list of coincident neutron events
    std::string CoincidenceFile;
    // Input file 2, light output to energy
    std::string LO2EnergyFile;
    // Output image 
    std::string OutImageFile;
    // Number of detectors
    u_int DetNum;
    // detector positions and associated uncertainties
    std::map<u_int, DetectorInfo> Detectors;
    // std::vector<DetectorInfo> Detectors;
    bool LO2Energy;
    double SigmaT;
    double Radius;
    std::vector<double> TrueAzimuth;
    std::vector<double> TrueElevation;
};

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t");

std::string removeComments(const std::string& str,
                 const std::string& begChar);

/**
 * @brief read input file and save the parameters
 * 
 * @param filepath path to the input file
 * @param setting struct storing the input parameters
 * @return int 
 */
int getInput(const std::string filepath, Parameters& setting);

#endif