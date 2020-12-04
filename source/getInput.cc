#include "getInput.h"

std::string trim(const std::string &str,
                 const std::string &whitespace)
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

std::string removeComments(const std::string &str,
                           const std::string &begChar)
{
    const auto strBegin = str.find_first_of(begChar);
    if (strBegin == std::string::npos)
        return str; // no comment

    if (strBegin == 0)
        return ""; // no content

    return str.substr(0, strBegin);
}

int getInput(const std::string filepath, Parameters &setting)
{
    std::ifstream fileptr;
    fileptr.open(filepath, std::ios::in);
    if (!fileptr.is_open())
    {
        std::cout << "can't open file " << filepath << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> lines;
    std::vector<std::string> newlines;
    std::string tag[] = {"CoincidencesFile:", "LO2EnergyFile:",
                         "OutputFile:", "DetNum:",
                         "CHANNEL:", "EnableLO2Energy:",
                         "SigmaT:", "Radius:",
                         "TrueAzimuth:", "TrueElevation:"};
    while (getline(fileptr, line))
    {
        lines.push_back(line);
    }
    u_int lineindex(0);
    while (lineindex < lines.size())
    {
        line = lines[lineindex];
        line = trim(line); //remove all leading and trailing " \t".
        if (line.length() > 1)
        {
            //remove all comments (starts with "#")
            line = removeComments(line, "#");
            // remove all white spaces
            line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
        }
        if (line.length() > 1)
        {
            //testout << line << std::endl;
            newlines.push_back(line);
        }
        lineindex++;
    }
    lines.clear();

    for (std::size_t i = 0; i < newlines.size(); i++)
    {
        line = newlines[i];
        std::size_t taglength(0);
        // Input file 1, list of coincident neutron events
        if (line.find(tag[0]) == 0)
        {
            taglength = tag[0].length();
            setting.CoincidenceFile = line.substr(taglength, std::string::npos);
        }
        // Input file 2, Light output to energy
        else if (line.find(tag[1]) == 0)
        {
            taglength = tag[1].length();
            setting.LO2EnergyFile = line.substr(taglength, std::string::npos);
        }
        // output file
        else if (line.find(tag[2]) == 0)
        {
            taglength = tag[2].length();
            setting.OutImageFile = line.substr(taglength, std::string::npos);
        }
        // Number of detectors
        else if (line.find(tag[3]) == 0)
        {
            taglength = tag[3].length();
            setting.DetNum = stoi(line.substr(taglength, std::string::npos));
        }
        // Detector info
        else if (line.find(tag[4]) == 0)
        {
            taglength = tag[4].length();
            // setting.DetNum = stoi(line.substr(taglength, std::string::npos));
            std::string strTemp = line.substr(taglength, std::string::npos);
            std::stringstream ss(strTemp);
            std::vector<std::string> substrs;
            while (ss.good())
            {
                std::string substr;
                getline(ss, substr, ',');
                substrs.push_back(substr);
            }
            DetectorInfo newDetector;
            // newDetector.ChannelNum = stoi(substrs[0]);
            newDetector.pos_x = stod(substrs[1]);
            newDetector.pos_y = stod(substrs[2]);
            newDetector.pos_z = stod(substrs[3]);
            newDetector.sigma_x = stod(substrs[4]);
            newDetector.sigma_y = stod(substrs[5]);
            newDetector.sigma_z = stod(substrs[6]);
            setting.Detectors[stoi(substrs[0])] = newDetector;
            // setting.Detectors.push_back(newDetector);
        }
        // Enable lightout to energy
        else if (line.find(tag[5]) == 0)
        {
            taglength = tag[5].length();
            setting.LO2Energy = stoi(line.substr(taglength, std::string::npos));
        }
        // Time resolution
        else if (line.find(tag[6]) == 0)
        {
            taglength = tag[6].length();
            setting.SigmaT = stod(line.substr(taglength, std::string::npos));
        }
        // Radius of projection sphere
        else if (line.find(tag[7]) == 0)
        {
            taglength = tag[7].length();
            setting.Radius = stod(line.substr(taglength, std::string::npos));
        }
        // true azimuthal angle and uncertainty in degrees, -180 ~180
        else if (line.find(tag[8]) == 0)
        {
            taglength = tag[8].length();
            std::string strTemp = line.substr(taglength, std::string::npos);
            std::stringstream ss(strTemp);
            while (ss.good())
            {
                std::string substr;
                getline(ss, substr, ',');
                // substr.erase(remove_if(substr.begin(), substr.end(), isspace), substr.end());
                setting.TrueAzimuth.push_back(stod(substr));
            }
        }
        // true elvation angle and uncertainty in degrees, -90 ~90
        else if (line.find(tag[9]) == 0)
        {
            taglength = tag[9].length();
            std::string strTemp = line.substr(taglength, std::string::npos);
            std::stringstream ss(strTemp);
            while (ss.good())
            {
                std::string substr;
                getline(ss, substr, ',');
                // substr.erase(remove_if(substr.begin(), substr.end(), isspace), substr.end());
                setting.TrueElevation.push_back(stod(substr));
            }
        }
        // add more options here.
    }
    return 0;
}