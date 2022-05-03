#include "cpimaging.h"

LO2EnergyTable::LO2EnergyTable(const std::string filepath)
{
    // read data from file
    std::ifstream fileptr;
    fileptr.open(filepath, std::ios::in);
    if (!fileptr.is_open())
    {
        // std::cout << "can't open file: " << filepath << std::endl;
        // exit(1);
        std::string errMessage = "can't open file: " + filepath;
        throw std::runtime_error(errMessage);
    }
    std::string line;
    // Read one line at a time into the variable line:
    Int_t lineIndex = 1;
    Float_t a=0,b=0;
    while(std::getline(fileptr, line))
    {
        if (lineIndex > 2)
        {
            std::stringstream lineStream(line);
            lineStream >> a >> b; 
            table[a] = b;
        }
        lineIndex ++;
    }
}

Float_t LO2EnergyTable::getEnergy(const Float_t LO)
{
    auto ub = table.upper_bound(LO); // find the first entry that > LO
    if (ub!=table.cbegin() && ub!=table.cend())
    {
        auto lb = std::prev(ub);
        Float_t value = ((LO-lb->first) * (ub->second) + (ub->first - LO) * (lb->second)) / (ub->first - lb->first);
        return value;
    }
    else
    {
        // all entries are <= LO because LO is too large, or all entries are > LO 
        throw std::runtime_error("LO is too small or too large.");
    }
}

// void mergeSort(const Parameters& setting, const std::vector< std::vector< Event > >& multiCHEvents, std::vector<neutronPair >& tuples)
// {
//     // Lightoutput to energy
//     LO2EnergyTable LO2Energy("input/LO2Energy.csv");

//     const int nChannels = multiCHEvents.size();
//     std::vector<Event> sortedEvents;
//     for (int i = 0; i < multiCHEvents.size(); i++)
//     {
//         sortedEvents.insert(sortedEvents.end(), multiCHEvents[i].begin(), multiCHEvents[i].end());
//     }
//     std::sort(sortedEvents.begin(),sortedEvents.end());
//     // std::vector<std::pair<Long64_t, Float_t> > pairs;
//     const Long64_t timeWindow = 1000 * setting.TimeWindow;
//     Long64_t delT=0;
//     for (int i = 0; i < sortedEvents.size()-1; i++)
//     {
//         if (sortedEvents[i].isBad || sortedEvents[i+1].isBad)
//         {
//             continue;
//         }
//         delT = sortedEvents[i+1].timeStampHeader - sortedEvents[i].timeStampHeader;
//         if (delT < timeWindow)
//         {
//             // tuples.push_back(std::make_tuple(delT, sortedEvents[i].energy, sortedEvents[i].CHnum, sortedEvents[i+1].CHnum));
//             tuples.push_back(neutronPair(float_t(delT)/1E3, LO2Energy.getEnergy(sortedEvents[i].energy/1000.0), sortedEvents[i].CHnum, sortedEvents[i+1].CHnum));
//         }
//     }
// }
void loadNeutronPairs(const Parameters settings, std::vector<neutronPair>& tuples)
{
    // Lightoutput to energy
    LO2EnergyTable LO2Energy(settings.LO2EnergyFile);
    // read data from file
    std::ifstream fileptr;
    fileptr.open(settings.CoincidenceFile, std::ios::in);
    if (!fileptr.is_open())
    {
        // std::cout << "can't open file: " << filepath << std::endl;
        // exit(1);
        std::string errMessage = "can't open file: " + settings.CoincidenceFile;
        throw std::runtime_error(errMessage);
    }
    std::string line;
    // Read one line at a time into the variable line:
    Int_t lineIndex = 1;
    Float_t delT=0,Energy0=0;
    UShort_t CH0=0, CH1=0;
    if (settings.LO2Energy)
    {
        while(std::getline(fileptr, line))
        {
            if (lineIndex > 1) // ignore first line
            {
                std::stringstream lineStream(line);
                lineStream >> delT >> Energy0 >> CH0 >> CH1;
                // if ((CH0 == 0 || CH1 == 4) ||(CH0 == 4 || CH1 == 0))
                // {
                //     continue;
                // }
                tuples.push_back(neutronPair(delT, LO2Energy.getEnergy(Energy0), CH0, CH1));
            }
            lineIndex ++;
        }
    }
    else
    {
        while(std::getline(fileptr, line))
        {
            if (lineIndex > 1) // ignore first line
            {
                std::stringstream lineStream(line);
                lineStream >> delT >> Energy0 >> CH0 >> CH1;
                // if ((CH0 == 0 || CH1 == 4) ||(CH0 == 4 || CH1 == 0))
                // {
                //     continue;
                // }
                tuples.push_back(neutronPair(delT, Energy0, CH0, CH1));
            }
            lineIndex ++;
        }
    }
}

void cpFBP(const Parameters settings, const std::vector<neutronPair>& tuples, std::vector<std::vector<Float_t>>& ImageFBP)
{
    // Altitude angle
    // int thetaBins = 180;
    int thetaBins = ImageFBP.size();
    Float_t dtheta = 179/float_t(thetaBins);
    Float_t thetaMin = -89;//
    Float_t thetaMax = 89;//
    // Azimuthal angle
    // int phiBins = 360;
    int phiBins = ImageFBP[0].size();
    Float_t dphi = 360/float_t(phiBins);
    Float_t phiMin = -180 + dphi / 2;//
    Float_t phiMax = 180 - dphi / 2;//

    // creata a sphere and pixelate it
    Float_t radius = settings.Radius; // cm
    float_t theta(thetaMin * M_PI / 180);
    float_t phi(phiMin * M_PI / 180);
    std::vector<std::vector<std::vector<Float_t>>> xbs(thetaBins, std::vector<std::vector<Float_t>>(phiBins, std::vector<Float_t>(3,0)));
    for (int i = 0; i < thetaBins; i++)
    {
        phi = phiMin * M_PI / 180;
        for (int j = 0; j < phiBins; j++)
        {
            xbs[i][j][0] = radius * std::cos(theta) * std::cos(phi);
            xbs[i][j][1] = radius * std::cos(theta) * std::sin(phi);
            xbs[i][j][2] = radius * std::sin(theta);
            phi += dphi * M_PI / 180;
        }
        theta += dtheta * M_PI / 180;
    }
        
    u_int CH0(0);
    u_int CH1(0);
    float_t ErgToF(0);
    float_t dx10,dy10,dz10;
    float_t dxb0,dyb0,dzb0;
    float_t dr10square(0);
    float_t drb0square(0);
    float_t dr10dotdrb0(0);
    float dx10dotsigma0(0);
    float dx10dotsigma1(0);

    float_t alpha(0);
    float_t beta(0);
    // sigma
    const float_t sigmat=settings.SigmaT; // ns
    float_t sigmabeta(0);
    float_t sigmaalpha(0);

    const float_t k(0.5*939.56542052/898.755719); // Mev / (cm/ns)^2
    float temp(0);
    float temp2(0);
    float temp3(0);
    int i = 0;

    std::cout << "Projecting cones onto the designated spherical surface..."<<std::endl;
    float progress=0.0;
    u_int currentConeNum(0);
    u_int totalConeNum(tuples.size());
    for (neutronPair const &event : tuples)
    {
        std::vector<std::vector<Float_t>> newCone(thetaBins, std::vector<Float_t>(phiBins, 0));
        // show progress bar
        int barWidth = 70;
        currentConeNum++;
        // if (currentConeNum % 10 == 0)
        {
            progress = float(currentConeNum) / float(totalConeNum);
            std::cout << "[";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
        }
        
        CH0 = event.CH0;
        CH1 = event.CH1;
        sigmaalpha=0;
        dx10 = settings.Detectors.at(event.CH1).pos_x - settings.Detectors.at(event.CH0).pos_x;
        dy10 = settings.Detectors.at(event.CH1).pos_y - settings.Detectors.at(event.CH0).pos_y;
        dz10 = settings.Detectors.at(event.CH1).pos_z - settings.Detectors.at(event.CH0).pos_z;
        dr10square = (dx10 * dx10 + dy10 * dy10 + dz10 * dz10); // cm^2
        // dx10 dot sigma_0
        dx10dotsigma0 = (dx10 * settings.Detectors.at(event.CH0).sigma_x + 
                         dy10 * settings.Detectors.at(event.CH0).sigma_y + 
                         dz10 * settings.Detectors.at(event.CH0).sigma_z); // cm^2
        // dx10 dot sigma_1
        dx10dotsigma1 = (dx10 * settings.Detectors.at(event.CH1).sigma_x + 
                         dy10 * settings.Detectors.at(event.CH1).sigma_y + 
                         dz10 * settings.Detectors.at(event.CH1).sigma_z); // cm^2
        // Energy of the scattered neutron, MeV
        ErgToF = k*(dr10square)/std::pow(event.delT,2); // MeV
        // cosine squared
        alpha = ErgToF / (ErgToF + event.Energy0);
        // uncertainty
        temp = alpha / (event.Energy0 + ErgToF);
        // partial alpha over partial t
        sigmaalpha += std::pow(2*event.Energy0 / event.delT*temp*sigmat,2);
        // partial alpha over x0,y0,z0
        sigmaalpha += std::pow(2*event.Energy0*temp*dx10dotsigma0/dr10square,2);
        // partial alpha over x1,y1,z1
        sigmaalpha += std::pow(2*event.Energy0*temp*dx10dotsigma1/dr10square,2);
        // sigmaalpha = std::sqrt(sigmaalpha);

        double pixelSum(0);
        for (int i =0; i < thetaBins; i++)
        {
            // theta = thetaMin + i * dtheta;
            sigmabeta=0;
            
            for (int j = 0; j < phiBins; j++)
            {
                dxb0 = xbs[i][j][0] - settings.Detectors.at(event.CH0).pos_x;
                dyb0 = xbs[i][j][1] - settings.Detectors.at(event.CH0).pos_y;
                dzb0 = xbs[i][j][2] - settings.Detectors.at(event.CH0).pos_z;
                drb0square = dxb0 * dxb0 + dyb0 * dyb0 + dzb0 * dzb0;
                dr10dotdrb0 = dx10 * dxb0 + dy10 * dyb0 + dz10 * dzb0;
                // if(dr10dotdrb0 < 0)
                {
                    // phi = phiMin + j * dphi;
                    beta = std::pow(dr10dotdrb0, 2) / (dr10square * drb0square);
                    // uncertainty
                    // partial beta over x1,y1,z1
                    temp = beta / dr10dotdrb0;
                    temp2 = beta / dr10square;
                    temp3 = beta / drb0square;
                    sigmabeta += (std::pow((2*dxb0*temp-2*dx10*temp2) * settings.Detectors.at(event.CH1).sigma_x,2)+
                                  std::pow((2*dyb0*temp-2*dy10*temp2) * settings.Detectors.at(event.CH1).sigma_y,2)+
                                  std::pow((2*dzb0*temp-2*dz10*temp2) * settings.Detectors.at(event.CH1).sigma_z,2));
                    // partial beta over x0,y0,z0
                    sigmabeta += (std::pow((-2*(dxb0+dx10)*temp+2*dx10*temp2+2*dxb0*temp3) * settings.Detectors.at(event.CH0).sigma_x,2)+
                                  std::pow((-2*(dyb0+dy10)*temp+2*dy10*temp2+2*dyb0*temp3) * settings.Detectors.at(event.CH0).sigma_y,2)+
                                  std::pow((-2*(dzb0+dz10)*temp+2*dz10*temp2+2*dzb0*temp3) * settings.Detectors.at(event.CH0).sigma_z,2));
                    // sigmabeta = std::sqrt(sigmabeta);

                    // pixel value
                    newCone[i][j] = std::exp(-1.0*(beta-alpha)*(beta-alpha) / (2*(sigmaalpha + sigmabeta)));
                    pixelSum += newCone[i][j];
                }
            }
        }

        for (int i =0; i < thetaBins; i++)
        {
            for (int j = 0; j < phiBins; j++)
            {
                ImageFBP[i][j] += newCone[i][j] / pixelSum;
            }
        }
    }
    std::cout << std::endl;
}

void getImage(const Parameters settings, TH2D *histo, TCanvas *canvas, const std::vector<std::vector<Float_t>>& imageFBP)
{
    // // canvas1->cd(4)->SetFrameFillColor(TColor::GetColorPalette(0));
    int phiBins = histo->GetNbinsX();
    int thetaBins = histo->GetNbinsY();
    // for (int i=240;i<280;i++)
    // {
    //     for (int j=100;j<140;j++)
    //     {
    //         histo->SetBinContent(i+1,j+1,imageFBP[j][i]);
    //     }
    // }
    for (int i=0;i<phiBins;i++)
    {
        for (int j=0;j<thetaBins;j++)
        {
            histo->SetBinContent(i+1,j+1,imageFBP[j][i]);
        }
    }
    gStyle->SetOptStat(0);
    histo->GetZaxis()->SetLabelSize(0.02);
    histo->Draw("z aitoff");
    
    // grid
    float conv=M_PI/180; // I am also aware of TMath::DegToRad() and TMath::Pi() which could be used there...
    float la, lo, x, y, z;
    
    const int Nl = 19; // Number of drawn latitudes
    const int NL = 19; // Number of drawn longitudes
    int M  = 30;

    TGraph  *latitudes[Nl];
    TGraph  *longitudes[NL];

    for (int j=0;j<Nl;++j) {
        latitudes[j]=new TGraph();
        la = -90+180/(Nl-1)*j;
        for (int i=0;i<M+1;++i) {
            lo = -180+360/M*i;
            z  = sqrt(1+cos(la*conv)*cos(lo*conv/2));
            x  = 180*cos(la*conv)*sin(lo*conv/2)/z;
            y  = 90*sin(la*conv)/z;
            latitudes[j]->SetPoint(i,x,y);
        }
    } 
    
    for (int j=0;j<NL;++j) {
        longitudes[j]=new TGraph();
        lo = -180+360/(NL-1)*j;
        for (int i=0;i<M+1;++i) {
            la = -90+180/M*i;
            z  = sqrt(1+cos(la*conv)*cos(lo*conv/2));
            x  = 180*cos(la*conv)*sin(lo*conv/2)/z;
            y  = 90*sin(la*conv)/z;
            longitudes[j]->SetPoint(i,x,y);
        }
    }
    // Draw the grid	
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000);
    pad2->SetFillColor(0);
    pad2->SetBorderSize(0);
    pad2->SetFrameBorderMode(0);
    pad2->SetFrameLineColor(0); 
    pad2->SetFrameBorderMode(0);
    pad2->Draw();
    pad2->cd();
    Double_t ymin = -89.5;
    Double_t ymax = 89.5;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = -180;
    Double_t xmax = 180;
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right

    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);

    for (int j=0;j<Nl;++j) latitudes[j]->Draw("l");
    for (int j=0;j<NL;++j) longitudes[j]->Draw("l");
    // // // canvas1->cd(4)->SetLogz();

    // // // draw source 1
    // // TEllipse* el1 = new TEllipse(86.5276, 0, 5., 5.);
    // // el1->SetFillColor(0);
    // // el1->SetFillStyle(0);
    // // el1->SetLineColor(4);
    // // el1->Draw("SAME");
    // draw source 2
    TEllipse* el2 = new TEllipse(settings.TrueAzimuth[0], settings.TrueElevation[0], 
                                 settings.TrueAzimuth[1], settings.TrueElevation[1]);
    el2->SetFillColor(0);
    el2->SetFillStyle(0);
    el2->SetLineColor(4);
    el2->Draw("SAME");
    canvas->Draw();
}