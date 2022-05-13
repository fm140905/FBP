add_library(getInput STATIC getInput.cc)

add_library(cpimaging STATIC cpimaging.cc)
target_link_libraries(cpimaging PUBLIC getInput ROOT::Gpad)
# target_link_libraries(loadEvents PUBLIC getInput ROOT::Core ROOT::Gpad 
                    #   ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt 
                    #   ROOT::MathCore ROOT::Matrix ROOT::MultiProc 
                    #   ROOT::Net ROOT::Physics ROOT::Postscript 
                    #   ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps 
                    #   ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree)