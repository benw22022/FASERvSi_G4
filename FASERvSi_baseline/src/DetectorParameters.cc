#include "DetectorParameters.hh"
#include "G4SystemOfUnits.hh"

DetectorParameters *DetectorParameters::me = 0;

DetectorParameters::DetectorParameters()
{
    fexpHall_x = 1*m;
    fexpHall_y = 1*m;
    fexpHall_z = 6*m;
    
    ftargetStartPosZ = 550*mm;
    fdetWidth = 250*mm;
    fdetHeight = 300*mm;

    ftungstenThickness = 0.9*mm;
    fSCTThickness = 0.5*mm; // 7.08*mm;

    fnumSCTLayers = 132;

}

DetectorParameters* DetectorParameters::Get()
{
  if (!me)
    me = new DetectorParameters();
  return me; 
}

