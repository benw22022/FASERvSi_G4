#include "DetectorParameters.hh"
#include "G4SystemOfUnits.hh"

DetectorParameters *DetectorParameters::me = 0;

DetectorParameters::DetectorParameters()
{
    fexpHall_x = 0.5*m;
    fexpHall_y = 0.5*m;
    fexpHall_z = 2*m;
    
    ftargetStartPosZ = -0.75*m;
    fdetWidth = 25*cm;
    fdetHeight = 30*cm;

    ftungstenThickness = 9*mm;
    fSCTThickness = 7.08*mm;

    fnumSCTLayers = 132;

}

DetectorParameters* DetectorParameters::Get()
{
  if (!me)
    me = new DetectorParameters();
  return me; 
}

