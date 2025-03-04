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
    fdetHeight = 300*2*mm;

    ftungstenThickness = 13.92*mm;
    // fSCTThickness = 0.5*mm; // 7.08*mm;
    fSCTThickness = 7.08*mm;
    fSCTSideThickness = 1.00*mm; // 1 mm separation between front and back side of the SCT module
    fnumSCTLayers = 50;

    fstripWidth = 0.08*mm; // 80 mircons
    fstripsPerSide = 768;
    fstripStereoAngle = 0.04; // 40 miliradians
    fmoduleWidth = 125*mm;
    fmoduleHeight = 31.25*mm;
    fsideThickness = 0.1*mm;
    fsideSeparation = 0.1*mm;


}

DetectorParameters* DetectorParameters::Get()
{
  if (!me)
    me = new DetectorParameters();
  return me; 
}

