#include "SCTModuleHit.hh"
#include "SCTModuleGeometry.hh"
#include "DetectorConstruction.hh"

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>

#include "G4SDManager.hh"
#include "G4RunManager.hh"


using LineEquation = std::pair<G4ThreeVector, G4ThreeVector>;

LineEquation getStripLineEquation(SCTModuleHit* hit) {
    // Get the rotation matrix of the strip volume
    // G4RotationMatrix* rotation = stripVolume->GetRotation();
    G4RotationMatrix rotation = hit->GetStripRotation();
    
    std::cout << "Rotation matrix of the strip in SpacePoint.hh: " << std::endl;
    std::cout << rotation << std::endl;

    // The strip is aligned along the local X-axis before rotation
    G4ThreeVector localDir(0, 1, 0);
    // Rotate to get the global direction
    G4ThreeVector globalDir = rotation* localDir;
    globalDir = globalDir.unit(); // Normalize

    // Get a point on the strip (the center of the volume)
    // G4ThreeVector pointOnStrip = stripVolume->GetTranslation();
    G4ThreeVector pointOnStrip = hit->GetStripCentre();

    // The line equation can be represented as: P(t) = pointOnStrip + t * globalDir
    return std::make_pair(pointOnStrip, globalDir);
}


std::set<SCTModuleHit> makeSpacePoints(SCTModuleHitCollection* hitCollection) {
  
    std::set<SCTModuleHit> spacePoints;

    // Get the detector construction to access geometry information
    auto *runManager = G4RunManager::GetRunManager();
    auto detector = (DetectorConstruction*) (runManager->GetUserDetectorConstruction());

    G4int nModules = detector->GetNlayers() * detector->GetModulesPerLayer();

    std::map<G4int,  // layer id  
    std::map<G4int, // module id
    std::pair<std::set<SCTModuleHit*>, std::set<SCTModuleHit*>>>> hitsInModule; // First: side1, Second: side2

    // Group hits by layer, module and side
    for (auto hit : *hitCollection->GetVector())
    {
        G4int layer_id = hit->GetLayerNumber();
        G4int module_id = hit->GetModuleNumber();
        G4int side = hit->GetStripSide();
        if (side == 0) {
            hitsInModule[layer_id][module_id].first.insert(hit);
            G4cout << "Hit on layer " << layer_id << ", module " << module_id << ", side " << side << ": " << *hit << G4endl; 
        } else if (side == 1) {
            hitsInModule[layer_id][module_id].second.insert(hit);
            G4cout << "Hit on layer " << layer_id << ", module " << module_id << ", side " << side << ": " << *hit << G4endl;   
        }
    }

    // Create space points from paired hits on both sides of the module
    for (G4int layer_id = 0; layer_id < detector->GetNlayers(); ++layer_id) {
        for (G4int module_id = 0; module_id < detector->GetModulesPerLayer(); ++module_id) {
            auto& side1Hits = hitsInModule[layer_id][module_id].first;
            auto& side2Hits = hitsInModule[layer_id][module_id].second;

            std::vector<LineEquation> stripLineEqs1;
            std::vector<LineEquation> stripLineEqs2;

           // Get the physical volumes for hits on each side
            for (auto* hit1 : side1Hits) {
                // stripLineEqs1.push_back(getStripLineEquation(hit1->GetPhysVol()));
                stripLineEqs1.push_back(getStripLineEquation(hit1));
            }
            for (auto* hit2 : side2Hits) {
                // stripLineEqs2.push_back(getStripLineEquation(hit2->GetPhysVol()));
                stripLineEqs2.push_back(getStripLineEquation(hit2));
            }
    
            for (const auto& top_strip: stripLineEqs1) {
                for (const auto& bottom_strip: stripLineEqs2)
                {
                    // Calculate the closest points between the two lines
                    G4ThreeVector p1 = top_strip.first;
                    G4ThreeVector d1 = top_strip.second;
                    G4ThreeVector p2 = bottom_strip.first;
                    G4ThreeVector d2 = bottom_strip.second;
                    
                    // Translate the two planes to be in midpoint of the two planes in z
                    G4double average_zpos = (p1.z() + p2.z())/2;
                    p1.setZ(average_zpos);
                    p2.setZ(average_zpos);
                    
                    G4ThreeVector r = p1 - p2;
                    G4double a = d1.dot(d1);
                    G4double b = d1.dot(d2);
                    G4double c = d2.dot(d2);
                    G4double e = d1.dot(r);
                    G4double f = d2.dot(r);

                    G4double denom = a * c - b * b;
                    if (std::abs(denom) < 1e-6) {
                        std::cout << "Lines are parallel, cannot form space point." << std::endl;
                        continue; // Lines are parallel, skip
                    }

                    G4double t1 = (b * f - c * e) / denom;
                    G4double t2 = (a * f - b * e) / denom;

                    G4ThreeVector closestPoint1 = p1 + t1 * d1;
                    G4ThreeVector closestPoint2 = p2 + t2 * d2;

                    // Midpoint between the closest points
                    G4ThreeVector spacePointPos = 0.5 * (closestPoint1 + closestPoint2);

                    // Create a new SCTModuleHit for the space point
                    SCTModuleHit spacePointHit;
                    spacePointHit.SetPosition(spacePointPos.x(), spacePointPos.y(), spacePointPos.z());
                    spacePointHit.SetLayerNumber(layer_id);
                    spacePointHit.SetModuleNumber(module_id);
                    // Additional properties can be set as needed

                    G4cout << "Created space point at " << spacePointPos << " for layer " << layer_id << ", module " << module_id << std::endl;

                    spacePoints.insert(spacePointHit);
                }
            }
        }
    }


    return spacePoints;
}
