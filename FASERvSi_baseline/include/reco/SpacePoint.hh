#pragma once
#include "SCTModuleHit.hh"
#include "SCTModuleGeometry.hh"
#include "DetectorConstruction.hh"

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"

using LineEquation = std::pair<G4ThreeVector, G4ThreeVector>;

LineEquation getStripLineEquation(SCTModuleHit* hit) {
    // Get the rotation matrix of the strip volume
    // G4RotationMatrix* rotation = stripVolume->GetRotation();
    G4RotationMatrix rotation = hit->GetStripRotation();
    
    // std::cout << "Rotation matrix of the strip in SpacePoint.hh: " << std::endl;
    // std::cout << rotation << std::endl;

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


static G4double cross2D(const G4TwoVector& v1, const G4TwoVector& v2) {
    return v1.x() * v2.y() - v1.y() * v2.x();
}


std::set<SCTModuleHit> makeSpacePoints(SCTModuleHitCollection* hitCollection) {
    
    G4cout << "Creating space points from hits..." << G4endl;
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
                // stripLineEqs1.push_back(getStripLineEquation(hit1));
                stripLineEqs1.push_back(hit1->GetStripEnds());

                G4cout << "Top: " << hit1->GetStripEnds().first << ", "  << hit1->GetStripEnds().second << std::endl;
                SCTModuleHit* spacePointHit1 = new SCTModuleHit();
                spacePointHit1->SetPosition(hit1->GetStripEnds().first[0], hit1->GetStripEnds().first[1], hit1->GetStripEnds().first[2]);
                spacePointHit1->SetLayerNumber(1*layer_id);
                spacePointHit1->SetModuleNumber(2*module_id);
                spacePointHit1->SetColour(G4Colour::White());
                spacePoints.insert(*spacePointHit1);
                std::cout << "Inserting space point hit at " << hit1->GetStripEnds().first[0] << ", " << hit1->GetStripEnds().first[1] << ", " << hit1->GetStripEnds().first[2] << std::endl;

                SCTModuleHit* spacePointHit2 = new SCTModuleHit();
                spacePointHit2->SetPosition(hit1->GetStripEnds().second[0], hit1->GetStripEnds().second[1], hit1->GetStripEnds().second[2]);
                spacePointHit2->SetLayerNumber(3*layer_id);
                spacePointHit2->SetModuleNumber(4*module_id);
                spacePointHit2->SetColour(G4Colour::White());
                spacePoints.insert(*spacePointHit2);
                std::cout << "Inserting space point hit at " << hit1->GetStripEnds().second[0] << ", " << hit1->GetStripEnds().second[1] << ", " << hit1->GetStripEnds().second[2] << std::endl;
                
            }
            for (auto* hit2 : side2Hits) {
                // stripLineEqs2.push_back(getStripLineEquation(hit2->GetPhysVol()));
                // stripLineEqs2.push_back(getStripLineEquation(hit2));

                stripLineEqs2.push_back(hit2->GetStripEnds());
                G4cout << "Bottom: " << hit2->GetStripEnds().first << ", "  << hit2->GetStripEnds().second << std::endl;
                SCTModuleHit* spacePointHit1 = new SCTModuleHit();
                spacePointHit1->SetPosition(hit2->GetStripEnds().first[0], hit2->GetStripEnds().first[1], hit2->GetStripEnds().first[2]);
                spacePointHit1->SetLayerNumber(5*layer_id);
                spacePointHit1->SetModuleNumber(6*module_id);
                spacePointHit1->SetColour(G4Colour::White());
                spacePoints.insert(*spacePointHit1);
                std::cout << "Inserting space point hit at " << hit2->GetStripEnds().first[0] << ", " << hit2->GetStripEnds().first[1] << ", " << hit2->GetStripEnds().first[2] << std::endl;

                SCTModuleHit* spacePointHit2 = new SCTModuleHit();
                spacePointHit2->SetPosition(hit2->GetStripEnds().second[0], hit2->GetStripEnds().second[1], hit2->GetStripEnds().second[2]);
                spacePointHit2->SetLayerNumber(7*layer_id);
                spacePointHit2->SetModuleNumber(8*module_id);
                spacePointHit2->SetIsReco(true);
                spacePointHit2->SetColour(G4Colour::White());
                spacePoints.insert(*spacePointHit2);
                std::cout << "Inserting space point hit at " << hit2->GetStripEnds().second[0] << ", " << hit2->GetStripEnds().second[1] << ", " << hit2->GetStripEnds().second[2] << std::endl;
            }
    
            G4double stripHalfLength = SCTModuleGeometry::stripLength() / 2;

            // https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect,
            for (const auto& top_strip : stripLineEqs1) {
                for (const auto& bottom_strip : stripLineEqs2) {
                    
                    //* Line segment is parameterised as the p -> p + s
                    //* p is the start of the strip and s is the 2vector that I add to p to get to the end of the strip
                    //* i.e. p + s gives me the coordinates to the other end of the script
                    G4TwoVector q = G4TwoVector(top_strip.first.x(), top_strip.first.y());         // top strip end 1
                    G4TwoVector s = G4TwoVector(top_strip.second.x(), top_strip.second.y()) - q;
                    G4cout << "q = " << q << " s = " << s << std::endl;

                    G4TwoVector p = G4TwoVector(bottom_strip.first.x(), bottom_strip.first.y());       // bottom strip end 1
                    G4TwoVector r = G4TwoVector(bottom_strip.second.x(), bottom_strip.second.y()) - p;
                    G4cout << "p = " << p << " r = " << r << std::endl;
                    
                    G4double t = cross2D(q - p, s) / cross2D(r, s);
                    G4double u = cross2D(q - p, r) / cross2D(r, s);

                    G4cout << "t = "  << t << ", u = " << u << std::endl;

                    if (cross2D(r, s) == 0) {
                        if (cross2D(q - p, r) !=0)
                        {
                        std::cout << "Lines are parallel, cannot form space point." << std::endl;
                        continue;
                        }
                        else{
                            std::cout << "Lines are collinear, cannot form space point." << std::endl;
                            continue;
                        }
                    }
                    else if (t < 0 || t > 1 || u < 0 || u > 1) {
                        std::cout << "No intersection within the strip bounds, skipping." << std::endl;
                        continue;
                    }

                    G4TwoVector intersection2D = p + t * r;
                    G4double xpos = intersection2D.x();
                    G4double ypos = intersection2D.y();
                    G4double zpos = (top_strip.first[2] +  bottom_strip.first[2]) / 2;

                    G4ThreeVector spacePointPos = G4ThreeVector(xpos, ypos, zpos);
                    SCTModuleHit spacePointHit;
                    spacePointHit.SetPosition(xpos, ypos, zpos);
                    spacePointHit.SetLayerNumber(layer_id);
                    spacePointHit.SetModuleNumber(module_id);
                    spacePointHit.SetColour(G4Colour::Blue());

                    G4cout << "Created space point at " << spacePointPos << " for layer " << layer_id << ", module " << module_id << std::endl;

                    spacePoints.insert(spacePointHit);
                }
            }
        }
    }


    return spacePoints;
}
