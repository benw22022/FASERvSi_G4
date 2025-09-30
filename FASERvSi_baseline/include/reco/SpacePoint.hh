#pragma once
#include "SCTModuleHit.hh"
#include "TruthHit.hh"
#include "SCTModuleGeometry.hh"
#include "DetectorConstruction.hh"

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"


namespace SpacePointUtils {

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

    // 2D cross product helper
    static G4double cross2D(const G4TwoVector& v1, const G4TwoVector& v2) {
        return v1.x() * v2.y() - v1.y() * v2.x();
    }


    void makeSpacePoints(SCTModuleHitCollection* inputHitCollection, SCTModuleHitCollection* outputHitCollection) {

        //* Get the detector construction to access geometry information
        auto *runManager = G4RunManager::GetRunManager();
        auto detector = (DetectorConstruction*) (runManager->GetUserDetectorConstruction());
        G4int nModules = detector->GetNlayers() * detector->GetModulesPerLayer();

        //* Use nested map to group hits by layer, module and side
        std::map<G4int,  // layer id  
        std::map<G4int,  // module id
        std::pair<std::vector<SCTModuleHit*>, std::vector<SCTModuleHit*>>>> hitsInModule; // First: side1, Second: side2

        //* Group hits by layer, module and side
        for (auto hit : *inputHitCollection->GetVector())
        {
            G4int layer_id = hit->GetLayerNumber();
            G4int module_id = hit->GetModuleNumber();
            G4int side = hit->GetStripSide();
            if (side == 0) {
                hitsInModule[layer_id][module_id].first.push_back(hit);
                // G4cout << "Hit on layer " << layer_id << ", module " << module_id << ", side " << side << ": " << *hit << G4endl; 
            } else if (side == 1) {
                hitsInModule[layer_id][module_id].second.push_back(hit);
                // G4cout << "Hit on layer " << layer_id << ", module " << module_id << ", side " << side << ": " << *hit << G4endl;   
            }
        }

        //* Create space points by computing the 2d intercept of each strip by treating as line segements
        //* Loop through every possible pair of strips to find all space points on module
        G4int nReconstructedSpacePoints{0};
        for (G4int layer_id = 0; layer_id < detector->GetNlayers(); ++layer_id) {
            for (G4int module_id = 0; module_id < detector->GetModulesPerLayer(); ++module_id) {
                auto& side1Hits = hitsInModule[layer_id][module_id].first;
                auto& side2Hits = hitsInModule[layer_id][module_id].second;

                if (side1Hits.empty() || side2Hits.empty()) {
                    // G4cout << "No hits on one side of module " << module_id << " in layer " << layer_id << ", skipping." << G4endl;
                    continue; // No hits on one side, skip
                }

                // Store as sets as you can get multiple hits on the same strip
                std::set<LineEquation> stripLineEqs1;
                std::set<LineEquation> stripLineEqs2;

            //* Get the ends of each strip that was hit
                for (auto* hit1 : side1Hits) {
                    // stripLineEqs1.push_back(hit1->GetStripEnds());
                    stripLineEqs1.insert(hit1->GetStripEnds());

                    // G4cout << "layer: " << layer_id << " module: " << module_id <<  "  Top: " << hit1->GetStripEnds().first << ", "  << hit1->GetStripEnds().second << std::endl;
                    // SCTModuleHit* spacePointHit1 = new SCTModuleHit();
                    // spacePointHit1->SetPosition(hit1->GetStripEnds().first[0], hit1->GetStripEnds().first[1], hit1->GetStripEnds().first[2]);
                    // spacePointHit1->SetLayerNumber(1*layer_id);
                    // spacePointHit1->SetModuleNumber(2*module_id);
                    // spacePointHit1->SetColour(G4Colour::White());
                    // spacePoints.insert(*spacePointHit1);
                    // std::cout << "Inserting space point hit at " << hit1->GetStripEnds().first[0] << ", " << hit1->GetStripEnds().first[1] << ", " << hit1->GetStripEnds().first[2] << std::endl;

                    // SCTModuleHit* spacePointHit2 = new SCTModuleHit();
                    // spacePointHit2->SetPosition(hit1->GetStripEnds().second[0], hit1->GetStripEnds().second[1], hit1->GetStripEnds().second[2]);
                    // spacePointHit2->SetLayerNumber(3*layer_id);
                    // spacePointHit2->SetModuleNumber(4*module_id);
                    // spacePointHit2->SetColour(G4Colour::White());
                    // spacePoints.insert(*spacePointHit2);
                    // std::cout << "Inserting space point hit at " << hit1->GetStripEnds().second[0] << ", " << hit1->GetStripEnds().second[1] << ", " << hit1->GetStripEnds().second[2] << std::endl;
                    
                }
                for (auto* hit2 : side2Hits) {
                    // stripLineEqs2.push_back(hit2->GetStripEnds());
                    stripLineEqs2.insert(hit2->GetStripEnds());
        
                    // G4cout << "layer: " << layer_id << " module: " << module_id <<  "  Bottom: " << hit2->GetStripEnds().first << ", "  << hit2->GetStripEnds().second << std::endl;
                    // SCTModuleHit* spacePointHit1 = new SCTModuleHit();
                    // spacePointHit1->SetPosition(hit2->GetStripEnds().first[0], hit2->GetStripEnds().first[1], hit2->GetStripEnds().first[2]);
                    // spacePointHit1->SetLayerNumber(5*layer_id);
                    // spacePointHit1->SetModuleNumber(6*module_id);
                    // spacePointHit1->SetColour(G4Colour::White());
                    // spacePoints.insert(*spacePointHit1);
                    // std::cout << "Inserting space point hit at " << hit2->GetStripEnds().first[0] << ", " << hit2->GetStripEnds().first[1] << ", " << hit2->GetStripEnds().first[2] << std::endl;

                    // SCTModuleHit* spacePointHit2 = new SCTModuleHit();
                    // spacePointHit2->SetPosition(hit2->GetStripEnds().second[0], hit2->GetStripEnds().second[1], hit2->GetStripEnds().second[2]);
                    // spacePointHit2->SetLayerNumber(7*layer_id);
                    // spacePointHit2->SetModuleNumber(8*module_id);
                    // spacePointHit2->SetIsReco(true);
                    // spacePointHit2->SetColour(G4Colour::White());
                    // spacePoints.insert(*spacePointHit2);
                    // std::cout << "Inserting space point hit at " << hit2->GetStripEnds().second[0] << ", " << hit2->GetStripEnds().second[1] << ", " << hit2->GetStripEnds().second[2] << std::endl;
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

                        G4TwoVector p = G4TwoVector(bottom_strip.first.x(), bottom_strip.first.y());       // bottom strip end 1
                        G4TwoVector r = G4TwoVector(bottom_strip.second.x(), bottom_strip.second.y()) - p;
                        
                        G4double denom = cross2D(r, s);
                        G4double t = cross2D(q - p, s) / denom;
                        G4double u = cross2D(q - p, r) / denom;
                        
                        // G4cout << "layer: " << layer_id << " module: " << module_id <<  "  [q = " << q << ", s = " << s <<",  t = "  << t << ", u = " << u << "]" << G4endl;

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
                            // std::cout << "No intersection within the strip bounds, skipping." << std::endl;
                            continue;
                        }

                        G4TwoVector intersection2D = p + t * r;
                        G4double xpos = intersection2D.x();
                        G4double ypos = intersection2D.y();
                        G4double zpos = (top_strip.first[2] +  bottom_strip.first[2]) / 2;

                        G4ThreeVector spacePointPos = G4ThreeVector(xpos, ypos, zpos);
                        SCTModuleHit* spacePointHit = new SCTModuleHit();
                        spacePointHit->SetPosition(xpos, ypos, zpos);
                        spacePointHit->SetLayerNumber(layer_id);
                        spacePointHit->SetModuleNumber(module_id);
                        spacePointHit->SetColour(G4Colour::Gray());
                        spacePointHit->SetIsReco(true);
                        // G4cout << "Created space point at " << spacePointPos << " for layer " << layer_id << ", module " << module_id << std::endl;

                        outputHitCollection->insert(spacePointHit);
                        nReconstructedSpacePoints++;
                    }
                }
            }
        }

        // G4cout << "Spacepoint maker: Made " << nReconstructedSpacePoints << " space points" << G4endl;
    }


    void truthMatchSpacePoints(SCTModuleHitCollection* spacePointCollection, TruthHitCollection* truthHitCollection) {
        //* For each space point, find the closest truth hit in the same layer and module
        std::vector<G4int> matchedTruthHitIndices;


        for (auto truthHit : *truthHitCollection->GetVector()) {
            G4int th_layer = truthHit->GetLayerNumber();
            G4int th_module = truthHit->GetModuleNumber();
            G4ThreeVector th_pos = G4ThreeVector(truthHit->GetX(), truthHit->GetY(), truthHit->GetZ());
            
            SCTModuleHit* bestMatch = nullptr;
            // Use distance squared since we don't need to bother with the sqrt when computing the Euclidean distance
            G4double bestDistance2 = pow(1 * mm, 2); //std::numeric_limits<G4double>::max();//pow(2 * mm, 2);
                
            for (auto spacePoint : *spacePointCollection->GetVector()) {
                G4int sp_layer = spacePoint->GetLayerNumber();
                G4int sp_module = spacePoint->GetModuleNumber();
                G4ThreeVector sp_pos = G4ThreeVector(spacePoint->GetX(), spacePoint->GetY(), spacePoint->GetZ());
                
                if (th_layer == sp_layer && th_module == sp_module) {
                    G4double distance2 = (sp_pos - th_pos).mag2();
                    if (distance2 < bestDistance2) {
                        bestDistance2 = distance2;
                        bestMatch = spacePoint;
                    }
                }
                
            if (bestMatch) {
                // Link the space point to the truth hit
                bestMatch->SetTrackID(truthHit->GetTrackID());
                bestMatch->SetParentID(truthHit->GetParentID());
                bestMatch->SetPDGID(truthHit->GetPDGID());
                bestMatch->SetEnergy(truthHit->GetEnergy());
                bestMatch->SetCharge(truthHit->GetCharge());
                bestMatch->SetPx(truthHit->GetPx());
                bestMatch->SetPy(truthHit->GetPy());
                bestMatch->SetPz(truthHit->GetPz());
                bestMatch->SetMass(truthHit->GetMass());
                bestMatch->SetTruthHitID(truthHit->GetTruthHitID());
                bestMatch->SetIsTruthMatched(true);
                bestMatch->SetColour(G4Colour::Blue());

                // G4cout << "Matching truth hit (trackID " << truthHit->GetTrackID() << ", layer " << th_layer << ", module " << th_module << ") to space point at (" 
                //        << bestMatch->GetX() << ", " << bestMatch->GetY() << ", " << bestMatch->GetZ() << ")" << G4endl;

                // Optionally, set a flag or store the matched TruthHit pointer
                // spacePoint->SetMatchedTruthHit(bestMatch);
            }
            
            }

        }


    //     matchedTruthHitIndices.reserve(spacePointCollection->GetVector()->size());
    //     for (auto spacePoint : *spacePointCollection->GetVector()) {
    //         G4int sp_layer = spacePoint->GetLayerNumber();
    //         G4int sp_module = spacePoint->GetModuleNumber();
    //         G4ThreeVector sp_pos = G4ThreeVector(spacePoint->GetX(), spacePoint->GetY(), spacePoint->GetZ());

    //         TruthHit* bestMatch = nullptr;
    //         G4double bestDistance2 = std::numeric_limits<G4double>::max();

    //         G4int th_index{0};
    //         for (auto truthHit : *truthHitCollection->GetVector()) {
    //             G4int th_layer = truthHit->GetLayerNumber();
    //             G4int th_module = truthHit->GetModuleNumber();
    //             G4ThreeVector th_pos = G4ThreeVector(truthHit->GetX(), truthHit->GetY(), truthHit->GetZ());

    //             if (th_layer == sp_layer && th_module == sp_module) {
    //                 G4double distance2 = (sp_pos - th_pos).mag2();
    //                 if (distance2 < bestDistance2) {
    //                     bestDistance2 = distance2;
    //                     bestMatch = truthHit;
    //                 }
    //             }
    //             th_index++;
    //         }

    //         if (bestMatch) {
    //             // Link the space point to the truth hit
    //             spacePoint->SetTrackID(bestMatch->GetTrackID());
    //             spacePoint->SetParentID(bestMatch->GetParentID());
    //             spacePoint->SetPDGID(bestMatch->GetPDGID());
    //             spacePoint->SetEnergy(bestMatch->GetEnergy());
    //             spacePoint->SetCharge(bestMatch->GetCharge());
    //             spacePoint->SetPx(bestMatch->GetPx());
    //             spacePoint->SetPy(bestMatch->GetPy());
    //             spacePoint->SetPz(bestMatch->GetPz());
    //             spacePoint->SetMass(bestMatch->GetMass());
    //             spacePoint->SetIsTruthMatched(true);
    //             // Optionally, set a flag or store the matched TruthHit pointer
    //             // spacePoint->SetMatchedTruthHit(bestMatch);
    //         }
    //     }
    }
} // namespace SpacePointUtils