
#ifndef G4d2oNeutronGun_h
#define G4d2oNeutronGun_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4d2oPrimaryGeneratorAction.hh"
#include "inputVariables.hh"
#include "G4d2oNeutrinoAlley.hh"
#include "G4d2oDetector.hh"

#include "TF1.h"
#include <RQ_OBJECT.h>

class G4ParticleGun;
class G4Event;

class G4d2oNeutronGun : public G4d2oPrimaryGeneratorAction
{
    RQ_OBJECT("G4d2oNeutronGun")

public:
    G4d2oNeutronGun();
    ~G4d2oNeutronGun();
    
public:
    virtual void GeneratePrimaries(G4Event* anEvent);
    
    virtual G4double GetSourceEnergy() { return sourceEnergy; }
    virtual G4ThreeVector GetOriginalPosition() { return SourcePosition; }
    
protected:
    // ---- Beam mode flags ----
    // Set ONE to true, the others to false
    static constexpr G4bool fTopBeamMode = false;    // true = beam from top (Z)
    static constexpr G4bool fSideBeamMode = true;    // true = beam from side (-X)
    static constexpr G4bool fVolumeMode = false;     // true = uniform inside D₂O

private:
    G4int totalEvents;
    G4double cthetaRange1, cthetaRange2;
    G4int theEventNum;
    G4int irand;

    // Outer H2O geometry (hardcoded)
    G4ThreeVector h2oGlobalCenter;
    G4double h2oGlobalRadius;
    G4double h2oGlobalHalfHeight;
    G4bool h2oVolumeFound;

    void GenerateVolumeNeutron();
    void GenerateTopBeamNeutron();
    void GenerateSideBeamNeutron();
};

#endif
