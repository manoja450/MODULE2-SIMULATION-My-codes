#include "G4d2oNeutronGun.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"

#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"

G4d2oNeutronGun::G4d2oNeutronGun()
    : G4d2oPrimaryGeneratorAction()
{
    G4cout << "\tConstructing G4d2oNeutronGun..." ;

    input = inputVariables::GetIVPointer();
    irand = input->GetRandomStatus();
    if (irand == 0) gRandom->SetSeed(1);
    if (irand == 1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    if (particleGun) {
        delete particleGun;
    }
    particleGun = new G4ParticleGun(1);
    G4String particleName = "neutron";
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));

    SourcePosition.setX(0.0);
    SourcePosition.setY(0.0);
    SourcePosition.setZ(0.0);

    G4double thetaAngle1 = 0.0;
    G4double thetaAngle2 = 180.0;
    cthetaRange1 = cos(thetaAngle1 * deg);
    cthetaRange2 = cos(thetaAngle2 * deg);

    // ---- Hardcode outer H2O geometry (from your detector) ----
    h2oGlobalCenter.setX(0.0);
    h2oGlobalCenter.setY(110.73384 * cm);
    h2oGlobalCenter.setZ(-80.74526 * cm);
    h2oGlobalRadius = 44.45 * cm;
    h2oGlobalHalfHeight = 95.3135 * cm;
    h2oVolumeFound = true;

    G4cout << "  H2O center: (" << h2oGlobalCenter.x()/cm << ", "
           << h2oGlobalCenter.y()/cm << ", " << h2oGlobalCenter.z()/cm << ") cm" << G4endl;
    G4cout << "  H2O radius: " << h2oGlobalRadius/cm << " cm" << G4endl;
    G4cout << "  H2O half-height: " << h2oGlobalHalfHeight/cm << " cm" << G4endl;
    G4cout << "  Top Beam mode:   " << (fTopBeamMode ? "ON" : "OFF") << G4endl;
    G4cout << "  Side Beam mode:  " << (fSideBeamMode ? "ON" : "OFF") << G4endl;
    G4cout << "  Volume mode:     " << (fVolumeMode ? "ON" : "OFF") << G4endl;
    G4cout << "done." << G4endl;
}

G4d2oNeutronGun::~G4d2oNeutronGun()
{
    G4cout << "Deleting G4d2oNeutronGun...";
    G4cout << "done." << G4endl;
}

void G4d2oNeutronGun::GenerateVolumeNeutron()
{
    // Original uniform generation inside D₂O
    G4double maxRadius = 33.655;
    G4double randvarn1 = maxRadius * sqrt(G4UniformRand());
    G4double randvarn2 = 2.0 * TMath::Pi() * G4UniformRand();
    G4double randvarn3 = 89.9275 * (-1.0 + 2.0 * G4UniformRand());

    G4double x = randvarn1 * cos(randvarn2);
    G4double y = 110.73384 + randvarn1 * sin(randvarn2);
    G4double z = -80.74526 + randvarn3;

    SourcePosition.setX(x * cm);
    SourcePosition.setY(y * cm);
    SourcePosition.setZ(z * cm);

    G4double ctheta = G4UniformRand() * (cthetaRange1 - cthetaRange2) + cthetaRange2;
    G4double phi = 2.0 * TMath::Pi() * G4UniformRand();
    G4double stheta = sqrt(1.0 - ctheta * ctheta);

    G4double px = stheta * cos(phi);
    G4double py = stheta * sin(phi);
    G4double pz = ctheta;

    initDir.setX(px);
    initDir.setY(py);
    initDir.setZ(pz);
    initDir = initDir.unit();

    particleGun->SetParticleMomentumDirection(initDir);
    sourceEnergy = 0.025 * eV;
    particleGun->SetParticleEnergy(sourceEnergy);
}

void G4d2oNeutronGun::GenerateTopBeamNeutron()
{
    // Use hardcoded H2O geometry
    G4double planeZ = h2oGlobalCenter.z() + h2oGlobalHalfHeight - 0.5 * cm; // ~14.07 cm
    G4double R = h2oGlobalRadius - 0.5 * cm;

    G4double r = R * sqrt(G4UniformRand());
    G4double phi = 2.0 * TMath::Pi() * G4UniformRand();

    G4double x = h2oGlobalCenter.x() + r * cos(phi);
    G4double y = h2oGlobalCenter.y() + r * sin(phi);

    SourcePosition.setX(x);
    SourcePosition.setY(y);
    SourcePosition.setZ(planeZ);

    // Downward hemisphere (pz < 0)
    G4double ctheta = -G4UniformRand();
    G4double stheta = sqrt(1.0 - ctheta * ctheta);
    G4double pphi = 2.0 * TMath::Pi() * G4UniformRand();

    initDir.setX(stheta * cos(pphi));
    initDir.setY(stheta * sin(pphi));
    initDir.setZ(ctheta);
    initDir = initDir.unit();

    particleGun->SetParticleMomentumDirection(initDir);
    sourceEnergy = 0.025 * eV;
    particleGun->SetParticleEnergy(sourceEnergy);
}

void G4d2oNeutronGun::GenerateSideBeamNeutron()
{
    G4double R = h2oGlobalRadius - 0.5 * cm;
    G4double offset = 1.0 * cm;
    G4double Xpos = h2oGlobalCenter.x() - (R - offset);

    G4double dx = Xpos - h2oGlobalCenter.x();
    G4double maxY = sqrt(R * R - dx * dx);

    G4double y = h2oGlobalCenter.y() + (2.0 * G4UniformRand() - 1.0) * maxY;
    G4double z = h2oGlobalCenter.z() + (2.0 * G4UniformRand() - 1.0) * (h2oGlobalHalfHeight - 0.5 * cm);

    SourcePosition.setX(Xpos);
    SourcePosition.setY(y);
    SourcePosition.setZ(z);

    // Forward hemisphere (px > 0)
    G4double ctheta = G4UniformRand();
    G4double stheta = sqrt(1.0 - ctheta * ctheta);
    G4double pphi = 2.0 * TMath::Pi() * G4UniformRand();

    initDir.setX(ctheta);
    initDir.setY(stheta * cos(pphi));
    initDir.setZ(stheta * sin(pphi));
    initDir = initDir.unit();

    particleGun->SetParticleMomentumDirection(initDir);
    sourceEnergy = 0.025 * eV;
    particleGun->SetParticleEnergy(sourceEnergy);
}

void G4d2oNeutronGun::GeneratePrimaries(G4Event* anEvent)
{
    theEventNum = anEvent->GetEventID();

    if (fTopBeamMode) {
        GenerateTopBeamNeutron();
    } else if (fSideBeamMode) {
        GenerateSideBeamNeutron();
    } else if (fVolumeMode) {
        GenerateVolumeNeutron();
    } else {
        G4cout << "WARNING: No beam mode selected! Using volume mode." << G4endl;
        GenerateVolumeNeutron();
    }

    particleGun->SetParticlePosition(SourcePosition);

    if (theEventNum < 10) {
        G4Navigator* Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        G4VPhysicalVolume* volume = Navigator->LocateGlobalPointAndSetup(SourcePosition);
        G4String volName = volume ? volume->GetName() : "UNKNOWN";
        G4String matName = volume ? volume->GetLogicalVolume()->GetMaterial()->GetName() : "UNKNOWN";

        G4cout << "Event " << theEventNum;
        if (fTopBeamMode) G4cout << " (TOP BEAM)";
        else if (fSideBeamMode) G4cout << " (SIDE BEAM)";
        else G4cout << " (VOLUME MODE)";
        G4cout << ": Neutron E = " << sourceEnergy / eV << " eV"
               << ", Pos = (" << SourcePosition.x() / cm
               << ", " << SourcePosition.y() / cm
               << ", " << SourcePosition.z() / cm << ") cm"
               << ", Volume = " << volName
               << ", Material = " << matName
               << G4endl;
    }

    particleGun->GeneratePrimaryVertex(anEvent);
}
