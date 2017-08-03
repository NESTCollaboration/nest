//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The NEST program is intended for use with the Geant4 software,   *
// * which is copyright of the Copyright Holders of the Geant4        *
// * Collaboration. This additional software is copyright of the NEST *
// * development team. As such, it is subject to the terms and        *
// * conditions of both the Geant4 License, included with your copy   *
// * of Geant4 and available at http://cern.ch/geant4/license, as     *
// * well as the NEST License included with the download of NEST and  *
// * available at http://nest.physics.ucdavis.edu/                    *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutions, nor the agencies providing financial support for   *
// * this work make any representation or warranty, express or        *
// * implied, regarding this software system, or assume any liability *
// * for its use. Please read the pdf license or view it online       *
// * before download for the full disclaimer and lack of liability.   *
// *                                                                  *
// * This code implementation is based on work by Peter Gumplinger    *
// * and his fellow collaborators on Geant4 and is distributed with   *
// * the express written consent of the Geant4 collaboration. By      *
// * using, copying, modifying, or sharing the software (or any work  *
// * based on the software) you agree to acknowledge use of both NEST *
// * and Geant4 in resulting scientific publications, and you         *
// * indicate your acceptance of all the terms and conditions of the  *
// * licenses, which must always be included with this code.          *
// ********************************************************************
//
//
////////////////////////////////////////////////////////////////////////

#include "G4ParticleTypes.hh" //lets you refer to G4OpticalPhoton, etc.
#include "G4EmProcessSubType.hh" //lets you call this process Scintillation
#include "G4Version.hh" //tells you what Geant4 version you are running
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
#include "G4ProductionCuts.hh"
#include "G4Electron.hh"
#include <cmath>
#include "NESTProc.hh"
#include "NESTStackingAction.hh"
#include "G4RandomDirection.hh"

using namespace NEST;

NESTProc::NESTProc(const G4String& processName,
		           G4ProcessType type)
      : G4VRestDiscreteProcess(processName, type), electron_cut_index(G4ProductionCuts::GetIndex(G4Electron::Definition()))
{
        
        SetProcessSubType(fScintillation);
	
        fTrackSecondariesFirst = false;
	
        if (verboseLevel>0) {
	  G4cout << GetProcessName() << " is created " << G4endl;
        }
}

NESTProc::~NESTProc(){} //destructor needed to avoid linker error



G4Track* MakePhoton(NEST::Vertex vertex, bool exciton) {
    // Determine polarization of new photon
    G4ParticleMomentum photonMomentum(G4RandomDirection());
    G4ThreeVector perp = photonMomentum.cross(G4RandomDirection());
    G4ThreeVector photonPolarization = perp.unit();

    G4double PhotMean = 6.97*eV; G4double PhotWidth = 0.23*eV;
    G4double sampledEnergy = G4RandGauss::shoot(PhotMean, PhotWidth);
    G4DynamicParticle* aQuantum =
            new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
            photonMomentum);
    aQuantum->SetPolarization(photonPolarization.x(),
            photonPolarization.y(),
            photonPolarization.z());
    aQuantum->SetKineticEnergy(sampledEnergy);
    //calculate time
    
    double tau1 = G4RandGauss::shoot(3.1*ns,.7*ns); //err from wgted avg.
    double tau3 = G4RandGauss::shoot(24.*ns,1.*ns); //ibid.
    //these singlet and triplet times may not be the ones you're
    //used to, but are the world average: Kubota 79, Hitachi 83 (2
    //data sets), Teymourian 11, Morikawa 89, and Akimov '02
    G4double aSecondaryTime = vertex.gett();
    double SingTripRatioX, SingTripRatioR;
    if(vertex.getSpecies()==-1){
        //disregard tauR from original model--it's very small for any electric field.
        SingTripRatioX = G4RandGauss::shoot(0.17,0.05);
        SingTripRatioR = G4RandGauss::shoot(0.8, 0.2);
    }
    else if(vertex.getSpecies()==4){
        SingTripRatioR = G4RandGauss::shoot(2.3,0.51);
        SingTripRatioX = SingTripRatioR;
    }
    else{//NR
        SingTripRatioR = G4RandGauss::shoot(7.8,1.5);
        SingTripRatioX = SingTripRatioR;
    }
    if(exciton){
        if (G4UniformRand() < SingTripRatioR / (1 + SingTripRatioR))
            aSecondaryTime -= tau1 * log(G4UniformRand());
        else aSecondaryTime -= tau3 * log(G4UniformRand());
    } else {
        if (G4UniformRand() < SingTripRatioX / (1 + SingTripRatioX))
            aSecondaryTime -= tau1 * log(G4UniformRand());
        else aSecondaryTime -= tau3 * log(G4UniformRand());
    }
    if ( aSecondaryTime < 0 ) aSecondaryTime = 0;
    
    G4ThreeVector pos(vertex.getPos()[0],vertex.getPos()[1],vertex.getPos()[2]);
    return new G4Track(aQuantum,aSecondaryTime,pos);
    
}




G4VParticleChange*
NESTProc::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)
{
    aParticleChange.Initialize(aTrack);
    //ready to pop out OP and TE?
    if(NESTStackingAction::theStackingAction->isUrgentEmpty()
            && !unmerged_steps.empty()
            && aStep.GetSecondary()->empty()){
    
        const double Efield = 300;      
        std::vector<NEST::Vertex> merged_steps = NEST::cluster(unmerged_steps);
        for (auto vertex : merged_steps){       
            auto result = fNESTcalc->FullCalculation(vertex.getSpecies(),vertex.getEnergy(),vertex.getDensity(),Efield);
            
            for(auto time : result.photon_times){
                G4Track* onePhoton = MakePhoton(vertex,time);
                aParticleChange.AddSecondary(onePhoton);
            }
            for(int ie =0; ie<result.quanta.electrons; ++ie){
                G4ThreeVector pos(vertex.getPos()[0],vertex.getPos()[1],vertex.getPos()[2]);
//                Analysis::GetInstance()->AddThermalElectron(pos);
            }

        }
        unmerged_steps.clear();
    }
   return G4VRestDiscreteProcess::AtRestDoIt(aTrack, aStep); 
    
}

G4VParticleChange*
NESTProc::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
// this is the most important function, where all light & charge yields happen!
{
    aParticleChange.Initialize(aTrack);
    if (aStep.GetTotalEnergyDeposit() <= 0) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);



    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    G4ParticleDefinition *pDef = aParticle->GetDefinition();
    G4String particleName = pDef->GetParticleName();
    const G4Material* preMaterial = aStep.GetPreStepPoint()->GetMaterial();
    const G4Material* postMaterial = aStep.GetPostStepPoint()->GetMaterial();
    G4MaterialPropertiesTable* preMatTable = NULL;
    G4MaterialPropertiesTable* postMatTable = NULL;
    if (preMaterial) {
        preMatTable = preMaterial->GetMaterialPropertiesTable();
    }
    if (postMaterial) {
        postMatTable = postMaterial->GetMaterialPropertiesTable();
    }



    // code for determining whether the present/next material is noble
    // element, or, in other words, for checking if either is a valid NEST
    // scintillating material, and save Z for later L calculation, or
    // return if no valid scintillators are found on this step, which is
    // protection against G4Exception or seg. fault/violation
    G4Element *ElementA = NULL, *ElementB = NULL;
    if (preMaterial) {
        const G4ElementVector* theElementVector1 =
                preMaterial->GetElementVector();
        ElementA = (*theElementVector1)[0];
    }
    if (postMaterial) {
        const G4ElementVector* theElementVector2 =
                postMaterial->GetElementVector();
        ElementB = (*theElementVector2)[0];
    }
    G4int z1, z2;
    G4bool NobleNow = false, NobleLater = false;
    if (ElementA) z1 = (G4int) (ElementA->GetZ());
    else z1 = -1;
    if (ElementB) z2 = (G4int) (ElementB->GetZ());
    else z2 = -1;
    if (z1 == 2 || z1 == 10 || z1 == 18 || z1 == 36 || z1 == 54) {
        NobleNow = true;

    } //end of atomic number check
    if (z2 == 2 || z2 == 10 || z2 == 18 || z2 == 36 || z2 == 54) {
        NobleLater = true;

    } //end of atomic number check

    // Set step and production cut limits


    G4double max_step = .5 * um;
    G4LogicalVolume* scintvol = NULL;
    if (NobleNow) {
        scintvol = aStep.GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
        if (!scintvol->GetUserLimits() || scintvol->GetUserLimits()->GetMaxAllowedStep(aTrack) != max_step) {
            G4UserLimits* stepLimit = new G4UserLimits(max_step);
            scintvol->SetUserLimits(stepLimit);

        }
        assert(scintvol->GetMaterialCutsCouple()->GetProductionCuts()->GetProductionCut(electron_cut_index) <= 500 * nm);
    }
    if (NobleLater) {
        scintvol = aStep.GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
        if (!scintvol->GetUserLimits() || scintvol->GetUserLimits()->GetMaxAllowedStep(aTrack) != max_step) {
            G4UserLimits* stepLimit = new G4UserLimits(max_step);
            scintvol->SetUserLimits(stepLimit);

        }
        assert(scintvol->GetMaterialCutsCouple()->GetProductionCuts()->GetProductionCut(electron_cut_index) <= 500 * nm);
    }


          if ( !NobleLater )  
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);



    // retrieval of the particle's position, time, attributes at both the 
    // beginning and the end of the current step along its track
    G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    G4ThreeVector x1 = pPostStepPoint->GetPosition();
    G4ThreeVector x0 = pPreStepPoint->GetPosition();
    G4double evtStrt = pPreStepPoint->GetGlobalTime();
    G4double t0 = pPreStepPoint->GetLocalTime();
    G4double t1 = pPostStepPoint->GetLocalTime();


    G4double Density = preMaterial->GetDensity() / (g / cm3);
    if (fNESTcalc == NULL) {//initialize fNESTcacl
        // retrieve scintillation-related material properties
        
        G4double nDensity = Density*AVO; //molar mass factor applied below



        double molarmass = 0;
        const G4ElementVector* elvector = preMaterial->GetElementVector();
        const G4double* fvector = preMaterial->GetFractionVector();
        for (size_t i=0; i<preMaterial->GetNumberOfElements(); i++){
            volatile double Ai = (*elvector)[i]->GetA();
            molarmass += Ai*(fvector[i])/g;
        }
        //TODO: pass molar mass to NEST, or does it need it?
        fNESTcalc = std::unique_ptr<NEST::NESTcalc>(new NEST::NESTcalc());
    }


    //TODO:: calculate species from track
    INTERACTION_TYPE step_species = beta;
    
    double step_E = aStep.GetTotalEnergyDeposit()/keV;
    std::array<double, 3> pos{x0[0], x0[1], x0[2]};
    NEST::Vertex newvertex(step_E, step_species, pos, aStep.GetPostStepPoint()->GetGlobalTime() / ns,Density);
    unmerged_steps.push_back(newvertex);

    //the end (exiting)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------
G4double NESTProc::GetMeanFreePath(const G4Track&,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;
	// what this does is enforce the G4S1Light physics process as always
	// happening, so in effect scintillation is a meta-process on top of
	// any and all other energy depositions which may occur, just like the
	// original G4Scintillation (disregard DBL_MAX, this function makes the
	// mean free path zero really, not infinite)

        return DBL_MAX; //a C-defined constant
}

// GetMeanLifeTime
// ---------------
G4double NESTProc::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
        *condition = Forced;
	// this function and this condition has the same effect as the above
        return DBL_MAX;
}

