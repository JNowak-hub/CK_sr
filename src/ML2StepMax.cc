#include "ML2StepMax.hh"
#include "ML2StepMaxMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
ML2StepMax::ML2StepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),MaxChargedStep(DBL_MAX)
{
  pMess = new ML2StepMaxMessenger(this);
}
 
/////////////////////////////////////////////////////////////////////////////
ML2StepMax::~ML2StepMax()
{ delete pMess; }

/////////////////////////////////////////////////////////////////////////////
G4bool ML2StepMax::IsApplicable(const G4ParticleDefinition& particle) 
{ 
  return (particle.GetPDGCharge() != 0.);
}

/////////////////////////////////////////////////////////////////////////////    
void ML2StepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}

/////////////////////////////////////////////////////////////////////////////
G4double ML2StepMax::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double ProposedStep = DBL_MAX;

  if((MaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() == "DetectorPhys"))
     ProposedStep = MaxChargedStep;

  return ProposedStep;
}

/////////////////////////////////////////////////////////////////////////////
G4VParticleChange* ML2StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

