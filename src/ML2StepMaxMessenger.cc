#include "ML2StepMaxMessenger.hh"
#include "ML2StepMax.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

/////////////////////////////////////////////////////////////////////////////
ML2StepMaxMessenger::ML2StepMaxMessenger(ML2StepMax* stepM)
:stepMax(stepM)
{ 
  StepMaxCmd = new G4UIcmdWithADoubleAndUnit("/Step/waterPhantomStepMax",this);
  StepMaxCmd->SetGuidance("Set max allowed step length");
  StepMaxCmd->SetParameterName("mxStep",false);
  StepMaxCmd->SetRange("mxStep>0.");
  StepMaxCmd->SetUnitCategory("Length");
}

/////////////////////////////////////////////////////////////////////////////
ML2StepMaxMessenger::~ML2StepMaxMessenger()
{
  delete StepMaxCmd;
}

/////////////////////////////////////////////////////////////////////////////
void ML2StepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == StepMaxCmd)
    { stepMax->SetMaxStep(StepMaxCmd->GetNewDoubleValue(newValue));}
}

