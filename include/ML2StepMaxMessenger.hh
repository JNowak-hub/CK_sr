#ifndef ML2StepMaxMessenger_h
#define ML2StepMaxMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class ML2StepMax;
class G4UIcmdWithADoubleAndUnit;

/////////////////////////////////////////////////////////////////////////////
class ML2StepMaxMessenger: public G4UImessenger
{
  public:
    ML2StepMaxMessenger(ML2StepMax*);
   ~ML2StepMaxMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ML2StepMax* stepMax;
    G4UIcmdWithADoubleAndUnit* StepMaxCmd;
};

#endif
