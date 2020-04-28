#ifndef LORASteppingAction_h
#define LORASteppingAction_h 1

#include "G4UserSteppingAction.hh"

class LORARunAction;
class LORAEventAction;


class LORASteppingAction : public G4UserSteppingAction
{
public:
   LORASteppingAction(LORARunAction*, LORAEventAction*);
   ~LORASteppingAction() {};
   
   void UserSteppingAction(const G4Step*);
   
private:
   LORARunAction*    runAction;
   LORAEventAction*  eventAction;
};


#endif
