// $Id: MySession.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), june 2003

#include "G4UIsession.hh"
#include <fstream>
#include <iostream>

#ifndef MySession_H
#define MySession_H 1

using namespace std;

 class MySession : public G4UIsession
 {
 public:
 G4int ReceiveG4cout(G4String coutString);
 G4int ReceiveG4cerr(G4String cerrString);
 };

 G4int MySession::ReceiveG4cout(G4String coutString)
 {
 ofstream logFile;
logFile.open("results",ios::app);
 logFile << coutString << flush;
 return 0;
 }

 G4int MySession::ReceiveG4cerr(G4String cerrString)
{
 ofstream logFile1;
 logFile1.open("errors",ios::app);
 logFile1 << cerrString << flush;
 return 0;
 }
 
#endif