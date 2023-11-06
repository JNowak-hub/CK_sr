//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"

// Geant4-DNA MODELS

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNABornExcitationModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNARuddIonisationModel.hh"

#include "G4DNAChargeDecrease.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"

#include "G4DNAChargeIncrease.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAMeltonAttachmentModel.hh"

#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"

#include "G4DNAElectronSolvation.hh"

//
#include "QGSP_BIC_HP.hh"
#include "G4PreCompoundModel.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VEmModel.hh"
#include "G4DummyModel.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4MesonConstructor.hh"

#include "G4ElectronCapture.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronicParameters.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4LundStringFragmentation.hh"
#include "G4CrossSectionElastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ChipsElasticModel.hh"

// seb
#include "G4PhysicsListHelper.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4DNAElastic.hh"
#include "G4DNAGenericIonsManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VUserPhysicsList()
{
  defaultCutValue = 1 * micrometer;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  // Construct all mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructBosons()
{
  G4Gamma::GammaDefinition();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructLeptons()
{
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructBarions()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4GenericIon::GenericIonDefinition();

  // Geant4 DNA new particles
  G4DNAGenericIonsManager *genericIonsManager;
  genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while ((*particleIterator)())
  {

    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // *********************************
    // 1) Processes for the World region
    // *********************************

    if (particleName == "e-")
    {

      // STANDARD msc is active in the world
      G4eMultipleScattering *msc = new G4eMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4eIonisation *eion = new G4eIonisation();
      eion->SetEmModel(new G4MollerBhabhaModel());
      ph->RegisterProcess(eion, particle);

      // DNA elastic is not active in the world
      G4DNAElastic *theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("e-_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("e-_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA attachment is not active in the world
      G4DNAAttachment *dnaatt = new G4DNAAttachment("e-_G4DNAAttachment");
      dnaatt->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaatt);

      // DNA vib. excitation is not active in the world
      G4DNAVibExcitation *dnavib =
          new G4DNAVibExcitation("e-_G4DNAVibExcitation");
      dnavib->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnavib);

      // THE FOLLOWING PROCESS WILL KILL ALL ELECTRONS BELOW A
      // SELECTED ENERY THRESHOLD
      // Capture of low-energy e-
      G4ElectronCapture *ecap = new G4ElectronCapture("Target", 7.4 * eV);
      pmanager->AddDiscreteProcess(ecap);
      // 7.4 eV is compatible with the validity range of the elastic model
    }
    else if (particleName == "proton")
    {
      G4ChipsElasticModel *elastic_chip = new G4ChipsElasticModel();
      // Inelastic scattering
      const G4double theFTFMin0 = 0.0 * GeV;
      const G4double theFTFMin1 = 3.0 * GeV;
      const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
      const G4double theBERTMin0 = 0.0 * GeV;
      const G4double theBERTMin1 = 19.0 * MeV;
      const G4double theBERTMax = 6.0 * GeV;
      const G4double theHPMin = 0.0 * GeV;
      const G4double theHPMax = 20.0 * MeV;

      G4FTFModel *theStringModel = new G4FTFModel;
      G4ExcitedStringDecay *theStringDecay = new G4ExcitedStringDecay(new G4LundStringFragmentation);
      theStringModel->SetFragmentationModel(theStringDecay);
      G4PreCompoundModel *thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
      G4GeneratorPrecompoundInterface *theCascade = new G4GeneratorPrecompoundInterface(thePreEquilib);

      G4TheoFSGenerator *theFTFModel0 = new G4TheoFSGenerator("FTFP");
      theFTFModel0->SetHighEnergyGenerator(theStringModel);
      theFTFModel0->SetTransport(theCascade);
      theFTFModel0->SetMinEnergy(theFTFMin0);
      theFTFModel0->SetMaxEnergy(theFTFMax);

      G4TheoFSGenerator *theFTFModel1 = new G4TheoFSGenerator("FTFP");
      theFTFModel1->SetHighEnergyGenerator(theStringModel);
      theFTFModel1->SetTransport(theCascade);
      theFTFModel1->SetMinEnergy(theFTFMin1);
      theFTFModel1->SetMaxEnergy(theFTFMax);

      G4CascadeInterface *theBERTModel0 = new G4CascadeInterface;
      theBERTModel0->SetMinEnergy(theBERTMin0);
      theBERTModel0->SetMaxEnergy(theBERTMax);

      G4CascadeInterface *theBERTModel1 = new G4CascadeInterface;
      theBERTModel1->SetMinEnergy(theBERTMin1);
      theBERTModel1->SetMaxEnergy(theBERTMax);

      G4VCrossSectionDataSet *theAntiNucleonData = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS);
      G4ComponentGGNuclNuclXsc *ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
      G4VCrossSectionDataSet *theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);
      G4VCrossSectionDataSet *theGGNNEl = new G4CrossSectionElastic(ggNuclNuclXsec);
      G4ComponentGGHadronNucleusXsc *ggHNXsec = new G4ComponentGGHadronNucleusXsc();
      G4VCrossSectionDataSet *theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
      G4VCrossSectionDataSet *theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);
      // Elastic scattering
      G4HadronElasticProcess *theElasticProcess = new G4HadronElasticProcess;
      theElasticProcess->AddDataSet(new G4BGGNucleonElasticXS(G4Proton::Proton()));
      theElasticProcess->RegisterMe(elastic_chip);
      pmanager->AddDiscreteProcess(theElasticProcess);
      // Inelastic scattering
      G4HadronInelasticProcess *theInelasticProcess =
          new G4HadronInelasticProcess("inelastic", G4Proton::Definition());
      theInelasticProcess->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
      theInelasticProcess->RegisterMe(theFTFModel1);
      theInelasticProcess->RegisterMe(theBERTModel0);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4WentzelVIModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD Coulomb scattering is active in the world
      G4CoulombScattering *pcou = new G4CoulombScattering();
      ph->RegisterProcess(pcou, particle);

      // STANDARD ionisation is active in the world
      G4hIonisation *hion = new G4hIonisation();
      hion->SetEmModel(new G4BraggModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("proton_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("proton_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("proton_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);
    }
    else if (particleName == "hydrogen")
    {

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("hydrogen_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("hydrogen_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("hydrogen_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
    else if (particleName == "GenericIon")
    {

      // WARNING : THIS IS NEEDED FOR STANDARD ALPHA G4ionIonisation PROCESS

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);
    }
    else if (particleName == "alpha")
    {

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("alpha_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("alpha_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("alpha_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);
    }
    else if (particleName == "alpha+")
    {

      // STANDARD msc is active in the world
      G4hMultipleScattering *msc = new G4hMultipleScattering();
      msc->SetEmModel(new G4UrbanMscModel());
      ph->RegisterProcess(msc, particle);

      // STANDARD ionisation is active in the world
      G4ionIonisation *hion = new G4ionIonisation();
      hion->SetEmModel(new G4BraggIonModel());
      hion->SetEmModel(new G4BetheBlochModel());
      ph->RegisterProcess(hion, particle);

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("alpha+_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("alpha+_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("alpha+_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge decrease is not active in the world
      G4DNAChargeDecrease *dnacd = new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
      dnacd->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnacd);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
    else if (particleName == "helium")
    {

      // DNA excitation is not active in the world
      G4DNAExcitation *dnaex = new G4DNAExcitation("helium_G4DNAExcitation");
      dnaex->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaex);

      // DNA ionisation is not active in the world
      G4DNAIonisation *dnaioni = new G4DNAIonisation("helium_G4DNAIonisation");
      dnaioni->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaioni);

      // DNA elastic is not active in the world
      G4DNAElastic *dnael = new G4DNAElastic("helium_G4DNAElastic");
      dnael->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnael);

      // DNA charge increase is not active in the world
      G4DNAChargeIncrease *dnaci = new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
      dnaci->SetEmModel(new G4DummyModel());
      pmanager->AddDiscreteProcess(dnaci);
    }
  }
  // **************************************
  // 2) Define processes for Target region
  // **************************************

  // STANDARD EM processes should be inactivated when
  // corresponding DNA processes are used
  // - STANDARD EM e- processes are inactivated below 1 MeV
  // - STANDARD EM proton & alpha processes are inactivated below
  //   standEnergyLimit
  G4double standEnergyLimit = 9.9 * MeV;
  //

  G4double massFactor = 1.0079 / 4.0026;
  G4EmConfigurator *em_config =
      G4LossTableManager::Instance()->EmConfigurator();

  G4VEmModel *mod;

  // *** e-

  // ---> STANDARD EM processes are inactivated below 1 MeV

  mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-", "msc", mod, "Target");

  mod = new G4MollerBhabhaModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("e-",
                             "eIoni",
                             mod,
                             "Target",
                             0.0,
                             100 * TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAChampionElasticModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
                             mod, "Target", 7.4 * eV, 1. * MeV);

  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
                             mod, "Target", 11. * eV, 1. * MeV);
  // Note: valid from 11 eV to 0.999.. MeV then switch to std models at
  // higher energies ; same for other models

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
                             mod, "Target", 9. * eV, 1. * MeV);

  mod = new G4DNAMeltonAttachmentModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
                             mod, "Target", 4. * eV, 13. * eV);

  mod = new G4DNASancheExcitationModel();
  em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
                             mod, "Target", 2. * eV, 100. * eV);

  // *** proton

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4WentzelVIModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("proton", "msc", mod, "Target", 0, 100 * TeV);

  mod = new G4eCoulombScatteringModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("proton", "CoulombScat", mod, "Target", 0, 100 * TeV);

  mod = new G4BraggModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("proton", "hIoni",
                             mod, "Target", 0.0, 2 * MeV, new G4IonFluctuations());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("proton", "hIoni",
                             mod, "Target", 2 * MeV, 100 * TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
                             mod, "Target", 100 * eV, 0.5 * MeV);

  mod = new G4DNABornIonisationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
                             mod, "Target", 0.5 * MeV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
                             mod, "Target", 10 * eV, 0.5 * MeV);

  mod = new G4DNABornExcitationModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
                             mod, "Target", 0.5 * MeV, 10 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("proton", "proton_G4DNAChargeDecrease",
                             mod, "Target", 100 * eV, 10 * MeV);

  // *** hydrogen

  // ---> NO STANDARD EM processes

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAElastic",
                             mod, "Target", 100. * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAIonisation",
                             mod, "Target", 100 * eV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAExcitation",
                             mod, "Target", 10 * eV, 0.5 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAChargeIncrease",
                             mod, "Target", 100 * eV, 10 * MeV);

  // *** alpha

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("alpha", "msc", mod, "Target", 0, 100 * TeV);

  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha", "ionIoni",
                             mod, "Target", 0.0, 2 * MeV / massFactor,
                             new G4IonFluctuations());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha", "ionIoni",
                             mod, "Target", 2 * MeV / massFactor, 100 * TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha", "alpha_G4DNAChargeDecrease",
                             mod, "Target", 1 * keV, 10 * MeV);

  // *** alpha+

  // ---> STANDARD EM processes inactivated below standEnergyLimit
  //      or below 1 MeV for scattering

  // Inactivate following STANDARD processes

  mod = new G4UrbanMscModel();
  mod->SetActivationLowEnergyLimit(1. * MeV);
  em_config->SetExtraEmModel("alpha+", "msc", mod, "Target", 0, 100 * TeV);

  mod = new G4BraggIonModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
                             mod, "Target", 0.0, 2 * MeV / massFactor,
                             new G4IonFluctuations());

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(standEnergyLimit);
  em_config->SetExtraEmModel("alpha+", "ionIoni",
                             mod, "Target", 2 * MeV / massFactor, 100 * TeV,
                             new G4UniversalFluctuation());

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeIncrease",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeDecreaseModel();
  em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeDecrease",
                             mod, "Target", 1 * keV, 10 * MeV);

  // *** helium

  // ---> NO STANDARD EM processes

  // ---> DNA processes activated

  mod = new G4DNAIonElasticModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAElastic",
                             mod, "Target", 100 * eV, 1. * MeV);

  mod = new G4DNARuddIonisationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAIonisation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNAMillerGreenExcitationModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAExcitation",
                             mod, "Target", 1 * keV, 10 * MeV);

  mod = new G4DNADingfelderChargeIncreaseModel();
  em_config->SetExtraEmModel("helium", "helium_G4DNAChargeIncrease",
                             mod, "Target", 1 * keV, 10 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructGeneral()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
