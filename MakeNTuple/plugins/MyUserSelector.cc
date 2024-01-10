// -*- C++ -*-
//
// Package:    UserCode/MyUserSelector
// Class:      MyUserSelector
// 
/**\class MyUserSelector MyUserSelector.cc UserCode/MyUserSelector/plugins/MyUserSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joao Pedro Gomes Pinheiro
//         Created:  Fri, 09 Apr 2021 16:59:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Math/Vector4D.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <TMath.h>
#include <TLorentzVector.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

//PAT
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

//ROOT
#include "TTree.h"
#include <TGraph.h>

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

using namespace std;
//using namespace reco;
using namespace edm;

//
// class declaration
//

class MyUserSelector : public edm::stream::EDFilter<> {
   public:
      explicit MyUserSelector(const edm::ParameterSet&);
      ~MyUserSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

         edm::Handle<pat::MuonCollection> muons;
         edm::EDGetTokenT<pat::MuonCollection> muonsToken;
         edm::Handle<pat::ElectronCollection> electrons;
         edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
         edm::Handle<reco::VertexCollection> vertices;
         edm::EDGetTokenT<reco::VertexCollection> verticesToken;

	double Z_vtx, pT_min, eta_max;

	// ----------------------------------- SELECTION VARIABLES -------------------------
	
	bool vtx_quality; // |z|<15.cm and vxtIsVald
	bool lepton_pair; // 2 lepton of opposite charge, pT > 38GeV, |eta|<2.4 and loose ID.
	bool select_event; // Final condition to select the event (vtx_quality && lepton_pair)

	// ---------------------------------------------------------------------------------

         //VERTEX INFO
         int nVtx;
         bool vtx_isValid, vtx_isFake;
         double vtx_z;
         double lepton_P_pT;
         double lepton_N_pT;
         double lepton_P_index;
         double lepton_N_index;
         double lepton_P_ismuon;
         double lepton_N_ismuon;
         double lepton_P_iselectron;
         double lepton_N_iselectron;
         double lepton_P_eta;
         double lepton_N_eta;
         double lepton_P_phi;
         double lepton_N_phi;
         double lepton_P_E;
         double lepton_N_E;
 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyUserSelector::MyUserSelector(const edm::ParameterSet& iConfig):
          muonsToken (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
        , electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
        , verticesToken (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
	, Z_vtx (iConfig.getParameter<double>("Z_vtx"))
	, pT_min(iConfig.getParameter<double>("pT_min"))
	, eta_max(iConfig.getParameter<double>("eta_max"))
{
   //now do what ever initialization is needed

}

MyUserSelector::~MyUserSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MyUserSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

        iEvent.getByToken(muonsToken, muons);
        iEvent.getByToken(electronsToken, electrons);
        iEvent.getByToken(verticesToken, vertices);

        nVtx = vertices->size();

        vtx_isValid = false;
        vtx_isFake = true;
        vtx_z = -999;

        if (nVtx>0){
                vtx_isValid = vertices->at(0).isValid();
                vtx_isFake = vertices->at(0).isFake();
                vtx_z = vertices->at(0).z();
        }
	lepton_P_pT = 0;
	lepton_N_pT = .0;
	lepton_P_index = 0.;
	lepton_N_index = 0.;
	lepton_P_ismuon = 0.;
	lepton_N_ismuon = 0.;
        lepton_P_iselectron = 0.;
        lepton_N_iselectron = 0.;
	lepton_P_eta = 0.;
	lepton_N_eta = 0.;
	lepton_P_phi = 0.;
        lepton_N_phi = 0.;
	lepton_P_E = 0.;
        lepton_N_E = 0.;
	for(size_t i = 0;i<electrons->size();i++){
	     if (electrons->at(i).pt()>pT_min
             && abs(electrons->at(i).eta())<eta_max
             && electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V2-loose")){
                if (electrons->at(i).charge() == 1){
			if (electrons->at(i).pt()>lepton_P_pT){
				lepton_P_pT = electrons->at(i).pt();
				lepton_P_index = i;
				lepton_P_ismuon = 0;
				lepton_P_iselectron = 1;
			}
		}
                if (electrons->at(i).charge() == -1){
                        if (electrons->at(i).pt()>lepton_N_pT){
                                lepton_N_pT = electrons->at(i).pt();
                                lepton_N_index = i;
                                lepton_N_ismuon = 0;
                                lepton_N_iselectron = 1;
                        }
                }
	    }
	}
        for(size_t i = 0;i<muons->size();i++){
           if(muons->at(i).pt()>pT_min &&
             abs(muons->at(i).eta())<eta_max &&
             muon::isLooseMuon(muons->at(i))){
                if (muons->at(i).charge() == 1){
                        if (muons->at(i).pt()>lepton_P_pT){
                                lepton_P_pT = muons->at(i).pt();
                                lepton_P_index = i;
                                lepton_P_ismuon = 1;
                                lepton_P_iselectron = 0;
                        }
                }
                if (muons->at(i).charge() == -1){
                        if (muons->at(i).pt()>lepton_N_pT){
                                lepton_N_pT = muons->at(i).pt();
                                lepton_N_index = i;
                                lepton_N_ismuon = 1;
                                lepton_N_iselectron = 0;
                        }
                }
            }
        }
	//std::cout << "lepton_P_pT: " << lepton_P_pT << std::endl;
	//std::cout << "lepton_N_pT: " << lepton_N_pT << std::endl;
	if (lepton_N_ismuon==1){
		//std::cout << "muons->at(lepton_N_index).pt(): " << muons->at(lepton_N_index).pt() << std::endl;
		lepton_N_pT  = muons->at(lepton_N_index).pt();
		lepton_N_phi = muons->at(lepton_N_index).phi();		
		lepton_N_eta = muons->at(lepton_N_index).eta();
                lepton_N_E   = muons->at(lepton_N_index).energy();
	}
	if (lepton_N_iselectron==1){
		//std::cout << "electrons->at(lepton_N_index).pt(): " << electrons->at(lepton_N_index).pt() << std::endl;
		lepton_N_pT  = electrons->at(lepton_N_index).pt();
                lepton_N_phi = electrons->at(lepton_N_index).phi();
                lepton_N_eta = electrons->at(lepton_N_index).eta();
                lepton_N_E   = electrons->at(lepton_N_index).energy();
	}
	if (lepton_P_ismuon==1){
		//std::cout << "muons->at(lepton_P_index).pt(): " << muons->at(lepton_P_index).pt() << std::endl;
                lepton_P_pT  = muons->at(lepton_P_index).pt();
                lepton_P_phi = muons->at(lepton_P_index).phi();
                lepton_P_eta = muons->at(lepton_P_index).eta();
                lepton_P_E   = muons->at(lepton_P_index).energy();
	}
        if (lepton_P_iselectron==1){
		//std::cout << "electrons->at(lepton_P_index).pt(): " << electrons->at(lepton_P_index).pt() << std::endl;
                lepton_P_pT  = electrons->at(lepton_P_index).pt();
                lepton_P_phi = electrons->at(lepton_P_index).phi();
                lepton_P_eta = electrons->at(lepton_P_index).eta();
                lepton_P_E   = electrons->at(lepton_P_index).energy();
	}
	ROOT::Math::PtEtaPhiEVector lepton_P(lepton_P_pT,lepton_P_eta,lepton_P_phi,lepton_P_E);
        ROOT::Math::PtEtaPhiEVector lepton_N(lepton_N_pT,lepton_N_eta,lepton_N_phi,lepton_N_E);	

	double mll;
	
	mll = (lepton_P + lepton_N).mass();
        if (mll>110.) {std::cout << "mll = " << mll << std::endl;}
		
	vtx_quality=false;
	if (vtx_isValid && abs(vtx_z)<Z_vtx) {vtx_quality=true;} // FIRST pre-SELECTION CRITERIA	

	lepton_pair=false;
	if (mll>110.){lepton_pair=true;}

	select_event = false;
	select_event = vtx_quality && lepton_pair;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return select_event;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MyUserSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MyUserSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyUserSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MyUserSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyUserSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyUserSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyUserSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MyUserSelector);
