// -*- C++ -*-
//
// Package:    MakeNTuple/MakeNTuple
// Class:      MakeNTuple
// 
/**\class MakeNTuple MakeNTuple.cc MakeNTuple/MakeNTuple/plugins/MakeNTuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Tue, 09 Apr 2019 16:31:41 GMT
//
// Modified by JOAO PEDRO (since 17/02/21)
// The backup of this analyzer is at /afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/SemiLep_Studies/backup


// system include files
#include <memory>
#include <iostream> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TRandom.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h" 
#include <vector>
#include <TMath.h>
#include <TH2F.h>
#include <TString.h>
//HLT
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


//PAT
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//ROOT
#include "TTree.h" 
#include <TGraph.h>

// P P S
/*
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHit.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHitContainer.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrack.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
*/
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include "RoccoR.h"
#include "TRandom3.h"

using namespace std; 
using namespace reco;
using namespace edm;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MakeNTuple : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit MakeNTuple(const edm::ParameterSet&);
		~MakeNTuple();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::Handle<pat::MuonCollection> muons;
		edm::EDGetTokenT<pat::MuonCollection> muonsToken;
		edm::Handle<pat::ElectronCollection> electrons;
		edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
		edm::Handle<pat::METCollection> MET;
		edm::EDGetTokenT<pat::METCollection> MetToken;
		edm::Handle<std::vector<pat::PackedCandidate>> PFCand;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> PFCandToken;
		edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
		edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSumInfoInputTag;
		edm::Handle<reco::VertexCollection> vertices;
		edm::EDGetTokenT<reco::VertexCollection> verticesToken;
		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;

		edm::LumiReWeighting *LumiWeights_;

		bool MC, MC_Signal, DATA;
		std::string tag;
		bool preVFP, postVFP;
		edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

		// FOR GEN INFO
		edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken;
		edm::Handle< std::vector<reco::GenParticle> > genHandle;
		edm::EDGetTokenT< edm::HepMCProduct > mcEventToken;
		edm::Handle< edm::HepMCProduct > EvtHandle ;
		void Get_t_and_xi(const TLorentzVector* p,double& t, double& xi);
		double GetWeightLeptons(double pt, double eta, TH2F* hist2D);
		const double ProtonMass = 0.93827;
		const double ProtonMassSQ = pow(ProtonMass,2);
		double fBeamEnergy;
		double fBeamMomentum;
		void set_BeamEnergy(double e) {fBeamEnergy=e;fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);};

		//FOR PU REWGT

		TTree* EventBranchs;
		int BX, Run, LumiSection, EventNum, xangle;
		double PUWeight, nPU;
		//VERTEX INFO
		int nVtx;
		bool vtx_isValid, vtx_isFake;
		double vtx_z;

		HLTConfigProvider hltConfig_;
		HLTPrescaleProvider hltPrescaleProvider_;

		std::vector<double> *ArmF_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thy_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thy_gen        = new std::vector<double> ();

		std::vector<double> *ProtCand_xi		= new std::vector<double> ();
		std::vector<double> *ProtCand_t			= new std::vector<double> ();
		std::vector<double> *ProtCand_ThX		= new std::vector<double> ();
		std::vector<double> *ProtCand_ThY		= new std::vector<double> ();
		std::vector<double> *ProtCand_rpid		= new std::vector<double> ();
		std::vector<double> *ProtCand_arm       	= new std::vector<double> ();
		std::vector<double> *ProtCand_ismultirp		= new std::vector<double> ();
                std::vector<double> *ProtCand_x                 = new std::vector<double> ();
                std::vector<double> *ProtCand_y                 = new std::vector<double> ();
                std::vector<double> *ProtCand_xn                = new std::vector<double> ();
                std::vector<double> *ProtCand_yn                = new std::vector<double> ();
                std::vector<double> *ProtCand_xf                = new std::vector<double> ();
                std::vector<double> *ProtCand_yf                = new std::vector<double> ();
		double ProtCand_weight;

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
		double lepton_N_px;
		double lepton_P_px;
                double lepton_N_py;
                double lepton_P_py;
                double lepton_N_pz;
                double lepton_P_pz;
                double lepton_N_vtxZ;
                double lepton_P_vtxZ;
                double lepton_N_charge;
                double lepton_P_charge;
                double lepton_N_PFBasedIso;
                double lepton_P_PFBasedIso;
                double lepton_N_isTight;
                double lepton_P_isTight;
                double lepton_N_isMedium;
                double lepton_P_isMedium;

		double prefweight;
		double prefweightup;
		double prefweightdown;

		int channel; // this variable will be set as:
			     // 0 if the final state leptons are both muons
			     // 1 if the final state leptons are both electrons
			     // 2 if the final state leptons are a muon and a electron

		double METPx, METPy, METPt, METphi;

		std::vector<double> *pfphi      = new std::vector<double> ();
		std::vector<double> *pfeta      = new std::vector<double> ();

		std::vector<int> *pffromPV      = new std::vector<int> ();
		std::vector<double> *pfdz      = new std::vector<double> ();
		std::vector<double> *pfpt      = new std::vector<double> ();

		std::vector<std::string> *HLT_name        = new std::vector<std::string> ();
		std::vector<bool> *HLT_pass        = new std::vector<bool> ();
		std::vector<int> *HLT_prescale        = new std::vector<int> ();
		std::vector<std::string> HLT_list = {"HLT_IsoMu24_v", "HLT_Ele27_WPTight_Gsf_v"};

		double lepton_N_gen_px;
                double lepton_N_gen_py;
                double lepton_N_gen_pz;
                double lepton_N_gen_pT;
                double lepton_N_gen_E;
                double lepton_N_gen_phi;
                double lepton_N_gen_eta;
                double lepton_P_gen_px;
                double lepton_P_gen_py;
                double lepton_P_gen_pz;
                double lepton_P_gen_pT;
                double lepton_P_gen_E;
                double lepton_P_gen_phi;
                double lepton_P_gen_eta;
		double lepton_P_gen_ismuon;
                double lepton_P_gen_iselectron;
                double lepton_N_gen_ismuon;
                double lepton_N_gen_iselectron;

		double neut_P_gen_px;
                double neut_P_gen_py;
                double neut_P_gen_pz;
                double neut_P_gen_pT;
                double neut_P_gen_E;
                double neut_P_gen_phi;
                double neut_P_gen_eta;
                double neut_N_gen_px;
                double neut_N_gen_py;
                double neut_N_gen_pz;
                double neut_N_gen_pT;
                double neut_N_gen_E;
                double neut_N_gen_phi;
                double neut_N_gen_eta;

		double Wp_gen_px;
                double Wp_gen_py;
                double Wp_gen_pz;
                double Wp_gen_pT;
                double Wp_gen_E;
                double Wp_gen_phi;
                double Wp_gen_eta;
                double Wm_gen_px;
                double Wm_gen_py;
                double Wm_gen_pz;
                double Wm_gen_pT;
                double Wm_gen_E;
                double Wm_gen_phi;
                double Wm_gen_eta;

                double muon_pT;
                double electron_pT;
                double muon_eta;
                double electron_eta;

		double weight1_id;
                double weight1_iso;
		double weight2_id;
		double weight2_iso;
		double muSF;
		double elSF;
		double weight_lep;
		double weight;

		double weight_multitrack_45_F_B;
		double weight_multitrack_45_F_C;
		double weight_multitrack_45_F_G;
		double weight_multitrack_45_F_H;
		double weight_multitrack_45_N_B;
		double weight_multitrack_45_N_C;
		double weight_multitrack_45_N_G;
		double weight_multitrack_45_N_H;
		double weight_multitrack_56_F_B;
		double weight_multitrack_56_F_C;
		double weight_multitrack_56_F_G;
		double weight_multitrack_56_F_H;
		double weight_multitrack_56_N_B;
		double weight_multitrack_56_N_C;
		double weight_multitrack_56_N_G;
		double weight_multitrack_56_N_H;

		double w_raddam_45_Multi_B_N;
		double w_raddam_45_Multi_C_N;
		double w_raddam_45_Multi_G_N;
		double w_raddam_45_Multi_B_F;
		double w_raddam_45_Multi_C_F;
		double w_raddam_45_Multi_G_F;
		double w_raddam_45_Multi;
		
		double weight_45_N;
		double weight_45_F;
		double weight_45;

		double w_raddam_56_Multi_B_N;
		double w_raddam_56_Multi_C_N;
		double w_raddam_56_Multi_G_N;
		double w_raddam_56_Multi_B_F;
		double w_raddam_56_Multi_C_F;
		double w_raddam_56_Multi_G_F;
		double w_raddam_56_Multi;

		double weight_56_N;
		double weight_56_F;
		double weight_56;

		double multiRP_weight_B;
		double multiRP_weight_C1;
		double multiRP_weight_C2;
		double multiRP_weight_D;
		double multiRP_weight_E;
		double multiRP_weight_F1;
		double multiRP_weight_F2;
		double multiRP_weight_F3;
		double multiRP_weight;

		double strips_weight_rddam_B;
		double strips_weight_rddam_C;
		double strips_weight_rddam_D;
		double strips_weight_rddam_E;
		double strips_weight_rddam_F;
		double strips_weight_rddam;

		double strips_weight_multitrack_B;
		double strips_weight_multitrack_C;
		double strips_weight_multitrack_D;
		double strips_weight_multitrack_E;
		double strips_weight_multitrack_F;
		double strips_weight_multitrack;
		
	
		double multiRP2018_weight_A;
		double multiRP2018_weight_B;
		double multiRP2018_weight_B1;
		double multiRP2018_weight_B2;
		double multiRP2018_weight_C;
		double multiRP2018_weight_D;
		double multiRP2018_weight_D1;
		double multiRP2018_weight_D2;
		
		double radiation2018_weight_A_N;
		double radiation2018_weight_B_N;
		double radiation2018_weight_B1_N;
		double radiation2018_weight_B2_N;
		double radiation2018_weight_C_N;
		double radiation2018_weight_D_N;
		double radiation2018_weight_D1_N;
		double radiation2018_weight_D2_N;
		
		double radiation2018_weight_A_F;
		double radiation2018_weight_B_F;
		double radiation2018_weight_B1_F;
		double radiation2018_weight_B2_F;
		double radiation2018_weight_C_F;
		double radiation2018_weight_D_F;
		double radiation2018_weight_D1_F;
		double radiation2018_weight_D2_F;
		double multiRP2018_weight;
		
		double weight2018_rddam_N;
		double weight2018_rddam_F;
		double weight2018_rddam;

		TFile *f_SFmuon_ID_preVFP   = new TFile("Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
	        TH2F  *h_SFmuon_ID_preVFP   = 0;

        	TFile *f_SFmuon_ISO_preVFP  = new TFile("Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root");
        	TH2F  *h_SFmuon_ISO_preVFP  = 0;

        	TFile *f_SFmuon_ID_postVFP  = new TFile("Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
        	TH2F  *h_SFmuon_ID_postVFP  = 0;

        	TFile *f_SFmuon_ISO_postVFP = new TFile("Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root");
        	TH2F  *h_SFmuon_ISO_postVFP = 0;

        	TFile *f_SFmuon_ID_2017     = new TFile("Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root");
        	TH2F  *h_SFmuon_ID_2017     = 0;
        	
		TFile *f_SFmuon_ISO_2017    = new TFile("Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root");
        	TH2F  *h_SFmuon_ISO_2017    = 0;

        	TFile *f_SFmuon_ID_2018     = new TFile("Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
        	TH2F  *h_SFmuon_ID_2018     = 0;

        	TFile *f_SFmuon_ISO_2018    = new TFile("Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
        	TH2F  *h_SFmuon_ISO_2018    = 0;

        	TFile *f_SFelectron_ID_preVFP   = new TFile("egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root");
        	TH2F  *h_SFelectron_ID_preVFP   = 0;

        	TFile *f_SFelectron_ID_postVFP  = new TFile("egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root");
        	TH2F  *h_SFelectron_ID_postVFP  = 0;

        	TFile *f_SFelectron_ID_2017     = new TFile("egammaEffi.txt_EGM2D_Tight_UL17.root");
        	TH2F  *h_SFelectron_ID_2017     = 0;

        	TFile *f_SFelectron_ID_2018     = new TFile("egammaEffi.txt_Ele_Tight_EGM2D.root");
        	TH2F  *h_SFelectron_ID_2018     = 0;

		TFile *eff_radiation = new TFile("pixelEfficiencies_radiation_reMiniAOD.root");
		TH2F *eff_radiation_A_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018A/h45_210_2018A_all_2D");
		TH2F *eff_radiation_A_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018A/h45_220_2018A_all_2D");
        	TH2F *eff_radiation_B1_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B1/h45_210_2018B1_all_2D");
        	TH2F *eff_radiation_B1_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B1/h45_220_2018B1_all_2D");
        	TH2F *eff_radiation_B2_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B2/h45_210_2018B2_all_2D");
        	TH2F *eff_radiation_B2_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B2/h45_220_2018B2_all_2D");
        	TH2F *eff_radiation_B_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B/h45_210_2018B_all_2D");
        	TH2F *eff_radiation_B_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B/h45_220_2018B_all_2D");
		TH2F *eff_radiation_C_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018C/h45_210_2018C_all_2D");
        	TH2F *eff_radiation_C_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018C/h45_220_2018C_all_2D");
        	TH2F *eff_radiation_D1_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D1/h45_210_2018D1_all_2D");
        	TH2F *eff_radiation_D1_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D1/h45_220_2018D1_all_2D");
        	TH2F *eff_radiation_D2_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D2/h45_210_2018D2_all_2D");
        	TH2F *eff_radiation_D2_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D2/h45_220_2018D2_all_2D");
        	TH2F *eff_radiation_D_45_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D/h45_210_2018D_all_2D");
        	TH2F *eff_radiation_D_45_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D/h45_220_2018D_all_2D");
        	TH2F *eff_radiation_A_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018A/h56_210_2018A_all_2D");
        	TH2F *eff_radiation_A_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018A/h56_220_2018A_all_2D");
        	TH2F *eff_radiation_B1_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B1/h56_210_2018B1_all_2D");
        	TH2F *eff_radiation_B1_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B1/h56_220_2018B1_all_2D");
        	TH2F *eff_radiation_B2_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B2/h56_210_2018B2_all_2D");
        	TH2F *eff_radiation_B2_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B2/h56_220_2018B2_all_2D");
        	TH2F *eff_radiation_B_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018B/h56_210_2018B_all_2D");
        	TH2F *eff_radiation_B_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018B/h56_220_2018B_all_2D");
        	TH2F *eff_radiation_C_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018C/h56_210_2018C_all_2D");
        	TH2F *eff_radiation_C_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018C/h56_220_2018C_all_2D");
        	TH2F *eff_radiation_D1_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D1/h56_210_2018D1_all_2D");
        	TH2F *eff_radiation_D1_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D1/h56_220_2018D1_all_2D");
        	TH2F *eff_radiation_D2_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D2/h56_210_2018D2_all_2D");
        	TH2F *eff_radiation_D2_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D2/h56_220_2018D2_all_2D");
        	TH2F *eff_radiation_D_56_N = (TH2F*)eff_radiation->Get("Pixel/2018/2018D/h56_210_2018D_all_2D");
        	TH2F *eff_radiation_D_56_F = (TH2F*)eff_radiation->Get("Pixel/2018/2018D/h56_220_2018D_all_2D");

        	TFile *eff_multiRP = new TFile("pixelEfficiencies_multiRP_reMiniAOD.root");
		TH2F *eff_multiRP_A_56_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018A/h56_220_2018A_all_2D");
		TH2F *eff_multiRP_B_56_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B/h56_220_2018B_all_2D");
		TH2F *eff_multiRP_B1_56_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B1/h56_220_2018B1_all_2D");
		TH2F *eff_multiRP_B2_56_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B2/h56_220_2018B2_all_2D");
		TH2F *eff_multiRP_C_56_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018C/h56_220_2018C_all_2D");
        	TH2F *eff_multiRP_D_56_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D/h56_220_2018D_all_2D");
		TH2F *eff_multiRP_D1_56_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D1/h56_220_2018D1_all_2D");
		TH2F *eff_multiRP_D2_56_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D2/h56_220_2018D2_all_2D");
        	TH2F *eff_multiRP_A_45_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018A/h45_220_2018A_all_2D");
        	TH2F *eff_multiRP_B_45_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B/h45_220_2018B_all_2D");
        	TH2F *eff_multiRP_B1_45_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B1/h45_220_2018B1_all_2D");
        	TH2F *eff_multiRP_B2_45_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018B2/h45_220_2018B2_all_2D");
        	TH2F *eff_multiRP_C_45_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018C/h45_220_2018C_all_2D");
        	TH2F *eff_multiRP_D_45_2018  = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D/h45_220_2018D_all_2D");
        	TH2F *eff_multiRP_D1_45_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D1/h45_220_2018D1_all_2D");
        	TH2F *eff_multiRP_D2_45_2018 = (TH2F*)eff_multiRP->Get("Pixel/2018/2018D2/h45_220_2018D2_all_2D");

        	TH2F *eff_multiRP_B_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017B/h45_220_2017B_all_2D");
        	TH2F *eff_multiRP_C1_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017C1/h45_220_2017C1_all_2D");
        	TH2F *eff_multiRP_C2_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017C2/h45_220_2017C2_all_2D");
        	TH2F *eff_multiRP_D_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017D/h45_220_2017D_all_2D");
        	TH2F *eff_multiRP_E_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017E/h45_220_2017E_all_2D");
        	TH2F *eff_multiRP_F1_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F1/h45_220_2017F1_all_2D");
        	TH2F *eff_multiRP_F2_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F2/h45_220_2017F2_all_2D");
        	TH2F *eff_multiRP_F3_45 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F3/h45_220_2017F3_all_2D");
        	TH2F *eff_multiRP_B_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017B/h56_220_2017B_all_2D");
        	TH2F *eff_multiRP_C1_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017C1/h56_220_2017C1_all_2D");
        	TH2F *eff_multiRP_C2_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017C2/h56_220_2017C2_all_2D");
        	TH2F *eff_multiRP_D_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017D/h56_220_2017D_all_2D");
        	TH2F *eff_multiRP_E_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017E/h56_220_2017E_all_2D");
        	TH2F *eff_multiRP_F1_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F1/h56_220_2017F1_all_2D");
        	TH2F *eff_multiRP_F2_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F2/h56_220_2017F2_all_2D");
        	TH2F *eff_multiRP_F3_56 = (TH2F*)eff_multiRP->Get("Pixel/2017/2017F3/h56_220_2017F3_all_2D");

        	TFile *eff_strips = new TFile("PreliminaryEfficiencies_July132020_1D2DMultiTrack.root");
        	TH2F *eff_strips_rddam_B_45 = (TH2F*)eff_strips->Get("Strips/2017/2017B/h45_2017B_all_2D");
        	TH2F *eff_strips_rddam_C_45 = (TH2F*)eff_strips->Get("Strips/2017/2017C/h45_2017C_all_2D");
        	TH2F *eff_strips_rddam_D_45 = (TH2F*)eff_strips->Get("Strips/2017/2017D/h45_2017D_all_2D");
        	TH2F *eff_strips_rddam_E_45 = (TH2F*)eff_strips->Get("Strips/2017/2017E/h45_2017E_all_2D");
        	TH2F *eff_strips_rddam_F_45 = (TH2F*)eff_strips->Get("Strips/2017/2017F/h45_2017F_all_2D");
        	TH2F *eff_strips_rddam_B_56 = (TH2F*)eff_strips->Get("Strips/2017/2017B/h56_2017B_all_2D");
        	TH2F *eff_strips_rddam_C_56 = (TH2F*)eff_strips->Get("Strips/2017/2017C/h56_2017C_all_2D");
        	TH2F *eff_strips_rddam_D_56 = (TH2F*)eff_strips->Get("Strips/2017/2017D/h56_2017D_all_2D");
        	TH2F *eff_strips_rddam_E_56 = (TH2F*)eff_strips->Get("Strips/2017/2017E/h56_2017E_all_2D");
        	TH2F *eff_strips_rddam_F_56 = (TH2F*)eff_strips->Get("Strips/2017/2017F/h56_2017F_all_2D");
        	TH1F *eff_strips_multitrack_B_45 = (TH1F*)eff_strips->Get("Strips/2017/2017B/h45multitrackeff_2017B_avg_RP3");
        	TH1F *eff_strips_multitrack_C_45 = (TH1F*)eff_strips->Get("Strips/2017/2017C/h45multitrackeff_2017C_avg_RP3");	
        	TH1F *eff_strips_multitrack_D_45 = (TH1F*)eff_strips->Get("Strips/2017/2017D/h45multitrackeff_2017D_avg_RP3");
        	TH1F *eff_strips_multitrack_E_45 = (TH1F*)eff_strips->Get("Strips/2017/2017E/h45multitrackeff_2017E_avg_RP3");
        	TH1F *eff_strips_multitrack_F_45 = (TH1F*)eff_strips->Get("Strips/2017/2017F/h45multitrackeff_2017F_avg_RP3");
        	TH1F *eff_strips_multitrack_B_56 = (TH1F*)eff_strips->Get("Strips/2017/2017B/h56multitrackeff_2017B_avg_RP103");
        	TH1F *eff_strips_multitrack_C_56 = (TH1F*)eff_strips->Get("Strips/2017/2017C/h56multitrackeff_2017C_avg_RP103");
        	TH1F *eff_strips_multitrack_D_56 = (TH1F*)eff_strips->Get("Strips/2017/2017D/h56multitrackeff_2017D_avg_RP103");
        	TH1F *eff_strips_multitrack_E_56 = (TH1F*)eff_strips->Get("Strips/2017/2017E/h56multitrackeff_2017E_avg_RP103");
        	TH1F *eff_strips_multitrack_F_56 = (TH1F*)eff_strips->Get("Strips/2017/2017F/h56multitrackeff_2017F_avg_RP103");

		TH1F *multitrack_45_F_B = (TH1F *)eff_strips->Get("Strips/2016/2016B/h45multitrackeff_2016B_avg_RP3");
		TH1F *multitrack_45_F_C = (TH1F *)eff_strips->Get("Strips/2016/2016C/h45multitrackeff_2016C_avg_RP3");
		TH1F *multitrack_45_F_G = (TH1F *)eff_strips->Get("Strips/2016/2016G/h45multitrackeff_2016G_avg_RP3");
		TH1F *multitrack_45_F_H = (TH1F *)eff_strips->Get("Strips/2016/2016H/h45multitrackeff_2016H_avg_RP3");
		TH1F *multitrack_45_N_B = (TH1F *)eff_strips->Get("Strips/2016/2016B/h45multitrackeff_2016B_avg_RP2");
		TH1F *multitrack_45_N_C = (TH1F *)eff_strips->Get("Strips/2016/2016C/h45multitrackeff_2016C_avg_RP2");
		TH1F *multitrack_45_N_G = (TH1F *)eff_strips->Get("Strips/2016/2016G/h45multitrackeff_2016G_avg_RP2");
		TH1F *multitrack_45_N_H = (TH1F *)eff_strips->Get("Strips/2016/2016H/h45multitrackeff_2016H_avg_RP2");
		TH1F *multitrack_56_F_B = (TH1F *)eff_strips->Get("Strips/2016/2016B/h56multitrackeff_2016B_avg_RP103");
		TH1F *multitrack_56_F_C = (TH1F *)eff_strips->Get("Strips/2016/2016C/h56multitrackeff_2016C_avg_RP103");
		TH1F *multitrack_56_F_G = (TH1F *)eff_strips->Get("Strips/2016/2016G/h56multitrackeff_2016G_avg_RP103");
		TH1F *multitrack_56_F_H = (TH1F *)eff_strips->Get("Strips/2016/2016H/h56multitrackeff_2016H_avg_RP103");
		TH1F *multitrack_56_N_B = (TH1F *)eff_strips->Get("Strips/2016/2016B/h56multitrackeff_2016B_avg_RP102");
		TH1F *multitrack_56_N_C = (TH1F *)eff_strips->Get("Strips/2016/2016C/h56multitrackeff_2016C_avg_RP102");
		TH1F *multitrack_56_N_G = (TH1F *)eff_strips->Get("Strips/2016/2016G/h56multitrackeff_2016G_avg_RP102");
		TH1F *multitrack_56_N_H = (TH1F *)eff_strips->Get("Strips/2016/2016H/h56multitrackeff_2016H_avg_RP102");
		TH2F *RadDam_45_F_B = (TH2F *)eff_strips->Get("Strips/2016/2016B/h45_2016B_RP3_all_2D");
		TH2F *RadDam_45_N_B = (TH2F *)eff_strips->Get("Strips/2016/2016B/h45_2016B_RP2_all_2D");
		TH2F *RadDam_56_F_B = (TH2F *)eff_strips->Get("Strips/2016/2016B/h56_2016B_RP103_all_2D");
		TH2F *RadDam_56_N_B = (TH2F *)eff_strips->Get("Strips/2016/2016B/h56_2016B_RP102_all_2D");
		TH2F *RadDam_45_F_C = (TH2F *)eff_strips->Get("Strips/2016/2016C/h45_2016C_RP3_all_2D");
		TH2F *RadDam_45_N_C = (TH2F *)eff_strips->Get("Strips/2016/2016C/h45_2016C_RP2_all_2D");
		TH2F *RadDam_56_F_C = (TH2F *)eff_strips->Get("Strips/2016/2016C/h56_2016C_RP103_all_2D");
		TH2F *RadDam_56_N_C = (TH2F *)eff_strips->Get("Strips/2016/2016C/h56_2016C_RP102_all_2D");
		TH2F *RadDam_45_F_G = (TH2F *)eff_strips->Get("Strips/2016/2016G/h45_2016G_RP3_all_2D");
		TH2F *RadDam_45_N_G = (TH2F *)eff_strips->Get("Strips/2016/2016G/h45_2016G_RP2_all_2D");
		TH2F *RadDam_56_F_G = (TH2F *)eff_strips->Get("Strips/2016/2016G/h56_2016G_RP103_all_2D");
		TH2F *RadDam_56_N_G = (TH2F *)eff_strips->Get("Strips/2016/2016G/h56_2016G_RP102_all_2D");


		edm::EDGetTokenT< double > prefweight_token;
                edm::EDGetTokenT< double > prefweightup_token;
                edm::EDGetTokenT< double > prefweightdown_token;

};

void MakeNTuple::Get_t_and_xi(const TLorentzVector* proton,double& t,double& xi) {
	set_BeamEnergy(13000/2.);
	t = 0.;
	xi = -1.;
	if (!proton) return;
	double mom    = proton->P();
	if (mom>fBeamMomentum) mom=fBeamMomentum;
	double energy = proton->E();
	double theta  = (proton->Pz()>0)?proton->Theta():TMath::Pi()-proton->Theta();
	t      = -2.*(ProtonMassSQ-fBeamEnergy*energy+fBeamMomentum*mom*cos(theta));
	xi     = (1.0-energy/fBeamEnergy);
}
double MakeNTuple::GetWeightLeptons(double pt, double eta, TH2F* hist2D){
	double weight=1.;
	if (pt<120)  weight = hist2D->GetBinContent(hist2D->GetXaxis()->FindBin(abs(eta)),hist2D->GetYaxis()->FindBin(pt));
	else         weight = hist2D->GetBinContent(hist2D->GetXaxis()->FindBin(abs(eta)),hist2D->GetYaxis()->FindBin(119.));
	return weight;
}
MakeNTuple::MakeNTuple(const edm::ParameterSet& iConfig):
	  muonsToken (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
	, electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
	, MetToken (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MET")))
	, PFCandToken (consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCand")))
	, PileupSumInfoInputTag (consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupSumInfoInputTag")))
	, verticesToken (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
	, recoProtonsSingleRPToken_(consumes<std::vector<reco::ForwardProton>>(iConfig.getParameter<edm::InputTag>("ppsRecoProtonSingleRPTag")))
	, recoProtonsMultiRPToken_(consumes<std::vector<reco::ForwardProton>>(iConfig.getParameter<edm::InputTag>("ppsRecoProtonMultiRPTag")))
	, LumiWeights_(0)
	, MC(iConfig.getParameter<bool>("MC"))
	, MC_Signal(iConfig.getParameter<bool>("MC_Signal"))
	, DATA(iConfig.getParameter<bool>("DATA"))
	, tag(iConfig.getParameter<std::string>("tag_"))
        , preVFP(iConfig.getParameter<bool>("preVFP_"))
        , postVFP(iConfig.getParameter<bool>("postVFP_"))
	, triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults")))
	, genToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genTag")))
	, mcEventToken (consumes<edm::HepMCProduct>(iConfig.getUntrackedParameter<edm::InputTag>("MCEvent",std::string(""))))
	, hltPrescaleProvider_(iConfig, consumesCollector(), *this)
	, prefweight_token (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb")))
	, prefweightup_token  (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp")))
	, prefweightdown_token (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown")))
{
	//now do what ever initialization is needed

	usesResource("TFileService");
	edm::Service<TFileService> fs;
	EventBranchs = fs->make<TTree>( "Events","Events" );
	EventBranchs->Branch("Run", &Run, "Run/I");
	EventBranchs->Branch("LumiSection", &LumiSection, "LumiSection/I");
	EventBranchs->Branch("BX", &BX, "BX/I");
	EventBranchs->Branch("xangle", &xangle, "xangle/I");

	EventBranchs->Branch("EventNum", &EventNum, "EventNum/I");
	EventBranchs->Branch("nVtx", &nVtx, "nVtx/I");

	EventBranchs->Branch("PUWeight",&PUWeight,"PUWeight/D");
	EventBranchs->Branch("nPU",&nPU,"nPU/D");

	const edm::InputTag& PileupSumInfoInputTags = edm::InputTag( "slimmedAddPileupInfo" ) ;

	LumiWeights_ = new edm::LumiReWeighting("PileupMC.root", "MyDataPileupHistogram.root", "input_Event/N_TrueInteractions", "pileup", PileupSumInfoInputTags);

	EventBranchs->Branch("ArmF_xi_gen","std::vector<double>",&ArmF_xi_gen);
        EventBranchs->Branch("ArmB_xi_gen","std::vector<double>",&ArmB_xi_gen);
        EventBranchs->Branch("ArmF_t_gen","std::vector<double>",&ArmF_t_gen);
        EventBranchs->Branch("ArmB_t_gen","std::vector<double>",&ArmB_t_gen);
        EventBranchs->Branch("ArmF_thx_gen","std::vector<double>",&ArmF_thx_gen);
        EventBranchs->Branch("ArmB_thx_gen","std::vector<double>",&ArmB_thx_gen);
        EventBranchs->Branch("ArmF_thy_gen","std::vector<double>",&ArmF_thy_gen);
        EventBranchs->Branch("ArmB_thy_gen","std::vector<double>",&ArmB_thy_gen);

	EventBranchs->Branch("ProtCand_xi","std::vector<double>",&ProtCand_xi);
	EventBranchs->Branch("ProtCand_t","std::vector<double>",&ProtCand_t);
        EventBranchs->Branch("ProtCand_ismultirp","std::vector<double>",&ProtCand_ismultirp);
	EventBranchs->Branch("ProtCand_weight",&ProtCand_weight,"ProtCand_weight/D");
        EventBranchs->Branch("ProtCand_arm","std::vector<double>",&ProtCand_arm);

        EventBranchs->Branch("lepton_P_pT", &lepton_P_pT, "lepton_P_pT/D");
        EventBranchs->Branch("lepton_N_pT", &lepton_N_pT, "lepton_N_pT/D");
        EventBranchs->Branch("lepton_P_index", &lepton_P_index, "lepton_P_index/D");
        EventBranchs->Branch("lepton_N_index", &lepton_N_index, "lepton_N_index/D");
        EventBranchs->Branch("lepton_P_ismuon", &lepton_P_ismuon, "lepton_P_ismuon/D");
        EventBranchs->Branch("lepton_N_ismuon", &lepton_N_ismuon, "lepton_N_ismuon/D"); 
        EventBranchs->Branch("lepton_P_iselectron", &lepton_P_iselectron, "lepton_P_iselectron/D");
        EventBranchs->Branch("lepton_N_iselectron", &lepton_N_iselectron, "lepton_N_iselectron/D");
        EventBranchs->Branch("lepton_P_eta", &lepton_P_eta, "lepton_P_eta/D");
        EventBranchs->Branch("lepton_N_eta", &lepton_N_eta, "lepton_N_eta/D");
        EventBranchs->Branch("lepton_P_phi", &lepton_P_phi, "lepton_P_phi/D");
        EventBranchs->Branch("lepton_N_phi", &lepton_N_phi, "lepton_N_phi/D");
        EventBranchs->Branch("lepton_P_E", &lepton_P_E, "lepton_P_E/D");
        EventBranchs->Branch("lepton_N_E", &lepton_N_E, "lepton_N_E/D");
        EventBranchs->Branch("lepton_P_px", &lepton_P_px, "lepton_P_px/D");
        EventBranchs->Branch("lepton_N_px", &lepton_N_px, "lepton_N_px/D");
        EventBranchs->Branch("lepton_P_py", &lepton_P_py, "lepton_P_py/D");
        EventBranchs->Branch("lepton_N_py", &lepton_N_py, "lepton_N_py/D");
        EventBranchs->Branch("lepton_P_pz", &lepton_P_pz, "lepton_P_pz/D");
        EventBranchs->Branch("lepton_N_pz", &lepton_N_pz, "lepton_N_pz/D");
        EventBranchs->Branch("lepton_P_vtxZ", &lepton_P_vtxZ, "lepton_P_vtxZ/D");
        EventBranchs->Branch("lepton_N_vtxZ", &lepton_N_vtxZ, "lepton_N_vtxZ/D");
        EventBranchs->Branch("lepton_P_charge", &lepton_P_charge, "lepton_P_charge/D");
        EventBranchs->Branch("lepton_N_charge", &lepton_N_charge, "lepton_N_charge/D");
        EventBranchs->Branch("lepton_P_PFBasedIso", &lepton_P_PFBasedIso, "lepton_P_PFBasedIso/D");
        EventBranchs->Branch("lepton_N_PFBasedIso", &lepton_N_PFBasedIso, "lepton_N_PFBasedIso/D");
        EventBranchs->Branch("lepton_P_isTight", &lepton_P_isTight, "lepton_P_isTight/D");
        EventBranchs->Branch("lepton_N_isTight", &lepton_N_isTight, "lepton_N_isTight/D");
        EventBranchs->Branch("lepton_P_isMedium", &lepton_P_isMedium, "lepton_P_isMedium/D");
        EventBranchs->Branch("lepton_N_isMedium", &lepton_N_isMedium, "lepton_N_isMedium/D");

        EventBranchs->Branch("lepton_P_gen_pT", &lepton_P_gen_pT, "lepton_P_gen_pT/D");
        EventBranchs->Branch("lepton_N_gen_pT", &lepton_N_gen_pT, "lepton_N_gen_pT/D");
        EventBranchs->Branch("lepton_P_gen_eta", &lepton_P_gen_pT, "lepton_P_gen_eta/D");
        EventBranchs->Branch("lepton_N_gen_eta", &lepton_N_gen_pT, "lepton_N_gen_eta/D");
        EventBranchs->Branch("lepton_P_gen_phi", &lepton_P_gen_phi, "lepton_P_gen_phi/D");
        EventBranchs->Branch("lepton_N_gen_phi", &lepton_N_gen_phi, "lepton_N_gen_phi/D");
        EventBranchs->Branch("lepton_P_gen_E", &lepton_P_gen_E, "lepton_P_gen_E/D");
        EventBranchs->Branch("lepton_N_gen_E", &lepton_N_gen_E, "lepton_N_gen_E/D");
        EventBranchs->Branch("lepton_P_gen_px", &lepton_P_gen_px, "lepton_P_gen_px/D");
        EventBranchs->Branch("lepton_N_gen_px", &lepton_N_gen_px, "lepton_N_gen_px/D");
        EventBranchs->Branch("lepton_P_gen_py", &lepton_P_gen_py, "lepton_P_gen_py/D");
        EventBranchs->Branch("lepton_N_gen_py", &lepton_N_gen_py, "lepton_N_gen_py/D");
        EventBranchs->Branch("lepton_P_gen_pz", &lepton_P_gen_pz, "lepton_P_gen_pz/D");
        EventBranchs->Branch("lepton_N_gen_pz", &lepton_N_gen_pz, "lepton_N_gen_pz/D");
        EventBranchs->Branch("lepton_P_gen_ismuon", &lepton_P_gen_ismuon, "lepton_P_gen_ismuon/D");
        EventBranchs->Branch("lepton_N_gen_ismuon", &lepton_N_gen_ismuon, "lepton_N_gen_ismuon/D");
        EventBranchs->Branch("lepton_P_gen_iselectron", &lepton_P_gen_ismuon, "lepton_P_gen_iselectron/D");
        EventBranchs->Branch("lepton_N_gen_iselectron", &lepton_N_gen_ismuon, "lepton_N_gen_iselectron/D");

        EventBranchs->Branch("neut_P_gen_pT", &neut_P_gen_pT, "neut_P_gen_pT/D");
        EventBranchs->Branch("neut_N_gen_pT", &neut_N_gen_pT, "neut_N_gen_pT/D");
        EventBranchs->Branch("neut_P_gen_eta", &neut_P_gen_eta, "neut_P_gen_eta/D");
        EventBranchs->Branch("neut_N_gen_eta", &neut_N_gen_eta, "neut_N_gen_eta/D");
        EventBranchs->Branch("neut_P_gen_phi", &neut_P_gen_phi, "neut_P_gen_phi/D");
        EventBranchs->Branch("neut_N_gen_phi", &neut_N_gen_phi, "neut_N_gen_phi/D");
        EventBranchs->Branch("neut_P_gen_E", &neut_P_gen_E, "neut_P_gen_E/D");
        EventBranchs->Branch("neut_N_gen_E", &neut_N_gen_E, "neut_N_gen_E/D");
        EventBranchs->Branch("neut_P_gen_px", &neut_P_gen_px, "neut_P_gen_px/D");
        EventBranchs->Branch("neut_N_gen_px", &neut_N_gen_px, "neut_N_gen_px/D");
        EventBranchs->Branch("neut_P_gen_py", &neut_P_gen_py, "neut_P_gen_py/D");
        EventBranchs->Branch("neut_N_gen_py", &neut_N_gen_py, "neut_N_gen_py/D");
        EventBranchs->Branch("neut_P_gen_pz", &neut_P_gen_pz, "neut_P_gen_pz/D");
        EventBranchs->Branch("neut_N_gen_pz", &neut_N_gen_pz, "neut_N_gen_pz/D");

        EventBranchs->Branch("Wp_gen_pT", &Wp_gen_pT, "Wp_gen_pT/D");
        EventBranchs->Branch("Wm_gen_pT", &Wm_gen_pT, "Wm_gen_pT/D");
        EventBranchs->Branch("Wp_gen_eta", &Wp_gen_eta, "Wp_gen_eta/D");
        EventBranchs->Branch("Wm_gen_eta", &Wm_gen_eta, "Wm_gen_eta/D");
        EventBranchs->Branch("Wp_gen_phi", &Wp_gen_phi, "Wp_gen_phi/D");
        EventBranchs->Branch("Wm_gen_phi", &Wm_gen_phi, "Wm_gen_phi/D");
        EventBranchs->Branch("Wp_gen_E", &Wp_gen_E, "Wp_gen_E/D");
        EventBranchs->Branch("Wm_gen_E", &Wm_gen_E, "Wm_gen_E/D");
        EventBranchs->Branch("Wp_gen_px", &Wp_gen_px, "Wp_gen_px/D");
        EventBranchs->Branch("Wm_gen_px", &Wm_gen_px, "Wm_gen_px/D");
        EventBranchs->Branch("Wp_gen_py", &Wp_gen_py, "Wp_gen_py/D");
        EventBranchs->Branch("Wm_gen_py", &Wm_gen_py, "Wm_gen_py/D");
        EventBranchs->Branch("Wp_gen_pz", &Wp_gen_pz, "Wp_gen_pz/D");
        EventBranchs->Branch("Wm_gen_pz", &Wm_gen_pz, "Wm_gen_pz/D");

	EventBranchs->Branch("prefweight",&prefweight, "prefweight/D");
	EventBranchs->Branch("weight",&weight, "weight/D");

	EventBranchs->Branch("pfphi","std::vector<double>",&pfphi);
	EventBranchs->Branch("pfeta","std::vector<double>",&pfeta);

	EventBranchs->Branch("pffromPV","std::vector<int>",&pffromPV);
	EventBranchs->Branch("pfdz","std::vector<double>",&pfdz);
	EventBranchs->Branch("pfpt","std::vector<double>",&pfpt);

	EventBranchs->Branch("METPx", &METPx, "METPx/D");
	EventBranchs->Branch("METPy", &METPy, "METPy/D");
	EventBranchs->Branch("METPt", &METPt, "METPt/D");
	EventBranchs->Branch("METphi", &METphi, "METphi/D");

	EventBranchs->Branch("channel",&channel,"channel/I");

	EventBranchs->Branch("HLT_name","std::vector<std::string>",&HLT_name);
	EventBranchs->Branch("HLT_pass","std::vector<bool>",&HLT_pass);
	EventBranchs->Branch("HLT_prescale","std::vector<int>",&HLT_prescale);

}


MakeNTuple::~MakeNTuple()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

	f_SFmuon_ID_preVFP->Close();
        f_SFmuon_ID_postVFP->Close();
        f_SFmuon_ID_2017->Close();
        f_SFmuon_ID_2018->Close();
        f_SFmuon_ISO_preVFP->Close();
        f_SFmuon_ISO_postVFP->Close();
        f_SFmuon_ISO_2017->Close();
        f_SFmuon_ISO_2018->Close();
        f_SFelectron_ID_preVFP->Close();
        f_SFelectron_ID_postVFP->Close();
        f_SFelectron_ID_2017->Close();
        f_SFelectron_ID_2018->Close();
	eff_strips->Close();
	eff_multiRP->Close();
	eff_radiation->Close();
        f_SFmuon_ID_preVFP->Delete();
        f_SFmuon_ID_postVFP->Delete();
        f_SFmuon_ID_2017->Delete();
        f_SFmuon_ID_2018->Delete();
        f_SFmuon_ISO_preVFP->Delete();
        f_SFmuon_ISO_postVFP->Delete();
        f_SFmuon_ISO_2017->Delete();
        f_SFmuon_ISO_2018->Delete();
        f_SFelectron_ID_preVFP->Delete();
        f_SFelectron_ID_postVFP->Delete();
        f_SFelectron_ID_2017->Delete();
        f_SFelectron_ID_2018->Delete();
        eff_strips->Delete();
        eff_multiRP->Delete();
        eff_radiation->Delete();
//                                                                                                                      	
}


//
// member functions
//

// ------------ method called for each event  ------------
	void
MakeNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
        //std::cout << "begin" << std::endl;
	f_SFmuon_ID_preVFP->GetObject("NUM_TightID_DEN_TrackerMuons_abseta_pt;1",h_SFmuon_ID_preVFP);
        f_SFmuon_ISO_preVFP->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt;1",h_SFmuon_ISO_preVFP);
        f_SFmuon_ID_postVFP->GetObject("NUM_TightID_DEN_TrackerMuons_abseta_pt;1",h_SFmuon_ID_postVFP);
        f_SFmuon_ISO_postVFP->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt;1",h_SFmuon_ISO_postVFP);
        f_SFmuon_ID_2017->GetObject("NUM_TightID_DEN_TrackerMuons_abseta_pt;1",h_SFmuon_ID_2017);
        f_SFmuon_ISO_2017->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt;1",h_SFmuon_ISO_2017);
        f_SFmuon_ID_2018->GetObject("NUM_TightID_DEN_TrackerMuons_abseta_pt;1",h_SFmuon_ID_2018);
        f_SFmuon_ISO_2018->GetObject("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt;1",h_SFmuon_ISO_2018);
        f_SFelectron_ID_preVFP->GetObject("EGamma_SF2D;1",h_SFelectron_ID_preVFP);
        f_SFelectron_ID_postVFP->GetObject("EGamma_SF2D;1",h_SFelectron_ID_postVFP);
        f_SFelectron_ID_2017->GetObject("EGamma_SF2D;1",h_SFelectron_ID_2017);
        f_SFelectron_ID_2018->GetObject("EGamma_SF2D;1",h_SFelectron_ID_2018);

	ArmF_xi_gen->clear();
	ArmB_xi_gen->clear();
	ArmF_t_gen->clear();
	ArmB_t_gen->clear();
	ArmF_thx_gen->clear();
	ArmB_thx_gen->clear();
	ArmF_thy_gen->clear();
	ArmB_thy_gen->clear();

	ProtCand_xi->clear();
	ProtCand_t->clear();
	ProtCand_ThX->clear();
	ProtCand_ThY->clear();
	ProtCand_rpid->clear();
	ProtCand_arm->clear();
	ProtCand_ismultirp->clear();
        ProtCand_x->clear();
        ProtCand_y->clear();
        ProtCand_xn->clear();
        ProtCand_yn->clear();
        ProtCand_xf->clear();
        ProtCand_yf->clear();

	pfphi->clear();
	pfeta->clear();
	pffromPV->clear();
	pfdz->clear();
	pfpt->clear();
	HLT_name->clear();
	HLT_pass->clear();
	HLT_prescale->clear();

        bool year2017 = (tag=="MC_2017") || (tag=="MC_signal_2017") || (tag == "data_2017");
        bool year2018 = (tag=="MC_2018") || (tag=="MC_signal_2018") || (tag == "data_2018");

	muSF = 1.;
	elSF = 1.;	
	channel = -1; // -1 is a non available value at the end of analysis, the allowed calues are 0, 1 or 2.
        lepton_P_pT = 0.; lepton_N_pT = 0.;
        lepton_P_index = -1.; lepton_N_index = -1.;
        lepton_P_ismuon = 0.; lepton_N_ismuon = 0.; lepton_P_iselectron = 0.; lepton_N_iselectron = 0.;
        lepton_P_eta = 0.; lepton_N_eta = 0.; lepton_P_phi = 0.; lepton_N_phi = 0.;
        lepton_P_E = 0.; lepton_N_E = 0.; lepton_P_px = 0.; lepton_N_px = 0.; 
	lepton_P_py = 0.; lepton_N_py = 0.; lepton_P_pz = 0.; lepton_N_pz = 0.;
        lepton_N_vtxZ = 0.; lepton_P_vtxZ = 0.; lepton_N_charge = 0.; lepton_P_charge = 0.;
        lepton_N_PFBasedIso = 0.; lepton_P_PFBasedIso = 0.; lepton_N_isTight = 0.; lepton_P_isTight = 0.;
        lepton_N_isMedium = 0.; lepton_P_isMedium = 0.; 
        lepton_P_gen_pT = 0.; lepton_N_gen_pT = 0.;
        lepton_P_gen_ismuon = 0.; lepton_N_gen_ismuon = 0.; lepton_P_gen_iselectron = 0.; lepton_N_gen_iselectron = 0.;
        lepton_P_gen_eta = 0.; lepton_N_gen_eta = 0.; lepton_P_gen_phi = 0.; lepton_N_gen_phi = 0.;
        lepton_P_gen_E = 0.; lepton_N_gen_E = 0.; lepton_P_gen_px = 0.; lepton_N_gen_px = 0.;
        lepton_P_gen_py = 0.; lepton_N_gen_py = 0.; lepton_P_gen_pz = 0.; lepton_N_gen_pz = 0.;
	prefweight=1.; prefweightup=1.; prefweightdown=1.;
        neut_P_gen_pT = 0.; neut_N_gen_pT = 0.;
        neut_P_gen_eta = 0.; neut_N_gen_eta = 0.; neut_P_gen_phi = 0.; neut_N_gen_phi = 0.;
        neut_P_gen_E = 0.; neut_N_gen_E = 0.; neut_P_gen_px = 0.; neut_N_gen_px = 0.;
        neut_P_gen_py = 0.; neut_N_gen_py = 0.; neut_P_gen_pz = 0.; neut_N_gen_pz = 0.;
        Wp_gen_eta = 0.; Wm_gen_eta = 0.; Wp_gen_phi = 0.; Wm_gen_phi = 0.;
        Wp_gen_E = 0.; Wm_gen_E = 0.; Wp_gen_px = 0.; Wm_gen_px = 0.;
        Wp_gen_py = 0.; Wm_gen_py = 0.; Wp_gen_pz = 0.; Wm_gen_pz = 0.;
	using namespace edm;
	using namespace std;
	//std::cout << "read tokens" << std::endl;
	iEvent.getByToken(muonsToken, muons);
	iEvent.getByToken(electronsToken, electrons);
	iEvent.getByToken(MetToken, MET);
	iEvent.getByToken(PFCandToken, PFCand);
	iEvent.getByToken(verticesToken, vertices);

	if (!year2018 || MC){
		edm::Handle< double > theprefweight;
		iEvent.getByToken(prefweight_token, theprefweight ) ;
		prefweight =(*theprefweight);

		edm::Handle< double > theprefweightup;
		iEvent.getByToken(prefweightup_token, theprefweightup ) ;
		prefweightup =(*theprefweightup);

		edm::Handle< double > theprefweightdown;
		iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
		prefweightdown =(*theprefweightdown);
	}
	if(MC){
		iEvent.getByToken(PileupSumInfoInputTag, PupInfo);
	}

	//	if(MC_Signal || DATA){
	edm::Handle<vector<reco::ForwardProton>> recoMultiRPProtons;
	iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
	edm::Handle<vector<reco::ForwardProton>> recoSingleRPProtons;
	iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);
	//	}

	// Event Info
	BX = iEvent.bunchCrossing();
	Run = iEvent.id().run();
	LumiSection = iEvent.luminosityBlock();
	EventNum = iEvent.id().event();
	//std::cout << "read vertex" << std::endl;
	nVtx = vertices->size();

	vtx_isValid = false;
	vtx_isFake = true;
	vtx_z = -999;

	if (nVtx>0){
		vtx_isValid = vertices->at(0).isValid();
		vtx_isFake = vertices->at(0).isFake();
		vtx_z = vertices->at(0).z();
	}

	// X A N G L E   I N F O
	/*	if(DATA){
		edm::ESHandle<LHCInfo> pSetup;
		const string label = "";
		iSetup.get<LHCInfoRcd>().get(label, pSetup);
		const LHCInfo* pInfo = pSetup.product();
		std::cout << pInfo->crossingAngle() << std::endl;
		xangle = pInfo->crossingAngle();
		}
		*/


	// F O R  P I L E U P  R E W E I G T I N G
	//std::cout << "pile-up reweight" << std::endl;
	if(MC){
		const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
		if (LumiWeights_) PUWeight = LumiWeights_->weight( (*iEventB) );
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		int npv = -1;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			int BXs = PVI->getBunchCrossing();
			if(BXs == 0) {
				npv = PVI->getTrueNumInteractions();
				continue;
			}
		}
		nPU = npv;
	}
	//std::cout << "read trigger info" << std::endl;
	// F O R  T R I G G E R  I N F O
	// The trigger information was included by Joao Pedro at signal samples
	edm::Handle<edm::TriggerResults> hltResults;
	iEvent.getByToken(triggerResultsToken_, hltResults);
	const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);
	for (size_t i=0; i<trigNames.size(); i++) {
	   for (size_t j =0; j<HLT_list.size(); j++) {
	      if (trigNames.triggerNames().at(i).find(HLT_list.at(j)) != std::string::npos){		
 	           HLT_name->push_back(trigNames.triggerNames().at(i));
		   //std::cout << "trigNames.triggerNames().at(i): " << trigNames.triggerNames().at(i) << std::endl;
	  	   HLT_pass->push_back(hltResults->accept(i));
		   HLT_prescale->push_back(hltPrescaleProvider_.prescaleValue(iEvent, iSetup,trigNames.triggerName(i)));
	      }
  	   }
	}



	// FOR GEN INFO	
	if (MC_Signal) {

           iEvent.getByToken( mcEventToken, EvtHandle );
           const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;

	   //iEvent.getByToken(genToken, genHandle);
	   std::vector<math::XYZTLorentzVector> protonCTPPS;
	   protonCTPPS.clear();

	   //for (auto&& genPart : *(genHandle.product())){
	   for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin(); i != Evt->particles_end();++i) {
	      int myId = (*i)->pdg_id();
	      if(myId==2212){
	 	//const pat::GenParticle& mcProton = (*genPart);
		HepMC::FourVector momentum=(*i)->momentum();
	        const HepMC::FourVector p((*i)->momentum());
		protonCTPPS.push_back(math::XYZTLorentzVector(p.x(),p.y(),p.z(),p.t()));
		double t,xi;
		double px = momentum.x();
		double py = momentum.y();
		double pz = momentum.z();
		double thx = atan(px/pz);
		double thy = atan(py/pz);
		double e = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
		TLorentzVector* proton = new TLorentzVector(px,py,pz,e);
				
		if((*i)->status() == 1) {
		   Get_t_and_xi(proton,t,xi);
		   if (pz > 5000.0) { // armId 0
		      std::cout << "xigen: " << xi << std::endl;
		      ArmF_xi_gen->push_back(xi);
		      ArmF_t_gen->push_back(t);
		      ArmF_thx_gen->push_back(thx);
		      ArmF_thy_gen->push_back(thy);
		   }
		   if (pz < 5000.0) { // armId 1
		      std::cout << "xigen: " << xi << std::endl;
		      ArmB_xi_gen->push_back(xi);
		      ArmB_t_gen->push_back(t);
		      ArmB_thx_gen->push_back(thx);
		      ArmB_thy_gen->push_back(thy);
		   }
		 }
	       }
	    
	       int pId = abs((*i)->pdg_id());
	     	
               if ( pId == 24){
		  for ( HepMC::GenVertex::particle_iterator dau  =(*i)->end_vertex()->particles_begin(HepMC::children); dau != (*i)->end_vertex()->particles_end(HepMC::children); ++dau ) {
		     if ( abs((*dau)->pdg_id()) == 24 ){
			for ( HepMC::GenVertex::particle_iterator ddau  =(*dau)->end_vertex()->particles_begin(HepMC::children); ddau != (*dau)->end_vertex()->particles_end(HepMC::children); ++ddau ) {
			    if ((*ddau)->pdg_id()==-13 || (*ddau)->pdg_id()==-11){
				if ( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5) > Wp_gen_pT){
				   Wp_gen_px = (*dau)->momentum().x();
		                   Wp_gen_py = (*dau)->momentum().y();
			           Wp_gen_pz = (*dau)->momentum().z();
			           Wp_gen_pT = pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5);
			           Wp_gen_E = pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y()+(*dau)->momentum().z()*(*dau)->momentum().z(),0.5);
		                   Wp_gen_phi = (*dau)->momentum().phi();
				   Wp_gen_eta = (*dau)->momentum().eta();
			        }
				if ( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5) > lepton_P_gen_pT){
			           if ((*ddau)->pdg_id()==-13) {lepton_P_gen_ismuon=1;}
				   if ((*ddau)->pdg_id()==-11) {lepton_P_gen_iselectron=1;} 
				   lepton_P_gen_px = (*ddau)->momentum().x();
                                   lepton_P_gen_py = (*ddau)->momentum().y();
                                   lepton_P_gen_pz = (*ddau)->momentum().z();
                                   lepton_P_gen_pT = pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5);
                                   lepton_P_gen_E=pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y()+(*ddau)->momentum().z()*(*ddau)->momentum().z(),0.5);
                                   lepton_P_gen_phi = (*ddau)->momentum().phi();
                                   lepton_P_gen_eta = (*ddau)->momentum().eta();
				}
				for(HepMC::GenVertex::particle_iterator ddau2=(*dau)->end_vertex()->particles_begin(HepMC::children);ddau2!=(*dau)->end_vertex()->particles_end(HepMC::children);++ddau2){
                                   if (abs((*ddau2)->pdg_id()) == 14 || abs((*ddau2)->pdg_id()) == 12){
				     if (pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y(),0.5) > neut_P_gen_pT){
                                        neut_P_gen_px = (*ddau2)->momentum().x();
                                        neut_P_gen_py = (*ddau2)->momentum().y();
                                        neut_P_gen_pz = (*ddau2)->momentum().z();
                                        neut_P_gen_pT = pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y(),0.5);
                                        neut_P_gen_E = pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y()+(*ddau2)->momentum().z()*(*ddau2)->momentum().z(),0.5);
                                        neut_P_gen_phi = (*ddau2)->momentum().phi();
                                        neut_P_gen_eta = (*ddau2)->momentum().eta();
				      }
                                   }
                                }
			     }
                            if ((*ddau)->pdg_id()==13 || (*ddau)->pdg_id()==11){
                                if ( pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5) > Wm_gen_pT){
                                   Wm_gen_px = (*dau)->momentum().x();
                                   Wm_gen_py = (*dau)->momentum().y();
                                   Wm_gen_pz = (*dau)->momentum().z();
                                   Wm_gen_pT = pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y(),0.5);
                                   Wm_gen_E = pow((*dau)->momentum().x()*(*dau)->momentum().x()+(*dau)->momentum().y()*(*dau)->momentum().y()+(*dau)->momentum().z()*(*dau)->momentum().z(),0.5);
                                   Wm_gen_phi = (*dau)->momentum().phi();
                                   Wm_gen_eta = (*dau)->momentum().eta();
                                }
                                if ( pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5) > lepton_N_gen_pT){
                                   if ((*ddau)->pdg_id()==13) {lepton_N_gen_ismuon=1;}
                                   if ((*ddau)->pdg_id()==11) {lepton_N_gen_iselectron=1;}
				   lepton_N_gen_px = (*ddau)->momentum().x();
                                   lepton_N_gen_py = (*ddau)->momentum().y();
                                   lepton_N_gen_pz = (*ddau)->momentum().z();
                                   lepton_N_gen_pT = pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y(),0.5);
                                   lepton_N_gen_E=pow((*ddau)->momentum().x()*(*ddau)->momentum().x()+(*ddau)->momentum().y()*(*ddau)->momentum().y()+(*ddau)->momentum().z()*(*ddau)->momentum().z(),0.5);
                                   lepton_N_gen_phi = (*ddau)->momentum().phi();
                                   lepton_N_gen_eta = (*ddau)->momentum().eta();
                                }
                                for(HepMC::GenVertex::particle_iterator ddau2=(*dau)->end_vertex()->particles_begin(HepMC::children);ddau2!=(*dau)->end_vertex()->particles_end(HepMC::children);++ddau2){
                                   if (abs((*ddau2)->pdg_id()) == 14 || abs((*ddau2)->pdg_id()) == 12){
                                     if (pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y(),0.5) > neut_N_gen_pT){
                                        neut_N_gen_px = (*ddau2)->momentum().x();
                                        neut_N_gen_py = (*ddau2)->momentum().y();
                                        neut_N_gen_pz = (*ddau2)->momentum().z();
                                        neut_N_gen_pT = pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y(),0.5);
                                        neut_N_gen_E = pow((*ddau2)->momentum().x()*(*ddau2)->momentum().x()+(*ddau2)->momentum().y()*(*ddau2)->momentum().y()+(*ddau2)->momentum().z()*(*ddau2)->momentum().z(),0.5);
                                        neut_N_gen_phi = (*ddau2)->momentum().phi();
                                        neut_N_gen_eta = (*ddau2)->momentum().eta();
                                      }
                                   }
                                }
                             }
		        } //for ( HepMC::GenVertex::particle_iterator ddau....  
		      } // if ( abs((*dau)->pdg_id()) == 24 ){ 
		    } //for ( HepMC::GenVertex::particle_iterator dau ...
	       } // if ( pId == 24){
	     } // for(HepMC::GenEvent::particle_const_iterator i=Evt->particles_begin()...
	 } // if (MC_Signal) {

	//std::cout << "read gen info" << std::endl;

        if (MC && !MC_Signal) {
	   edm::Handle<reco::GenParticleCollection> genColl; 
	   iEvent.getByToken(genToken, genColl);
	   
	   for (auto&& genPart : *(genColl.product())) { 
	   //std::cout << genPart.pdgId() << std::endl;
 	   
	      if (genPart.status()==1) {
	        if (genPart.pdgId() == 13 || genPart.pdgId() == 11) { 
		  if (pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5) > lepton_N_gen_pT){
		     if (genPart.pdgId()==13) {lepton_N_gen_ismuon=1;lepton_P_gen_iselectron=0;}
                     if (genPart.pdgId()==11) {lepton_P_gen_ismuon=0;lepton_N_gen_iselectron=1;}
		     lepton_N_gen_px = genPart.px();
                     lepton_N_gen_py = genPart.py();
		     lepton_N_gen_pz = genPart.pz();	
                     lepton_N_gen_pT = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5);
                     lepton_N_gen_E = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5);
                     lepton_N_gen_phi = genPart.phi();
                     lepton_N_gen_eta = genPart.eta();
		   }
	         } 
                if (genPart.pdgId() == -13 || genPart.pdgId() == -11) {
                  if (pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5) > lepton_P_gen_pT){
                     if (genPart.pdgId()==-13) {lepton_P_gen_ismuon=1;lepton_P_gen_iselectron=0;}
                     if (genPart.pdgId()==-11) {lepton_P_gen_ismuon=0;lepton_P_gen_iselectron=1;}
                     lepton_P_gen_px = genPart.px();
                     lepton_P_gen_py = genPart.py();
                     lepton_P_gen_pz = genPart.pz();
                     lepton_P_gen_pT = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5);
                     lepton_P_gen_E = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5);
                     lepton_P_gen_phi = genPart.phi();
                     lepton_P_gen_eta = genPart.eta();
                   }
                 } 
	       if ((genPart.pdgId() == 12) || (genPart.pdgId() == 14)) {
                  if (pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5) > neut_P_gen_pT){
		     neut_P_gen_px = genPart.px();
                     neut_P_gen_py = genPart.py();
                     neut_P_gen_pz = genPart.pz();
                     neut_P_gen_pT = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5);
                     neut_P_gen_E = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5);
                     neut_P_gen_phi = genPart.phi();
                     neut_P_gen_eta = genPart.eta();
                  }
	       }
               if ((genPart.pdgId() == -12) || (genPart.pdgId() == -14)) {
                  if (pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5) > neut_N_gen_pT){
                     neut_N_gen_px = genPart.px();
                     neut_N_gen_py = genPart.py();
                     neut_N_gen_pz = genPart.pz();
                     neut_N_gen_pT = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py(),0.5);
                     neut_N_gen_E = pow(genPart.px()*genPart.px()+genPart.py()*genPart.py()+genPart.pz()*genPart.pz(),0.5);
                     neut_N_gen_phi = genPart.phi();
                     neut_N_gen_eta = genPart.eta();
                  }
               }
	     }  
	   }// for (auto&& genPart : *(genColl.product())) 
	} // if (MC && !MC_Signal) 

        double L_A = 0.;
	double L_B = 0.;
        double L_C = 0.;
        double L_D = 0.;
        double L_E = 0.;
        double L_F = 0.;
        double xn(-999.);
        double yn(-999.);
        double xf(-999.);
        double yf(-999.);

	if (year2017){
        	L_B = 2.360904801;
        	L_C = 5.313012839 + 3.264135878;
        	L_D = 4.074723964;
        	L_E = 8.958810514;
        	L_F = 1.708478656 + 7.877903151 +3.632463163;
	}
	
	double weight_arm0 =1.;
	double weight_arm1 =1.;
	//if (DATA || MC_Signal){
	//std::cout << recoMultiRPProtons->size() << std::endl;

	//std::cout << "read pps info" << std::endl;
	for(size_t i = 0;i<recoMultiRPProtons->size();i++){
		if(recoMultiRPProtons->at(i).validFit()){
			CTPPSDetId rpId((*recoMultiRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
			//			int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
			int armId = rpId.arm();
			ProtCand_xi->push_back(recoMultiRPProtons->at(i).xi());
			//std::cout << recoMultiRPProtons->at(i).time() << std::endl;
			ProtCand_t->push_back(recoMultiRPProtons->at(i).time());
			ProtCand_ThX->push_back(recoMultiRPProtons->at(i).thetaX());
			ProtCand_ThY->push_back(recoMultiRPProtons->at(i).thetaY());
			ProtCand_rpid->push_back(-999);
			std::cout << "xirec: " << recoMultiRPProtons->at(i).xi() << std::endl;
			std::cout << "armId: " << armId << std::endl;
			ProtCand_arm->push_back(armId);
			ProtCand_ismultirp->push_back(1);
                        ProtCand_x->push_back(-999.);
                        ProtCand_y->push_back(-999.);
	                int cont = 0;
                        for (const auto &tr : recoMultiRPProtons->at(i).contributingLocalTracks()) {
                            if(cont==0){
                              ProtCand_xn->push_back(tr->getX());
                              ProtCand_yn->push_back(tr->getY());
			      xn = tr->getX();
			      yn = tr->getY();
                            }
                           if(cont==1){
                              ProtCand_xf->push_back(tr->getX());
                              ProtCand_yf->push_back(tr->getY());
                              xf = tr->getX();
                              yf = tr->getY();
                           }
                          cont++;
                       }	
	
		       if (preVFP || postVFP){
			  weight_multitrack_45_F_B = multitrack_45_F_B->GetBinContent(1);
			  weight_multitrack_45_F_C = multitrack_45_F_C->GetBinContent(1);
			  weight_multitrack_45_F_G = multitrack_45_F_G->GetBinContent(1);
			  weight_multitrack_45_F_H = multitrack_45_F_H->GetBinContent(1);
			  weight_multitrack_45_N_B = multitrack_45_N_B->GetBinContent(1);
			  weight_multitrack_45_N_C = multitrack_45_N_C->GetBinContent(1);
			  weight_multitrack_45_N_G = multitrack_45_N_G->GetBinContent(1);
			  weight_multitrack_45_N_H = multitrack_45_N_H->GetBinContent(1);
			  weight_multitrack_56_F_B = multitrack_56_F_B->GetBinContent(1);
			  weight_multitrack_56_F_C = multitrack_56_F_C->GetBinContent(1);
			  weight_multitrack_56_F_G = multitrack_56_F_G->GetBinContent(1);
			  weight_multitrack_56_F_H = multitrack_56_F_H->GetBinContent(1);
			  weight_multitrack_56_N_B = multitrack_56_N_B->GetBinContent(1);
			  weight_multitrack_56_N_C = multitrack_56_N_C->GetBinContent(1);
			  weight_multitrack_56_N_G = multitrack_56_N_G->GetBinContent(1);
			  if (armId==0 && recoMultiRPProtons->at(i).xi() < 0.111 && recoMultiRPProtons->at(i).xi()  > 0.04){ 
	w_raddam_45_Multi_B_N = RadDam_45_N_B->GetBinContent(RadDam_45_N_B->GetXaxis()->FindBin(xn),RadDam_45_N_B->GetYaxis()->FindBin(yn));
	w_raddam_45_Multi_C_N = RadDam_45_N_C->GetBinContent(RadDam_45_N_B->GetXaxis()->FindBin(xn),RadDam_45_N_B->GetYaxis()->FindBin(yn));
	w_raddam_45_Multi_G_N = RadDam_45_N_G->GetBinContent(RadDam_45_N_B->GetXaxis()->FindBin(xn),RadDam_45_N_B->GetYaxis()->FindBin(yn));
	w_raddam_45_Multi_B_F = RadDam_45_F_B->GetBinContent(RadDam_45_F_B->GetXaxis()->FindBin(xf),RadDam_45_F_B->GetYaxis()->FindBin(yf));
	w_raddam_45_Multi_C_F = RadDam_45_F_C->GetBinContent(RadDam_45_F_B->GetXaxis()->FindBin(xf),RadDam_45_F_B->GetYaxis()->FindBin(yf));
	w_raddam_45_Multi_G_F = RadDam_45_F_G->GetBinContent(RadDam_45_F_B->GetXaxis()->FindBin(xf),RadDam_45_F_B->GetYaxis()->FindBin(yf));
	w_raddam_45_Multi = ((4.55/9.79)*w_raddam_45_Multi_B_N+(1.59/9.79)*w_raddam_45_Multi_C_N+(3.65/9.79)*w_raddam_45_Multi_G_N)*((4.55/9.79)*w_raddam_45_Multi_B_F+(1.59/9.79)*w_raddam_45_Multi_C_F+(3.65/9.79)*w_raddam_45_Multi_G_F);
        weight_45_N = weight_multitrack_45_N_B*(4.55/9.79) + weight_multitrack_45_N_C*(1.59/9.79) + weight_multitrack_45_N_G*(3.65/9.79);
        weight_45_F = weight_multitrack_45_F_B*(4.55/9.79) + weight_multitrack_45_F_C*(1.59/9.79) + weight_multitrack_45_F_G*(3.65/9.79);
        weight_45 = (weight_45_N+weight_45_F)/2;
	weight_arm0 = w_raddam_45_Multi*weight_45;
			  }
                          else if (armId==1 && recoMultiRPProtons->at(i).xi() < 0.104 && recoMultiRPProtons->at(i).xi()  > 0.04){
        w_raddam_56_Multi_B_N = RadDam_56_N_B->GetBinContent(RadDam_56_N_B->GetXaxis()->FindBin(xn),RadDam_56_N_B->GetYaxis()->FindBin(yn));
        w_raddam_56_Multi_C_N = RadDam_56_N_C->GetBinContent(RadDam_56_N_B->GetXaxis()->FindBin(xn),RadDam_56_N_B->GetYaxis()->FindBin(yn));
        w_raddam_56_Multi_G_N = RadDam_56_N_G->GetBinContent(RadDam_56_N_B->GetXaxis()->FindBin(xn),RadDam_56_N_B->GetYaxis()->FindBin(yn));
        w_raddam_56_Multi_B_F = RadDam_56_F_B->GetBinContent(RadDam_56_F_B->GetXaxis()->FindBin(xf),RadDam_56_F_B->GetYaxis()->FindBin(yf));
        w_raddam_56_Multi_C_F = RadDam_56_F_C->GetBinContent(RadDam_56_F_B->GetXaxis()->FindBin(xf),RadDam_56_F_B->GetYaxis()->FindBin(yf));
        w_raddam_56_Multi_G_F = RadDam_56_F_G->GetBinContent(RadDam_56_F_B->GetXaxis()->FindBin(xf),RadDam_56_F_B->GetYaxis()->FindBin(yf));
        w_raddam_56_Multi = ((4.55/9.79)*w_raddam_56_Multi_B_N+(1.59/9.79)*w_raddam_56_Multi_C_N+(3.65/9.79)*w_raddam_56_Multi_G_N)*((4.55/9.79)*w_raddam_56_Multi_B_F+(1.59/9.79)*w_raddam_56_Multi_C_F+(3.65/9.79)*w_raddam_56_Multi_G_F);
        weight_56_N = weight_multitrack_56_N_B*(4.55/9.79) + weight_multitrack_56_N_C*(1.59/9.79) + weight_multitrack_56_N_G*(3.65/9.79);
        weight_56_F = weight_multitrack_56_F_B*(4.55/9.79) + weight_multitrack_56_F_C*(1.59/9.79) + weight_multitrack_56_F_G*(3.65/9.79);
        weight_56 = (weight_56_N+weight_56_F)/2;
        weight_arm1 = w_raddam_56_Multi*weight_56;
                          }
		       }

		       if (year2017){
           	          bool fid_45 = false;
               	          bool fid_56 = false;
          	          if(recoMultiRPProtons->at(i).xi() >= .02 &&
                             xn >= 1.995 &&
                             xn <= 24.334 &&
                             yn >= -10.098 &&
                             yn <= 4.298) fid_45 = true;
                          if(recoMultiRPProtons->at(i).xi()  >= .03 &&
                             xn >= 2.422 &&
                             xn <= 24.620 &&
                             yn >= -9.698 &&
                             yn <= 4.698) fid_56 = true;
			  if (fid_45 && armId == 0 && recoMultiRPProtons->at(i).thetaX() < .02 && recoMultiRPProtons->at(i).thetaY() < .02){
        			multiRP_weight_B = eff_multiRP_B_45->GetBinContent(eff_multiRP_B_45->FindBin(xf,yf));
        			multiRP_weight_C1 = eff_multiRP_C1_45->GetBinContent(eff_multiRP_C1_45->FindBin(xf, yf));
       		 		multiRP_weight_C2 = eff_multiRP_C2_45->GetBinContent(eff_multiRP_C2_45->FindBin(xf, yf));
        			multiRP_weight_D = eff_multiRP_D_45->GetBinContent(eff_multiRP_D_45->FindBin(xf, yf));
        			multiRP_weight_E = eff_multiRP_E_45->GetBinContent(eff_multiRP_E_45->FindBin(xf, yf));
        			multiRP_weight_F1 = eff_multiRP_F1_45->GetBinContent(eff_multiRP_F1_45->FindBin(xf, yf));
        			multiRP_weight_F2 = eff_multiRP_F2_45->GetBinContent(eff_multiRP_F2_45->FindBin(xf, yf));
        			multiRP_weight_F3 = eff_multiRP_F3_45->GetBinContent(eff_multiRP_F3_45->FindBin(xf, yf));
        			strips_weight_rddam_B = eff_strips_rddam_B_45->GetBinContent(eff_strips_rddam_B_45->FindBin(xn,yn));
        			strips_weight_rddam_C = eff_strips_rddam_C_45->GetBinContent(eff_strips_rddam_C_45->FindBin(xn,yn));
        			strips_weight_rddam_D = eff_strips_rddam_D_45->GetBinContent(eff_strips_rddam_D_45->FindBin(xn,yn));
        			strips_weight_rddam_E = eff_strips_rddam_E_45->GetBinContent(eff_strips_rddam_E_45->FindBin(xn,yn));
        			strips_weight_rddam_F = eff_strips_rddam_F_45->GetBinContent(eff_strips_rddam_F_45->FindBin(xn,yn));
        			strips_weight_multitrack_B=eff_strips_multitrack_B_45->GetBinContent(eff_strips_multitrack_B_45->FindBin(xn));
        			strips_weight_multitrack_C=eff_strips_multitrack_C_45->GetBinContent(eff_strips_multitrack_C_45->FindBin(xn));
        			strips_weight_multitrack_D=eff_strips_multitrack_D_45->GetBinContent(eff_strips_multitrack_D_45->FindBin(xn));
        			strips_weight_multitrack_E=eff_strips_multitrack_E_45->GetBinContent(eff_strips_multitrack_E_45->FindBin(xn));
        			strips_weight_multitrack_F=eff_strips_multitrack_F_45->GetBinContent(eff_strips_multitrack_F_45->FindBin(xn));
        			multiRP_weight = (multiRP_weight_B*L_B + multiRP_weight_C1*5.313012839 + multiRP_weight_C2*3.264135878 + multiRP_weight_D*L_D + multiRP_weight_E*L_E + multiRP_weight_F1*1.708478656 + multiRP_weight_F2*7.877903151 + multiRP_weight_F3*3.632463163)/(L_B+5.313012839+3.264135878+L_D+L_E+L_F);
        			strips_weight_rddam = (strips_weight_rddam_B*L_B + strips_weight_rddam_C*L_C + strips_weight_rddam_D*L_D + strips_weight_rddam_E*L_E + strips_weight_rddam_F*L_F)/(L_B+L_C+L_D+L_E+L_F);
        			strips_weight_multitrack = (strips_weight_multitrack_B*L_B + strips_weight_multitrack_C*L_C + strips_weight_multitrack_D*L_D+strips_weight_multitrack_E*L_E + strips_weight_multitrack_F*L_F)/(L_B+L_C+L_D+L_E+L_F);
        			weight_arm0 = multiRP_weight*strips_weight_rddam*strips_weight_multitrack;
		         }
			 else if (fid_56 && armId == 1 && recoMultiRPProtons->at(i).thetaX() < .02 && recoMultiRPProtons->at(i).thetaY() < .02){
                                multiRP_weight_B = eff_multiRP_B_56->GetBinContent(eff_multiRP_B_56->FindBin(xn,yn));
                                multiRP_weight_C1 = eff_multiRP_C1_56->GetBinContent(eff_multiRP_C1_56->FindBin(xn, yn));
                                multiRP_weight_C2 = eff_multiRP_C2_56->GetBinContent(eff_multiRP_C2_56->FindBin(xn, yn));
                                multiRP_weight_D = eff_multiRP_D_56->GetBinContent(eff_multiRP_D_56->FindBin(xn, yn));
                                multiRP_weight_E = eff_multiRP_E_56->GetBinContent(eff_multiRP_E_56->FindBin(xn, yn));
                                multiRP_weight_F1 = eff_multiRP_F1_56->GetBinContent(eff_multiRP_F1_56->FindBin(xn, yn));
                                multiRP_weight_F2 = eff_multiRP_F2_56->GetBinContent(eff_multiRP_F2_56->FindBin(xn, yn));
                                multiRP_weight_F3 = eff_multiRP_F3_56->GetBinContent(eff_multiRP_F3_56->FindBin(xn, yn));
                                strips_weight_rddam_B = eff_strips_rddam_B_56->GetBinContent(eff_strips_rddam_B_56->FindBin(xn,yn));
                                strips_weight_rddam_C = eff_strips_rddam_C_56->GetBinContent(eff_strips_rddam_C_56->FindBin(xn,yn));
                                strips_weight_rddam_D = eff_strips_rddam_D_56->GetBinContent(eff_strips_rddam_D_56->FindBin(xn,yn));
                                strips_weight_rddam_E = eff_strips_rddam_E_56->GetBinContent(eff_strips_rddam_E_56->FindBin(xn,yn));
                                strips_weight_rddam_F = eff_strips_rddam_F_56->GetBinContent(eff_strips_rddam_F_56->FindBin(xn,yn));
                                strips_weight_multitrack_B=eff_strips_multitrack_B_56->GetBinContent(eff_strips_multitrack_B_56->FindBin(xn));
                                strips_weight_multitrack_C=eff_strips_multitrack_C_56->GetBinContent(eff_strips_multitrack_C_56->FindBin(xn));
                                strips_weight_multitrack_D=eff_strips_multitrack_D_56->GetBinContent(eff_strips_multitrack_D_56->FindBin(xn));
                                strips_weight_multitrack_E=eff_strips_multitrack_E_56->GetBinContent(eff_strips_multitrack_E_56->FindBin(xn));
                                strips_weight_multitrack_F=eff_strips_multitrack_F_56->GetBinContent(eff_strips_multitrack_F_56->FindBin(xn));
                                multiRP_weight = (multiRP_weight_B*L_B + multiRP_weight_C1*5.313012839 + multiRP_weight_C2*3.264135878 + multiRP_weight_D*L_D + multiRP_weight_E*L_E + multiRP_weight_F1*1.708478656 + multiRP_weight_F2*7.877903151 + multiRP_weight_F3*3.632463163)/(L_B+5.313012839+3.264135878+L_D+L_E+L_F);
                                strips_weight_rddam = (strips_weight_rddam_B*L_B + strips_weight_rddam_C*L_C + strips_weight_rddam_D*L_D + strips_weight_rddam_E*L_E + strips_weight_rddam_F*L_F)/(L_B+L_C+L_D+L_E+L_F);
                                strips_weight_multitrack = (strips_weight_multitrack_B*L_B + strips_weight_multitrack_C*L_C + strips_weight_multitrack_D*L_D+strips_weight_multitrack_E*L_E + strips_weight_multitrack_F*L_F)/(L_B+L_C+L_D+L_E+L_F);
                                std::cout << multiRP_weight*strips_weight_rddam*strips_weight_multitrack << std::endl;
				weight_arm1 = multiRP_weight*strips_weight_rddam*strips_weight_multitrack;
			 }
		     }	
		     else if (year2018){
		   	     if (armId == 0){
				multiRP2018_weight_A  = eff_multiRP_A_45_2018->GetBinContent(eff_multiRP_A_45_2018->FindBin(xf,yf));	
                                multiRP2018_weight_B = eff_multiRP_B_45_2018->GetBinContent(eff_multiRP_B_45_2018->FindBin(xf,yf));
                                multiRP2018_weight_B1 = eff_multiRP_B1_45_2018->GetBinContent(eff_multiRP_B1_45_2018->FindBin(xf,yf));
                                multiRP2018_weight_B2 = eff_multiRP_B2_45_2018->GetBinContent(eff_multiRP_B2_45_2018->FindBin(xf,yf));
                                multiRP2018_weight_C  = eff_multiRP_C_45_2018->GetBinContent(eff_multiRP_C_45_2018->FindBin(xf,yf));	
                                multiRP2018_weight_D = eff_multiRP_D_45_2018->GetBinContent(eff_multiRP_D_45_2018->FindBin(xf,yf));
			        multiRP2018_weight_D1 = eff_multiRP_D1_45_2018->GetBinContent(eff_multiRP_D1_45_2018->FindBin(xf,yf));
                                multiRP2018_weight_D2 = eff_multiRP_D2_45_2018->GetBinContent(eff_multiRP_D2_45_2018->FindBin(xf,yf));
			     
                                radiation2018_weight_A_N  = eff_radiation_A_45_N->GetBinContent(eff_radiation_A_45_N->FindBin(xn,yn));
                                radiation2018_weight_B_N  = eff_radiation_B_45_N->GetBinContent(eff_radiation_B_45_N->FindBin(xn,yn));
                                radiation2018_weight_B1_N = eff_radiation_B1_45_N->GetBinContent(eff_radiation_B1_45_N->FindBin(xn,yn));
                                radiation2018_weight_B2_N = eff_radiation_B2_45_N->GetBinContent(eff_radiation_B2_45_N->FindBin(xn,yn));
                                radiation2018_weight_C_N  = eff_radiation_C_45_N->GetBinContent(eff_radiation_C_45_N->FindBin(xn,yn));
                                radiation2018_weight_D_N  = eff_radiation_D_45_N->GetBinContent(eff_radiation_D_45_N->FindBin(xn,yn));
                                radiation2018_weight_D1_N = eff_radiation_D1_45_N->GetBinContent(eff_radiation_D1_45_N->FindBin(xn,yn));
                                radiation2018_weight_D2_N = eff_radiation_D2_45_N->GetBinContent(eff_radiation_D2_45_N->FindBin(xn,yn));

                                radiation2018_weight_A_F  = eff_radiation_A_45_F->GetBinContent(eff_radiation_A_45_F->FindBin(xf,yf));
                                radiation2018_weight_B_F = eff_radiation_B_45_F->GetBinContent(eff_radiation_B_45_F->FindBin(xf,yf));
                                radiation2018_weight_B1_F = eff_radiation_B1_45_F->GetBinContent(eff_radiation_B1_45_F->FindBin(xf,yf));
                                radiation2018_weight_B2_F = eff_radiation_B2_45_F->GetBinContent(eff_radiation_B2_45_F->FindBin(xf,yf));
                                radiation2018_weight_C_F  = eff_radiation_C_45_F->GetBinContent(eff_radiation_C_45_F->FindBin(xf,yf));
                                radiation2018_weight_D_F = eff_radiation_D_45_F->GetBinContent(eff_radiation_D_45_F->FindBin(xf,yf));
                                radiation2018_weight_D1_F = eff_radiation_D1_45_F->GetBinContent(eff_radiation_D1_45_F->FindBin(xf,yf));

				weight2018_rddam_N = (radiation2018_weight_A_N+radiation2018_weight_B_N+radiation2018_weight_C_N+radiation2018_weight_D_N)/(L_A+L_B+L_C+L_E+L_D);
                                weight2018_rddam_F = (radiation2018_weight_A_F+radiation2018_weight_B_F+radiation2018_weight_C_F+radiation2018_weight_D_F)/(L_A+L_B+L_C+L_E+L_D);
				weight2018_rddam = (weight2018_rddam_F+weight2018_rddam_N)/2.;
				weight_arm0 = multiRP2018_weight*weight2018_rddam;
			     }
			     else if (armId ==1){
                                multiRP2018_weight_A  = eff_multiRP_A_56_2018->GetBinContent(eff_multiRP_A_56_2018->FindBin(xf,yf));   
                                multiRP2018_weight_B = eff_multiRP_B_56_2018->GetBinContent(eff_multiRP_B_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_B1 = eff_multiRP_B1_56_2018->GetBinContent(eff_multiRP_B1_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_B2 = eff_multiRP_B2_56_2018->GetBinContent(eff_multiRP_B2_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_C  = eff_multiRP_C_56_2018->GetBinContent(eff_multiRP_C_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_D = eff_multiRP_D_56_2018->GetBinContent(eff_multiRP_D_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_D1 = eff_multiRP_D1_56_2018->GetBinContent(eff_multiRP_D1_56_2018->FindBin(xf,yf));
                                multiRP2018_weight_D2 = eff_multiRP_D2_56_2018->GetBinContent(eff_multiRP_D2_56_2018->FindBin(xf,yf));

                                radiation2018_weight_A_N  = eff_radiation_A_56_N->GetBinContent(eff_radiation_A_56_N->FindBin(xn,yn));
                                radiation2018_weight_B_N  = eff_radiation_B_56_N->GetBinContent(eff_radiation_B_56_N->FindBin(xn,yn));
                                radiation2018_weight_B1_N = eff_radiation_B1_56_N->GetBinContent(eff_radiation_B1_56_N->FindBin(xn,yn));
                                radiation2018_weight_B2_N = eff_radiation_B2_56_N->GetBinContent(eff_radiation_B2_56_N->FindBin(xn,yn));
                                radiation2018_weight_C_N  = eff_radiation_C_56_N->GetBinContent(eff_radiation_C_56_N->FindBin(xn,yn));
                                radiation2018_weight_D_N  = eff_radiation_D_56_N->GetBinContent(eff_radiation_D_56_N->FindBin(xn,yn));
                                radiation2018_weight_D1_N = eff_radiation_D1_56_N->GetBinContent(eff_radiation_D1_56_N->FindBin(xn,yn));
                                radiation2018_weight_D2_N = eff_radiation_D2_56_N->GetBinContent(eff_radiation_D2_56_N->FindBin(xn,yn));

                                radiation2018_weight_A_F  = eff_radiation_A_56_F->GetBinContent(eff_radiation_A_56_F->FindBin(xf,yf));
                                radiation2018_weight_B_F = eff_radiation_B_56_F->GetBinContent(eff_radiation_B_56_F->FindBin(xf,yf));
                                radiation2018_weight_B1_F = eff_radiation_B1_56_F->GetBinContent(eff_radiation_B1_56_F->FindBin(xf,yf));
                                radiation2018_weight_B2_F = eff_radiation_B2_56_F->GetBinContent(eff_radiation_B2_56_F->FindBin(xf,yf));
                                radiation2018_weight_C_F  = eff_radiation_C_56_F->GetBinContent(eff_radiation_C_56_F->FindBin(xf,yf));
                                radiation2018_weight_D_F = eff_radiation_D_56_F->GetBinContent(eff_radiation_D_56_F->FindBin(xf,yf));
                                radiation2018_weight_D1_F = eff_radiation_D1_56_F->GetBinContent(eff_radiation_D1_56_F->FindBin(xf,yf));
                                radiation2018_weight_D2_F = eff_radiation_D2_56_F->GetBinContent(eff_radiation_D2_56_F->FindBin(xf,yf));
                                multiRP2018_weight = (multiRP2018_weight_A*L_A+multiRP2018_weight_B*L_B + multiRP2018_weight_C*L_C + multiRP2018_weight_D*L_D)/(L_A+L_B+L_C+L_D);
                                weight2018_rddam_N = (radiation2018_weight_A_N+radiation2018_weight_B_N+radiation2018_weight_C_N+radiation2018_weight_D_N)/(L_A+L_B+L_C+L_E+L_D);
                                weight2018_rddam_F = (radiation2018_weight_A_F+radiation2018_weight_B_F+radiation2018_weight_C_F+radiation2018_weight_D_F)/(L_A+L_B+L_C+L_E+L_D);
                                weight2018_rddam = (weight2018_rddam_F+weight2018_rddam_N)/2.;
                                weight_arm1 = multiRP2018_weight*weight2018_rddam;

			     }
		     }   
	      }
	}

	ProtCand_weight = weight_arm0*weight_arm1;
	if (ProtCand_xi->size()<2){
		ProtCand_xi->push_back(-1);
		ProtCand_xi->push_back(-1);
		ProtCand_ismultirp->push_back(0);
		ProtCand_ismultirp->push_back(0);
	}
	for(size_t i = 0;i<recoSingleRPProtons->size();i++){
		if(recoSingleRPProtons->at(i).validFit()){
			CTPPSDetId rpId((*recoSingleRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
			int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
			//                      int armId = rpId.arm();
			//ProtCand_xi->push_back(recoSingleRPProtons->at(i).xi());
			//                      ProtCand_t->push_back(recoSingleRPProtons->at(i).t());
			ProtCand_ThX->push_back(recoSingleRPProtons->at(i).thetaX());
			ProtCand_ThY->push_back(recoSingleRPProtons->at(i).thetaY());
			ProtCand_rpid->push_back(decRPId);
			ProtCand_arm->push_back(-999);
			//ProtCand_ismultirp->push_back(0);
                        for (const auto &tr : recoSingleRPProtons->at(i).contributingLocalTracks()) {
                            ProtCand_x->push_back(tr->getX());
                            ProtCand_y->push_back(tr->getY());
                        }
                        ProtCand_xn->push_back(-999.);
                        ProtCand_yn->push_back(-999.);
                        ProtCand_xf->push_back(-999.);
                        ProtCand_yf->push_back(-999.);
			//ProtCand_weight->push_back(1.);
		}
	}
	//}
	

	/*
	   if (MC && !MC_Signal){
	   double strips_xi_arm0_N;
	   double strips_xi_arm0_F;
	   double strips_xi_arm1_N;
	   double strips_xi_arm1_F;
	   float a, b, c, d;


	   Int_t ncols;
	   Int_t nlines = 0;
	   FILE *fp = fopen("RPs_xi.txt","r");

	   int rand_no = rand() % 557673;

	   while (nlines<rand_no+1) {

	   ncols = fscanf(fp,"%f %f %f %f",&a, &b, &c, &d);


	   if (ncols < 0) break;

	   if (nlines==rand_no){
	   strips_xi_arm0_N = a;
	   strips_xi_arm0_F = b;
	   strips_xi_arm1_N = c;
	   strips_xi_arm1_F = d;
	   }

	   nlines++;

	   }

	   fclose(fp);
	   if (strips_xi_arm0_N>0.005) {ArmF_N_xi->push_back(strips_xi_arm0_N);}
	   if (strips_xi_arm0_F>0.005) {ArmF_F_xi->push_back(strips_xi_arm0_F);}
	   if (strips_xi_arm1_N>0.005) {ArmB_N_xi->push_back(strips_xi_arm1_N);}
	   if (strips_xi_arm1_F>0.005) {ArmB_F_xi->push_back(strips_xi_arm1_F);}

	   }
	   */

	TRandom3 rnd(1234);
	RoccoR rc;
	//RoccoR rcA("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016aUL.txt"); 
	//RoccoR rcB("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016bUL.txt"); 

	//rcA.init(edm::FileInPath("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/plugins/RoccoR2016aUL.txt").fullPath());
        //rcB.init(edm::FileInPath("/afs/cern.ch/user/j/jgomespi/private/workspace/master_thesis/FullLep_Studies/2016_analysis/Ntuple_build/CMSSW_10_6_20/src/MakeNTuple/MakeNTuple/RoccoR2016bUL.txt").fullPath());

        for(size_t i = 0;i<electrons->size();i++){
             if (electrons->at(i).pt()>20.
             && abs(electrons->at(i).eta())<2.4
             && electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V2-loose")){
                if (electrons->at(i).charge() == 1){
                        if (electrons->at(i).pt()>lepton_P_pT){
                                lepton_P_pT = electrons->at(i).pt();
                                lepton_P_index = i;
                                lepton_P_ismuon = 0;
                                lepton_P_iselectron = 1;
				//std::cout << "electron_P_pT inside loop = electrons->at(" << i << ")" << " = " << lepton_P_pT << std::endl;
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
           if(muons->at(i).pt()>20. &&
             abs(muons->at(i).eta())<2.4 &&
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
	//if (lepton_P_index==-1 || lepton_N_index==-1){continue;}
        std::cout << "lepton_P_pT: " << lepton_P_pT << std::endl;
        std::cout << "lepton_N_pT: " << lepton_N_pT << std::endl;

        if (preVFP)  rc.init("RoccoR2016aUL.txt");
        if (postVFP) rc.init("RoccoR2016bUL.txt");
        if (year2017)    rc.init("RoccoR2017UL.txt");
        if (year2018)    rc.init("RoccoR2018UL.txt");

        if (lepton_N_ismuon==1){
                float bestdr(9999.);
                float bestmatch_pt(-1);

                if (MC && MC_Signal) {
                   iEvent.getByToken( mcEventToken, EvtHandle );
                   const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
                   for(HepMC::GenEvent::particle_const_iterator j=Evt->particles_begin(); j != Evt->particles_end();++j) {
                      if (std::abs((*j)->pdg_id()) != 13 || (*j)->status() != 1) continue;
                      HepMC::FourVector momentum=(*j)->momentum();
                      const HepMC::FourVector p((*j)->momentum());
                      double px = momentum.x(),  py = momentum.y(), pz = momentum.z(), e = sqrt(px*px+py*py+pz*pz);
                      TLorentzVector* GenPart = new TLorentzVector(px,py,pz,e);
                      float dr = deltaR2(muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(), GenPart->Eta(),GenPart->Phi());
                      if (dr < 0.1 && dr < bestdr) {
                         bestmatch_pt = GenPart->Pt();
                         bestdr = dr;
                      }
                    }
                 }


                 if (MC && !MC_Signal) {
                    edm::Handle<reco::GenParticleCollection> genColl;
                    iEvent.getByToken(genToken, genColl);
                    for (auto&& genPart : *(genColl.product())) {
                        if (genPart.status()!=1 || std::abs(genPart.pdgId()) != 13) continue;
                        float dr = deltaR2(muons->at(lepton_N_index).eta(),muons->at(lepton_N_index).phi(), genPart.eta(),genPart.phi());
                        if (dr < 0.1 && dr < bestdr) {
                           bestmatch_pt = genPart.pt();
                           bestdr = dr;
                        }
                    }
                 }
                 if (MC && preVFP){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_N_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),muons->at(lepton_N_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && postVFP){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_N_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),muons->at(lepton_N_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && year2017){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_N_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),muons->at(lepton_N_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && year2018){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_N_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi(),muons->at(lepton_N_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (DATA){muSF = rc.kScaleDT(muons->at(lepton_N_index).charge(), muons->at(lepton_N_index).pt(), muons->at(lepton_N_index).eta(), muons->at(lepton_N_index).phi());}

                 std::cout << "muons->at(lepton_N_index).pt(): " << muons->at(lepton_N_index).pt() << std::endl;
                 lepton_N_pT  = muons->at(lepton_N_index).pt()*muSF;
                 lepton_N_phi = muons->at(lepton_N_index).phi()*muSF;
                 lepton_N_eta = muons->at(lepton_N_index).eta()*muSF;
                 lepton_N_E   = muons->at(lepton_N_index).energy()*muSF;
                 lepton_N_px  = muons->at(lepton_N_index).px()*muSF;
                 lepton_N_py  = muons->at(lepton_N_index).py()*muSF;
                 lepton_N_pz  = muons->at(lepton_N_index).pz()*muSF;
                 lepton_N_vtxZ  = muons->at(lepton_N_index).vertex().z();
                 lepton_N_charge  = muons->at(lepton_N_index).charge();
                 lepton_N_PFBasedIso  = (muons->at(lepton_N_index).pfIsolationR04().sumChargedHadronPt + max(0., muons->at(lepton_N_index).pfIsolationR04().sumNeutralHadronEt + muons->at(lepton_N_index).pfIsolationR04().sumPhotonEt - 0.5*muons->at(lepton_N_index).pfIsolationR04().sumPUPt))/lepton_N_pT;
                 lepton_N_isTight  = muon::isTightMuon(muons->at(lepton_N_index), vertices->at(0));
                 lepton_N_isMedium  = muon::isMediumMuon(muons->at(lepton_N_index));
        }
        if (lepton_N_iselectron==1){
                 std::cout << "electrons->at(lepton_N_index).pt(): " << electrons->at(lepton_N_index).pt() << std::endl;
                 elSF = electrons->at(lepton_N_index).userFloat("ecalTrkEnergyPostCorr")/electrons->at(lepton_N_index).userFloat("ecalTrkEnergyPreCorr");
		 lepton_N_pT  = electrons->at(lepton_N_index).pt()*elSF;
                 lepton_N_phi = electrons->at(lepton_N_index).phi()*elSF;
                 lepton_N_eta = electrons->at(lepton_N_index).eta()*elSF;
                 lepton_N_E   = electrons->at(lepton_N_index).energy()*elSF;
                 lepton_N_px  = electrons->at(lepton_N_index).px()*elSF;
                 lepton_N_py  = electrons->at(lepton_N_index).py()*elSF;
                 lepton_N_pz  = electrons->at(lepton_N_index).pz()*elSF;
                 lepton_N_vtxZ  = electrons->at(lepton_N_index).vertex().z();
                 lepton_N_charge  = electrons->at(lepton_N_index).charge();
		 GsfElectron::PflowIsolationVariables pfIso = electrons->at(lepton_N_index).pfIsolationVariables();
                 double isoChargedHadrons = pfIso.sumChargedHadronPt;
                 double isoNeutralHadrons = pfIso.sumNeutralHadronEt;
                 double isoPhotons = pfIso.sumPhotonEt;
                 // The following line provides the isolation sum for particles from PU needed if deltaBeta pile-up correction is used
                 double isoChargedHadFromPileup = pfIso.sumPUPt;
                 lepton_N_PFBasedIso = (isoChargedHadrons + max(0., isoNeutralHadrons + isoPhotons - isoChargedHadFromPileup))/lepton_N_pT;
                 lepton_N_isTight = electrons->at(lepton_N_index).electronID("cutBasedElectronID-Fall17-94X-V2-tight");
                 lepton_N_isMedium = electrons->at(lepton_N_index).electronID("cutBasedElectronID-Fall17-94X-V2-medium");
        }
        if (lepton_P_ismuon==1){
                float bestdr(9999.);
                float bestmatch_pt(-1);

                if (MC && MC_Signal) {
                   iEvent.getByToken( mcEventToken, EvtHandle );
                   const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
                   for(HepMC::GenEvent::particle_const_iterator j=Evt->particles_begin(); j != Evt->particles_end();++j) {
                      if (std::abs((*j)->pdg_id()) != 13 || (*j)->status() != 1) continue;
                      HepMC::FourVector momentum=(*j)->momentum();
                      const HepMC::FourVector p((*j)->momentum());
                      double px = momentum.x(),  py = momentum.y(), pz = momentum.z(), e = sqrt(px*px+py*py+pz*pz);
                      TLorentzVector* GenPart = new TLorentzVector(px,py,pz,e);
                      float dr = deltaR2(muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(), GenPart->Eta(),GenPart->Phi());
                      if (dr < 0.1 && dr < bestdr) {
                         bestmatch_pt = GenPart->Pt();
                         bestdr = dr;
                      }
                    }
                 }

                 if (MC && !MC_Signal) {
                    edm::Handle<reco::GenParticleCollection> genColl;
                    iEvent.getByToken(genToken, genColl);
                    for (auto&& genPart : *(genColl.product())) {
                        if (genPart.status()!=1 || std::abs(genPart.pdgId()) != 13) continue;
                        float dr = deltaR2(muons->at(lepton_P_index).eta(),muons->at(lepton_P_index).phi(), genPart.eta(),genPart.phi());
                        if (dr < 0.1 && dr < bestdr) {
                           bestmatch_pt = genPart.pt();
                           bestdr = dr;
                        }
                    }
                 }
                 if (MC && preVFP){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_P_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(), muons->at(lepton_P_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && postVFP){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_P_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(), muons->at(lepton_P_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && year2017){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_P_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(),muons->at(lepton_P_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (MC && year2018){
                    if (bestmatch_pt>0){muSF = rc.kSpreadMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(),bestmatch_pt);}
                    else if (muons->at(lepton_P_index).innerTrack().isNonnull()){muSF = rc.kSmearMC(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi(), muons->at(lepton_P_index).innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());}
                 }
                 if (DATA){muSF = rc.kScaleDT(muons->at(lepton_P_index).charge(), muons->at(lepton_P_index).pt(), muons->at(lepton_P_index).eta(), muons->at(lepton_P_index).phi());}
                 std::cout << "muons->at(lepton_P_index).pt(): " << muons->at(lepton_P_index).pt() << std::endl;
                 lepton_P_pT  = muons->at(lepton_P_index).pt()*muSF;
                 lepton_P_phi = muons->at(lepton_P_index).phi()*muSF;
                 lepton_P_eta = muons->at(lepton_P_index).eta()*muSF;
                 lepton_P_E   = muons->at(lepton_P_index).energy()*muSF;
                 lepton_P_px  = muons->at(lepton_P_index).px()*muSF;
                 lepton_P_py  = muons->at(lepton_P_index).py()*muSF;
                 lepton_P_pz  = muons->at(lepton_P_index).pz()*muSF;
                 lepton_P_vtxZ  = muons->at(lepton_P_index).vertex().z();
                 lepton_P_charge  = muons->at(lepton_P_index).charge();
                 lepton_P_PFBasedIso  = (muons->at(lepton_P_index).pfIsolationR04().sumChargedHadronPt + max(0., muons->at(lepton_P_index).pfIsolationR04().sumNeutralHadronEt + muons->at(lepton_P_index).pfIsolationR04().sumPhotonEt - 0.5*muons->at(lepton_P_index).pfIsolationR04().sumPUPt))/lepton_P_pT;
                 lepton_P_isTight  = muon::isTightMuon(muons->at(lepton_P_index), vertices->at(0));
                 lepton_P_isMedium  = muon::isMediumMuon(muons->at(lepton_P_index));
        }
       if (lepton_P_iselectron==1){
                 std::cout << "electrons->at(lepton_P_index).pt(): " << electrons->at(lepton_P_index).pt() << std::endl;
                 elSF = electrons->at(lepton_P_index).userFloat("ecalTrkEnergyPostCorr")/electrons->at(lepton_P_index).userFloat("ecalTrkEnergyPreCorr");
                 lepton_P_pT  = electrons->at(lepton_P_index).pt()*elSF;
                 lepton_P_phi = electrons->at(lepton_P_index).phi()*elSF;
                 lepton_P_eta = electrons->at(lepton_P_index).eta()*elSF;
                 lepton_P_E   = electrons->at(lepton_P_index).energy()*elSF;
                 lepton_P_px  = electrons->at(lepton_P_index).px()*elSF;
                 lepton_P_py  = electrons->at(lepton_P_index).py()*elSF;
                 lepton_P_pz  = electrons->at(lepton_P_index).pz()*elSF;
                 lepton_P_vtxZ  = electrons->at(lepton_P_index).vertex().z();
                 lepton_P_charge  = electrons->at(lepton_P_index).charge();
	 	 GsfElectron::PflowIsolationVariables pfIso = electrons->at(lepton_P_index).pfIsolationVariables();
                 double isoChargedHadrons = pfIso.sumChargedHadronPt;
                 double isoNeutralHadrons = pfIso.sumNeutralHadronEt;
                 double isoPhotons = pfIso.sumPhotonEt;
                 // The following line provides the isolation sum for particles from PU needed if deltaBeta pile-up correction is used
                 double isoChargedHadFromPileup = pfIso.sumPUPt;
                 lepton_P_PFBasedIso = (isoChargedHadrons + max(0., isoNeutralHadrons + isoPhotons - isoChargedHadFromPileup))/lepton_P_pT;
                 lepton_P_isTight = electrons->at(lepton_P_index).electronID("cutBasedElectronID-Fall17-94X-V2-tight");
                 lepton_P_isMedium = electrons->at(lepton_P_index).electronID("cutBasedElectronID-Fall17-94X-V2-medium");
        }
        //ROOT::Math::PtEtaPhiEVector lepton_P(lepton_P_pT,lepton_P_eta,lepton_P_phi,lepton_P_E);
        //ROOT::Math::PtEtaPhiEVector lepton_N(lepton_N_pT,lepton_N_eta,lepton_N_phi,lepton_N_E);

 	if (lepton_P_ismuon==1 && lepton_N_ismuon==1){channel = 0;}
        if (lepton_P_iselectron==1 && lepton_N_iselectron==1){channel = 1;}
        if (lepton_P_ismuon==1 && lepton_N_iselectron==1){channel = 2;}
        if (lepton_P_iselectron==1 && lepton_N_ismuon==1){channel = 2;}
	std::cout << "channel: " << channel << std::endl;

	// Scale Factor for leptons
        weight1_id  = 1.;
	weight1_iso = 1.;
	weight2_id  = 1.;
	weight2_iso = 1.;
	
	if (channel == 0){
	   if (tag=="MC_preVFP" || tag=="MC_signal_preVFP"){
		weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ID_preVFP);
		weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ID_preVFP);		
                weight1_iso = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ISO_preVFP);
                weight2_iso = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ISO_preVFP);
	   }
	   else if (tag=="MC_postVFP" || tag=="MC_signal_postVFP"){
	 	weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ID_postVFP);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ID_postVFP);
                weight1_iso = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ISO_postVFP);
                weight2_iso = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ISO_postVFP);

	   }
           if (tag=="MC_2017" || tag=="MC_signal_2017"){
		weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ID_2017);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ID_2017);
                weight1_iso = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ISO_2017);
                weight2_iso = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ISO_2017);

           }
           else if (tag=="MC_2018" || tag=="MC_signal_2018"){
		weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ID_2018);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ID_2018);
                weight1_iso = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFmuon_ISO_2018);
                weight2_iso = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFmuon_ISO_2018);
           }
	}

        if (channel == 1){
           if (tag=="MC_preVFP" || tag=="MC_signal_preVFP"){
                weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFelectron_ID_preVFP);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFelectron_ID_preVFP);
           }
           else if (tag=="MC_postVFP" || tag=="MC_signal_postVFP"){
                weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFelectron_ID_postVFP);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFelectron_ID_postVFP);
           }
           if (tag=="MC_2017" || tag=="MC_signal_2017"){
                weight1_id  = 0.991*GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFelectron_ID_2017);
                weight2_id  = 0.991*GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFelectron_ID_2017);
           }
           else if (tag=="MC_2018" || tag=="MC_signal_2018"){
                weight1_id  = GetWeightLeptons(lepton_P_pT,lepton_P_eta,h_SFelectron_ID_2018);
                weight2_id  = GetWeightLeptons(lepton_N_pT,lepton_N_eta,h_SFelectron_ID_2018);
           }
        }

        if (channel == 2){
           
	   if (lepton_P_ismuon==1){
		muon_pT = lepton_P_pT; 
	   	electron_pT = lepton_N_pT;
		muon_eta = lepton_P_eta;
                electron_eta = lepton_N_eta;
	   }
	   else {
	        muon_pT = lepton_N_pT;   
                electron_pT = lepton_P_pT;
                muon_eta = lepton_N_eta;
                electron_eta = lepton_P_eta;	
	   }
	   if (tag=="MC_preVFP" || tag=="MC_signal_preVFP"){
                weight1_id  = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ID_preVFP);
                weight2_id  = GetWeightLeptons(electron_pT,electron_eta,h_SFelectron_ID_preVFP);
                weight1_iso = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ISO_preVFP);
           }
           else if (tag=="MC_postVFP" || tag=="MC_signal_postVFP"){
                weight1_id  = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ID_postVFP);
                weight2_id  = GetWeightLeptons(electron_pT,electron_eta,h_SFelectron_ID_postVFP);
                weight1_iso = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ISO_postVFP);

           }
           if (tag=="MC_2017" || tag=="MC_signal_2017"){
                weight1_id  = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ID_2017);
                weight2_id  = 0.991*GetWeightLeptons(electron_pT,electron_eta,h_SFelectron_ID_2017);
                weight1_iso = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ISO_2017);
           }
           else if (tag=="MC_2018" || tag=="MC_signal_2018"){
                weight1_id  = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ID_2018);
                weight2_id  = GetWeightLeptons(electron_pT,electron_eta,h_SFelectron_ID_2018);
                weight1_iso = GetWeightLeptons(muon_pT,muon_eta,h_SFmuon_ISO_2018);
           }
        }

	weight_lep = weight1_id*weight1_iso*weight2_id*weight2_iso;
	weight = weight_lep;
	
	METPx = (MET->front()).px();
        METPy = (MET->front()).py();
        METPt = (MET->front()).pt();
        METphi = (MET->front()).phi();

	int npfVtx = 0;
	for (size_t pf = 0; pf < PFCand->size(); pf++) {
		if (abs(PFCand->at(pf).pdgId()) == 211 || abs(PFCand->at(pf).pdgId()) == 11 || abs(PFCand->at(pf).pdgId()) == 13){
//			if (PFCand->at(pf).fromPV(0)>1) {
//				if (abs(PFCand->at(pf).dz())<0.15){
					npfVtx++;
					pfphi->push_back(PFCand->at(pf).phiAtVtx());
					pfeta->push_back(PFCand->at(pf).eta());
					pffromPV->push_back(PFCand->at(pf).fromPV(0));
					pfdz->push_back(PFCand->at(pf).dz());
					pfpt->push_back(PFCand->at(pf).pt());
//				}
//			}
		}
	}
/*	f_SFmuon_ID_preVFP->Close();
        f_SFmuon_ID_postVFP->Close();
        f_SFmuon_ID_2017->Close();
        f_SFmuon_ID_2018->Close();
        f_SFmuon_ISO_preVFP->Close();
        f_SFmuon_ISO_postVFP->Close();
        f_SFmuon_ISO_2017->Close();
        f_SFmuon_ISO_2018->Close();
        f_SFelectron_ID_preVFP->Close();
        f_SFelectron_ID_postVFP->Close();
        f_SFelectron_ID_2017->Close();
        f_SFelectron_ID_2018->Close();
//	eff_strips->Close();
//	eff_multiRP->Close();
	eff_radiation->Close();
        f_SFmuon_ID_preVFP->Delete();
        f_SFmuon_ID_postVFP->Delete();
        f_SFmuon_ID_2017->Delete();
        f_SFmuon_ID_2018->Delete();
        f_SFmuon_ISO_preVFP->Delete();
        f_SFmuon_ISO_postVFP->Delete();
        f_SFmuon_ISO_2017->Delete();
        f_SFmuon_ISO_2018->Delete();
        f_SFelectron_ID_preVFP->Delete();
        f_SFelectron_ID_postVFP->Delete();
        f_SFelectron_ID_2017->Delete();
        f_SFelectron_ID_2018->Delete();
//        eff_strips->Delete();
//        eff_multiRP->Delete();
        eff_radiation->Delete();
*/
	//delete rc;
	EventBranchs->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
MakeNTuple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
MakeNTuple::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeNTuple);
