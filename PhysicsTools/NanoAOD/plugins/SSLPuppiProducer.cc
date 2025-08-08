// Original Author:  Yuji Li
//         Created:  Wed, 13 Mar 2024 09:19:58 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
//
// class declaration
//
const double pi = 3.14159265358979323846;
class SSLPuppiProducer : public TritonEDProducer<> {
public:
  explicit SSLPuppiProducer(const edm::ParameterSet&);
  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override;
  ~SSLPuppiProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float deltaRcut = 0.8;
  float jetRadius_ = 0.8;
  int64_t npf;
  TH1D* h_gnn;
  TH1D* h_pf;
  TH1D* h_puppi;
  TH1D* h_pf_ak8;
  TH1D* h_puppi_ak8;
  TH1D* h_sd_ak8;
  TH1D* h_chs;

  TH1D* h_gnn_pt;
  TH1D* h_pf_pt;
  TH1D* h_puppi_pt;
  TH1D* h_chs_pt;

  TH1D* h_gnn_pt_diff;
  TH1D* h_pf_pt_diff;
  TH1D* h_puppi_pt_diff;
  TH1D* h_pf_pt_diff_ak8;
  TH1D* h_puppi_pt_diff_ak8;
  TH1D* h_sd_pt_diff_ak8;
  TH1D* h_chs_pt_diff;

  TH1D* h_gnn_mass;
  TH1D* h_pf_mass;
  TH1D* h_puppi_mass;
  TH1D* h_puppi_mass_ak8;
  TH1D* h_sd_mass_ak8;
  TH1D* h_chs_mass;


  TH1D* h_LV_score;
  TH1D* h_PU_score;
  TH1D* h_NE_score;
  TH2D* h2gnn; 
  TH2D* h2pf; 
  TH2D* h2puppi; 
  TH2D* h2chs;
  std::vector<float> mass_reso_gnn, mass_reso_pf, mass_reso_puppi, mass_reso_chs,pt_reso_gnn, pt_reso_pf, pt_reso_chs,pt_reso_puppi;
  std::vector<float> mass_truth_gnn_evt, pt_truth_gnn_evt,mass_truth_puppi_evt, mass_truth_chs_evt,pt_truth_puppi_evt,mass_truth_pf_evt,pt_truth_chs_evt, pt_truth_pf_evt;
  std::vector<float> mass_reso_gnn_evt, mass_reso_pf_evt, mass_reso_puppi_evt,mass_reso_chs_evt, pt_reso_gnn_evt, pt_reso_pf_evt,pt_reso_chs_evt,  pt_reso_puppi_evt;
  //TTree* tree_;
private:  
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genParticleSrc_;
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pf_token_;
  edm::EDGetTokenT<std::vector<pat::Jet>> ak8JetsToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>> GENJetsToken_;  
  //edm::EDGetTokenT<std::vector<pat::Jet>> ak8SDJetsToken_;
  //edm::EDGetTokenT<std::vector<reco::GenJet>> GENSDJetsToken_;
  const unsigned int max_n_pf_;
  
  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
// static data member definitions
//

//
// constructors and destructor
//
float quantile(const std::vector<float>& input, float q) {
  std::vector<float> sorted = input;
  std::sort(sorted.begin(), sorted.end());
  float pos = q * (sorted.size() - 1);
  size_t lower = static_cast<size_t>(pos);
  size_t upper = lower + 1;
  float weight = pos - lower;
  if (upper >= sorted.size()) {
      return sorted[lower];
  }
  return sorted[lower] * (1 - weight) + sorted[upper] * weight;
}

float getResol(const std::vector<float>& input) {
  return (quantile(input, 0.84) - quantile(input, 0.16)) / 2;
}

float median(const std::vector<float>& input) {
  std::vector<float> sorted = input;
  std::sort(sorted.begin(), sorted.end());
  size_t size = sorted.size();
  if (size % 2 == 0) {
      return (sorted[size / 2 - 1] + sorted[size / 2]) / 2;
  } else {
      return sorted[size / 2];
  }
}


SSLPuppiProducer::SSLPuppiProducer(const edm::ParameterSet& cfg)
    :TritonEDProducer<>(cfg),
    genParticleSrc_(mayConsume<std::vector<pat::PackedGenParticle>>(cfg.getParameter<edm::InputTag>("genParticleSrc"))),
    pf_token_(consumes<std::vector<pat::PackedCandidate>>(cfg.getParameter<edm::InputTag>("pf_src"))),
    ak8JetsToken_(consumes<std::vector<pat::Jet>>(cfg.getParameter<edm::InputTag>("JetsAK8"))),
    GENJetsToken_(consumes<std::vector<reco::GenJet>>(cfg.getParameter<edm::InputTag>("GENJetsAK8"))),
    //ak8SDJetsToken_(consumes<std::vector<pat::Jet>>(cfg.getParameter<edm::InputTag>("JetsAK8SD"))),
    //GENSDJetsToken_(consumes<std::vector<reco::GenJet>>(cfg.getParameter<edm::InputTag>("GENSDJetsAK8"))),
    max_n_pf_(cfg.getParameter<unsigned int>("max_n_pf")) {
  //register your products
  /*edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("myTree", "Example Tree");
  tree_->Branch("mass_reso_gnn", &mass_reso_gnn_evt);
  tree_->Branch("mass_reso_pf", &mass_reso_pf_evt);
  tree_->Branch("mass_reso_puppi", &mass_reso_puppi_evt);
  tree_->Branch("mass_reso_chs", &mass_reso_chs_evt);
  tree_->Branch("mass_truth_gnn", &mass_truth_gnn_evt);
  tree_->Branch("mass_truth_pf", &mass_truth_pf_evt);
  tree_->Branch("mass_truth_puppi", &mass_truth_puppi_evt);
  tree_->Branch("mass_truth_chs", &mass_truth_chs_evt);
  tree_->Branch("pt_reso_gnn", &pt_reso_gnn_evt);
  tree_->Branch("pt_reso_pf", &pt_reso_pf_evt);
  tree_->Branch("pt_reso_puppi", &pt_reso_puppi_evt);
  tree_->Branch("pt_reso_chs", &pt_reso_chs_evt);
  tree_->Branch("pt_truth_gnn", &pt_truth_gnn_evt);
  tree_->Branch("pt_truth_pf", &pt_truth_pf_evt);
  tree_->Branch("pt_truth_puppi", &pt_truth_puppi_evt);
  tree_->Branch("pt_truth_chs", &pt_truth_chs_evt);*/
  /*produces<std::vector<float>>("SSLscore");
  produces<std::vector<float>>("pfeta");
  produces<std::vector<float>>("pfphi");
  produces<std::vector<float>>("pfpuppipt");
  produces<std::vector<float>>("geneta");
  produces<std::vector<float>>("genphi");
  produces<std::vector<float>>("genpt"); */
  produces<std::vector<float>>("massdiff");

  
  h_gnn = new TH1D("mass_diff_GNN", "Mass diff", 40, -1, 1);
  h_puppi = new TH1D("mass_diff_PUPPI", "Mass diff", 40, -1, 1);
  h_chs = new TH1D("mass_diff_CHS", "Mass diff", 40, -1, 1);
  h_pf = new TH1D("mass_diff_PF", "Mass diff", 40, -1, 1);
  h_pf_ak8 = new TH1D("mass_diff_PF_ak8", "Mass diff", 40, -1, 1);
  h_puppi_ak8 = new TH1D("mass_diff_PUPPI_ak8", "Mass diff", 40, -1, 1);
  h_sd_ak8 = new TH1D("mass_diff_sd_ak8", "Mass diff", 40, -1, 1);

  h_gnn_pt_diff = new TH1D("mass_diff_GNN", "Mass diff", 40, -1, 1);
  h_puppi_pt_diff = new TH1D("mass_diff_PUPPI", "Mass diff", 40, -1, 1);
  h_puppi_pt_diff_ak8 = new TH1D("mass_diff_PUPPI_ak8", "Mass diff", 40, -1, 1);
  h_chs_pt_diff = new TH1D("mass_diff_CHS", "Mass diff", 40, -1, 1);
  h_pf_pt_diff = new TH1D("mass_diff_PF", "Mass diff", 40, -1, 1);
  h_pf_pt_diff_ak8 = new TH1D("mass_diff_PF_ak8", "Mass diff", 40, -1, 1);
  h_sd_pt_diff_ak8 = new TH1D("mass_diff_sd_ak8", "Mass diff", 40, -1, 1);

  h_gnn_pt = new TH1D("Pt_GNN", "Pt ", 40, 0, 1200);
  h_puppi_pt = new TH1D("Pt_PUPPI", "Pt", 40, 0, 1200);
  h_chs_pt = new TH1D("Pt_CHS", "Pt", 40, 0, 1200);
  h_pf_pt = new TH1D("Pt_PF", "Pt", 40, 0, 1200);

  h_gnn_mass = new TH1D("mass_GNN", "Mass ", 40, 0, 400);
  h_puppi_mass = new TH1D("mass_PUPPI", "Mass", 40, 0, 400);
  h_puppi_mass_ak8 = new TH1D("mass_PUPPI_ak8", "Mass", 40, 0, 400);
  h_sd_mass_ak8 = new TH1D("mass_sd_ak8", "Mass", 40, 0, 400);
  h_chs_mass = new TH1D("mass_CHS", "Mass", 40, 0, 400);
  h_pf_mass = new TH1D("mass_PF", "Mass", 40, 0, 400);

  h_LV_score = new TH1D("LV SSL", "LV SSL", 40, 0, 1);
  h_PU_score = new TH1D("PU SSL", "PU SSL", 40, 0, 1);
  h_NE_score = new TH1D("NE SSL", "NE SSL", 40, 0, 1);

  h2gnn = new TH2D("h2gnn", "2D Histogram", 40, 10, 170, 20, 0, 5);
  h2puppi = new TH2D("h2gnn", "2D Histogram", 40, 10, 170, 20, 0, 5);
  h2pf = new TH2D("h2gnn", "2D Histogram", 40, 10, 170, 20, 0, 5);
  h2chs = new TH2D("h2gnn", "2D Histogram", 40, 10, 170, 20, 0, 5);
  mass_reso_gnn.clear();mass_reso_pf.clear();mass_reso_puppi.clear();mass_reso_chs.clear();
  pt_reso_gnn.clear();pt_reso_pf.clear();pt_reso_puppi.clear();pt_reso_chs.clear();
  /* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
}

SSLPuppiProducer::~SSLPuppiProducer() {
   /*TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
   h_pf->SetDirectory(0);
   h_puppi->SetDirectory(0);
   h_pf_ak8->SetDirectory(0);
   h_puppi_ak8->SetDirectory(0);
   h_sd_ak8->SetDirectory(0);
   h_gnn->SetDirectory(0);

   h_pf->SetLineColor(kRed);
   h_puppi->SetLineColor(kGreen);
   h_pf_ak8->SetLineColor(kPink);
   h_puppi_ak8->SetLineColor(kOrange);
   h_sd_ak8->SetLineColor(kBlack);
   h_gnn->SetLineColor(kBlue);

   h_pf->SetLineWidth(2);
   h_puppi_ak8->SetLineWidth(2);
   h_sd_ak8->SetLineWidth(2);
   h_pf_ak8->SetLineWidth(2);
   h_puppi->SetLineWidth(2);
   h_gnn->SetLineWidth(2);

   h_gnn->GetXaxis()->SetTitle("mass diff");
   h_gnn->GetYaxis()->SetTitle("A.U.");

   h_pf->Scale(1./h_pf->Integral(),"width");
   h_puppi->Scale(1./h_puppi->Integral(),"width");
   h_pf_ak8->Scale(1./h_pf_ak8->Integral(),"width");
   h_puppi_ak8->Scale(1./h_puppi_ak8->Integral(),"width");
   h_sd_ak8->Scale(1./h_sd_ak8->Integral(),"width");
   h_gnn->Scale(1./h_gnn->Integral(),"width");

   h_gnn->Draw("h");
   h_puppi->Draw("h,same");
   h_pf->Draw("h,same");
   h_puppi_ak8->Draw("h,same");
   //h_sd_ak8->Draw("h,same");
   //h_pf_ak8->Draw("h,same");

   TLegend* legend = new TLegend(0.1, 0.65, 0.45, 0.8);
   legend->AddEntry(h_pf,Form("PF, #mu= %.2f, #sigma= %.2f", h_pf->GetMean(),h_pf->GetRMS()), "l");
   legend->AddEntry(h_puppi, Form("PUPPI, #mu= %.2f, #sigma= %.2f", h_puppi->GetMean(),h_puppi->GetRMS()), "l");
   //legend->AddEntry(h_pf_ak8, "PF", "l");
   legend->AddEntry(h_puppi_ak8, Form("PUPPI_AK8, #mu= %.2f, #sigma= %.2f", h_puppi_ak8->GetMean(),h_puppi_ak8->GetRMS()), "l");
   //legend->AddEntry(h_sd_ak8, Form("softdrop_AK8, #mu= %.2f, #sigma= %.2f", h_sd_ak8->GetMean(),h_sd_ak8->GetRMS()), "l");
   legend->AddEntry(h_gnn, Form("SSL, #mu= %.2f, #sigma= %.2f", h_gnn->GetMean(),h_gnn->GetRMS()), "l");
   legend->Draw("same");
   
   canvas->SaveAs("hist_mass_diff.png");
   std::cout<<"PF Mean:"<<median(mass_reso_pf)<<" PF stdv:"<<getResol(mass_reso_pf)<<std::endl;
   std::cout<<"PUPPI Mean:"<<median(mass_reso_puppi)<<" PUPPI stdv:"<<getResol(mass_reso_puppi)<<std::endl;
   std::cout<<"SSL Mean:"<<median(mass_reso_gnn)<<" SSL stdv:"<<getResol(mass_reso_gnn)<<std::endl;

   std::cout<<"PF Mean:"<<h_pf->GetMean()<<" PF RMS:"<<h_pf->GetRMS()<<std::endl;
   std::cout<<"PUPPI Mean:"<<h_puppi->GetMean()<<" PUPPI RMS:"<<h_puppi->GetRMS()<<std::endl;
   std::cout<<"PF_ak8 Mean:"<<h_pf_ak8->GetMean()<<" PF_ak8 RMS:"<<h_pf_ak8->GetRMS()<<std::endl;
   std::cout<<"PUPPI_ak8 Mean:"<<h_puppi_ak8->GetMean()<<" PUPPI_ak8 RMS:"<<h_puppi_ak8->GetRMS()<<std::endl;
   std::cout<<"SSL Mean:"<<h_gnn->GetMean()<<" SSL RMS:"<<h_gnn->GetRMS()<<std::endl;
   delete canvas;
   delete h_gnn;
   delete h_pf;
   delete h_pf_ak8;
   delete h_puppi_ak8;
   delete h_puppi;

   TCanvas* canvasp = new TCanvas("canvasp", "Canvasp", 800, 600);
   h_pf_pt_diff->SetDirectory(0);
   h_puppi_pt_diff->SetDirectory(0);
   h_pf_pt_diff_ak8->SetDirectory(0);
   h_puppi_pt_diff_ak8->SetDirectory(0);
   h_gnn_pt_diff->SetDirectory(0);

   h_pf_pt_diff->SetLineColor(kRed);
   h_puppi_pt_diff->SetLineColor(kGreen);
   h_pf_pt_diff_ak8->SetLineColor(kPink);
   h_puppi_pt_diff_ak8->SetLineColor(kOrange);
   h_gnn_pt_diff->SetLineColor(kBlue);

   h_pf_pt_diff->SetLineWidth(2);
   h_puppi_pt_diff->SetLineWidth(2);
   h_pf_pt_diff_ak8->SetLineWidth(2);
   h_puppi_pt_diff_ak8->SetLineWidth(2);
   h_gnn_pt_diff->SetLineWidth(2);

   h_gnn_pt_diff->GetXaxis()->SetTitle("pT diff");
   h_gnn_pt_diff->GetYaxis()->SetTitle("A.U.");

   h_pf_pt_diff->Scale(1./h_pf_pt_diff->Integral(),"width");
   h_puppi_pt_diff->Scale(1./h_puppi_pt_diff->Integral(),"width");
   h_pf_pt_diff_ak8->Scale(1./h_pf_pt_diff_ak8->Integral(),"width");
   h_puppi_pt_diff_ak8->Scale(1./h_puppi_pt_diff_ak8->Integral(),"width");
   h_gnn_pt_diff->Scale(1./h_gnn_pt_diff->Integral(),"width");

   h_gnn_pt_diff->Draw("h");
   h_puppi_pt_diff->Draw("h,same");
   h_pf_pt_diff->Draw("h,same");
   h_puppi_pt_diff_ak8->Draw("h,same");
   //h_pf_pt_diff_ak8->Draw("h,same");
   

   TLegend* legendp = new TLegend(0.2, 0.65, 0.35, 0.8);
   legendp->AddEntry(h_pf_pt_diff, "PF", "l");
   legendp->AddEntry(h_puppi_pt_diff, "PUPPI", "l");
   //egendp->AddEntry(h_pf_pt_diff_ak8, "PF", "l");
   legendp->AddEntry(h_puppi_pt_diff_ak8, "PUPPI_AK8", "l");
   legendp->AddEntry(h_gnn_pt_diff, "SSL", "l");
   legendp->Draw("same");
   
   canvasp->SaveAs("hist_pt_diff.png");

   TCanvas* canvas1 = new TCanvas("canvas1", "Canvas1", 800, 600);
   h_pf_pt->SetDirectory(0);
   h_puppi_pt->SetDirectory(0);
   h_gnn_pt->SetDirectory(0);

   h_pf_pt->SetLineColor(kRed);
   h_puppi_pt->SetLineColor(kGreen);
   h_gnn_pt->SetLineColor(kBlue);

   h_pf_pt->SetLineWidth(2);
   h_puppi_pt->SetLineWidth(2);
   h_gnn_pt->SetLineWidth(2);

   h_gnn_pt->GetXaxis()->SetTitle("mass diff");
   h_gnn_pt->GetYaxis()->SetTitle("A.U.");

   h_pf_pt->Scale(1./h_pf_pt->Integral(),"width");
   h_puppi_pt->Scale(1./h_puppi_pt->Integral(),"width");
   h_gnn_pt->Scale(1./h_gnn_pt->Integral(),"width");

   h_puppi_pt->Draw("h");
   h_gnn_pt->Draw("h,same");
   h_pf_pt->Draw("h,same");
   

   TLegend* legend1 = new TLegend(0.65, 0.65, 0.8, 0.8);
   legend1->AddEntry(h_pf_pt, "PF", "l");
   legend1->AddEntry(h_puppi_pt, "PUPPI", "l");
   legend1->AddEntry(h_gnn_pt, "SSL", "l");
   legend1->Draw("same");
   
   canvas1->SaveAs("hist_pt.png");
   
   std::cout<<"Pt diff:"<<std::endl;
   std::cout<<"PF Mean:"<<median(pt_reso_pf)<<" PF stdv:"<<getResol(pt_reso_pf)<<std::endl;
   std::cout<<"PUPPI Mean:"<<median(pt_reso_puppi)<<" PUPPI stdv:"<<getResol(pt_reso_puppi)<<std::endl;
   std::cout<<"SSL Mean:"<<median(pt_reso_gnn)<<" SSL stdv:"<<getResol(pt_reso_gnn)<<std::endl;
   std::cout<<"PF Mean:"<<h_pf_pt_diff->GetMean()<<" PF RMS:"<<h_pf_pt_diff->GetRMS()<<std::endl;
   std::cout<<"PUPPI Mean:"<<h_puppi_pt_diff->GetMean()<<" PUPPI RMS:"<<h_puppi_pt_diff->GetRMS()<<std::endl;
   std::cout<<"PF_ak8 Mean:"<<h_pf_pt_diff_ak8->GetMean()<<" PF_ak8 RMS:"<<h_pf_pt_diff_ak8->GetRMS()<<std::endl;
   std::cout<<"PUPPI_ak8 Mean:"<<h_puppi_pt_diff_ak8->GetMean()<<" PUPPI_ak8 RMS:"<<h_puppi_pt_diff_ak8->GetRMS()<<std::endl;
   std::cout<<"SSL Mean:"<<h_gnn_pt_diff->GetMean()<<" SSL RMS:"<<h_gnn_pt_diff->GetRMS()<<std::endl;

   TCanvas* canvas2 = new TCanvas("canvas2", "Canvas2", 800, 600);
   h_pf_mass->SetDirectory(0);
   h_puppi_mass->SetDirectory(0);
   h_puppi_mass_ak8->SetDirectory(0);
   h_gnn_mass->SetDirectory(0);
   

   h_pf_mass->SetLineColor(kRed);
   h_puppi_mass->SetLineColor(kGreen);
   h_puppi_mass_ak8->SetLineColor(kOrange);
   h_gnn_mass->SetLineColor(kBlue);

   h_pf_mass->SetLineWidth(2);
   h_puppi_mass->SetLineWidth(2);
   h_puppi_mass_ak8->SetLineWidth(2);
   h_gnn_mass->SetLineWidth(2);

   h_gnn_mass->GetXaxis()->SetTitle("mass diff");
   h_gnn_mass->GetYaxis()->SetTitle("A.U.");

   h_pf_mass->Scale(1./h_pf_mass->Integral(),"width");
   h_puppi_mass->Scale(1./h_puppi_mass->Integral(),"width");
   h_puppi_mass_ak8->Scale(1./h_puppi_mass_ak8->Integral(),"width");
   h_gnn_mass->Scale(1./h_gnn_mass->Integral(),"width");

   h_gnn_mass->Draw("h");
   h_puppi_mass->Draw("h,same");
   h_puppi_mass_ak8->Draw("h,same");
   h_pf_mass->Draw("h,same");
   

   TLegend* legend2 = new TLegend(0.65, 0.65, 0.8, 0.8);
   legend2->AddEntry(h_pf_mass, "PF", "l");
   legend2->AddEntry(h_puppi_mass, "PUPPI", "l");
   legend2->AddEntry(h_puppi_mass_ak8, "PUPPI_AK8", "l");
   legend2->AddEntry(h_gnn_mass, "SSL", "l");
   legend2->Draw("same");
   
   canvas2->SaveAs("hist_mass.png");

  TCanvas* canvas3 = new TCanvas("canvas3", "Canvas3", 800, 600);
   h_LV_score->SetDirectory(0);
   h_PU_score->SetDirectory(0);
   h_NE_score->SetDirectory(0);

   h_LV_score->SetLineColor(kGreen);
   h_PU_score->SetLineColor(kBlue);
   h_NE_score->SetLineColor(kRed);

   h_LV_score->SetLineWidth(2);
   h_PU_score->SetLineWidth(2);
   h_NE_score->SetLineWidth(2);

   h_LV_score->GetXaxis()->SetTitle("SSL weight");
   h_PU_score->GetYaxis()->SetTitle("A.U.");

   h_LV_score->Scale(1./h_LV_score->Integral(),"width");
   h_PU_score->Scale(1./h_PU_score->Integral(),"width");
   h_NE_score->Scale(1./h_NE_score->Integral(),"width");

   h_NE_score->Draw("h");
   h_LV_score->Draw("h,same");
   h_PU_score->Draw("h,same");
   
   

   TLegend* legend3 = new TLegend(0.65, 0.65, 0.8, 0.8);
   legend3->AddEntry(h_LV_score, "LV chg", "l");
   legend3->AddEntry(h_PU_score, "PU chg", "l");
   legend3->AddEntry(h_NE_score, "neu ", "l");
   legend3->Draw("same");
   
   canvas3->SaveAs("hist_SSL_weight.png");


   int nBinsXgnn = h2gnn->GetNbinsX();
   double* xValsgnn = new double[nBinsXgnn];
   double* yValsgnn = new double[nBinsXgnn];
   int nBinsXpf = h2pf->GetNbinsX();
   double* xValspf = new double[nBinsXpf];
   double* yValspf = new double[nBinsXpf];
   int nBinsXpuppi = h2puppi->GetNbinsX();
   double* xValspuppi = new double[nBinsXpuppi];
   double* yValspuppi = new double[nBinsXpuppi];

   for (int ignn = 1; ignn <= nBinsXgnn; ++ignn) {
    TH1D* projYgnn = h2gnn->ProjectionY("projYgnn", ignn, ignn);
    xValsgnn[ignn - 1] = h2gnn->GetXaxis()->GetBinCenter(ignn);
    yValsgnn[ignn - 1] = projYgnn->GetRMS();
    delete projYgnn;
   }

   for (int ipf = 1; ipf <= nBinsXpf; ++ipf) {
    TH1D* projYpf = h2pf->ProjectionY("projYpf", ipf, ipf);
    xValspf[ipf - 1] = h2pf->GetXaxis()->GetBinCenter(ipf);
    yValspf[ipf - 1] = projYpf->GetRMS();
    delete projYpf;
   }

   for (int ipuppi = 1; ipuppi <= nBinsXpuppi; ++ipuppi) {
    TH1D* projYipuppi = h2puppi->ProjectionY("projYipuppi", ipuppi, ipuppi);
    xValspuppi[ipuppi - 1] = h2puppi->GetXaxis()->GetBinCenter(ipuppi);
    yValspuppi[ipuppi - 1] = projYipuppi->GetRMS();
    delete projYipuppi;
   }

    TGraph* gr1 = new TGraph(nBinsXpuppi, xValspuppi, yValspuppi);
    TGraph* gr2 = new TGraph(nBinsXgnn, xValsgnn, yValsgnn);
    TGraph* gr3 = new TGraph(nBinsXpf, xValspf, yValspf);

    
    gr1->SetMarkerColor(kRed);
    gr1->SetLineColor(kRed);
    gr2->SetMarkerColor(kGreen);
    gr2->SetLineColor(kGreen);
    gr3->SetMarkerColor(kBlue);
    gr3->SetLineColor(kBlue);

    TCanvas* c2d = new TCanvas("c2d", "Three TGraphs", 800, 600);
    c2d->SetLogy();
    gr1->Draw("AP");
    gr2->Draw("P SAME");
    gr3->Draw("P SAME");
    c2d->SaveAs("mass_vs_reso.png");*/
  
}
//
// member functions
//

// ------------ method called to produce the data  ------------
void SSLPuppiProducer::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput){
  client_->setBatchSize(1);
  auto const& pfs = iEvent.get(pf_token_);
  npf = 0;
  auto& input_0 = iInput.at("INPUT0");
   
  auto& input_1 = iInput.at("INPUT1");
  //input_0.setShape(1, 10);
  //input_1.setShape(0, 2);
  
  int num_node=0; int num_edge=0;
  std::vector<float> pf_eta, pf_phi, pf_pt;
  size_t i_pf = 0;
  int dim;

  for (const auto& pf_count : pfs){
    
    if (abs(pf_count.eta())>2.5) continue;
    num_node++;
  }
  input_0.setShape(0, num_node);
  npf = num_node;
  auto pfnode = input_0.allocate<float>();
  auto& vpfnode = (*pfnode)[0];
  vpfnode.clear();
  for (const auto& pf : pfs){
    dim = 0;
    if (abs(pf.eta())>2.5) continue;
    vpfnode.push_back(pf.eta());
    vpfnode.push_back(pf.phi());
    vpfnode.push_back(pf.pt());
    dim += 3;
    pf_eta.push_back(pf.eta()); pf_phi.push_back(pf.phi());pf_pt.push_back(pf.pt());
    if(pf.pt()<0.1) std::cout<<pf.pt()<<std::endl;
    if (pf.charge()!=0){
      vpfnode.push_back(1); vpfnode.push_back(0);vpfnode.push_back(0);dim += 3;
    }
    if (pf.charge()==0){
      if (pf.pdgId()==22){
        vpfnode.push_back(0);vpfnode.push_back(1);vpfnode.push_back(0);dim += 3;
      }
      else {
        vpfnode.push_back(0);vpfnode.push_back(0);vpfnode.push_back(1);dim += 3;
      }
    }
    if (pf.charge()!=0){
      if (pf.puppiWeight()==1){
        vpfnode.push_back(0); vpfnode.push_back(1);vpfnode.push_back(0);dim += 3;
      }
      else{
        vpfnode.push_back(1); vpfnode.push_back(0);vpfnode.push_back(0);dim += 3;
      }
    }
    if (pf.charge()==0){
      vpfnode.push_back(0); vpfnode.push_back(0);vpfnode.push_back(1);dim += 3;
    }
    vpfnode.push_back(0);dim += 1;
    i_pf++;
    if (i_pf == max_n_pf_)  break;
    //std::cout<<"node dim: "<<dim<<std::endl;
  }
  
  if (pf_eta.size()==0) return;
  for (unsigned m_count=0; m_count<pf_eta.size(); m_count++){
    for (unsigned n_count=0; n_count<pf_eta.size();n_count++){
      float dphi_, deta_, dR_;
      dphi_ = fabs(pf_phi[m_count]-pf_phi[n_count]);
      deta_ = fabs(pf_eta[m_count]-pf_eta[n_count]);
      if(dphi_>pi){
        float temp_;temp_ = std::ceil((dphi_ - pi)/(2*pi))*(2*pi);
        dphi_ = dphi_ - temp_;
      } 
      dR_ = sqrt(dphi_*dphi_ + deta_*deta_);
      if(m_count==n_count) continue;
      if(dR_<deltaRcut) num_edge++;
    }
  }
  input_1.setShape(1, num_edge);
  auto pfedge = input_1.allocate<long>();
  auto& vpfedge = (*pfedge)[0];
  vpfedge.clear();
  std::vector<long> edgeline1;std::vector<long> edgeline2;
  edgeline1.clear();edgeline2.clear();
  for (unsigned m=0; m<pf_eta.size(); m++){
    for (unsigned n=0; n<pf_eta.size();n++){
      float dphi, deta, dR;
      dphi = fabs(pf_phi[m]-pf_phi[n]);
      deta = fabs(pf_eta[m]-pf_eta[n]);
      if(dphi>pi){
        float temp;temp = std::ceil((dphi - pi)/(2*pi))*(2*pi);
        dphi = dphi - temp;
      } 
      dR = sqrt(dphi*dphi + deta*deta);
      
      if(dR<deltaRcut){
        long m_long, n_long;
        m_long = m; n_long = n;
        if(m==n) continue;
        //vpfedge.push_back(m_long);vpfedge.push_back(n_long);
        edgeline1.push_back(m_long);edgeline2.push_back(n_long);
      }
    }
  } 

  for(unsigned e1=0; e1<edgeline1.size(); e1++){
    vpfedge.push_back(edgeline1[e1]);
  }

  for(unsigned  e2=0; e2<edgeline2.size(); e2++){
    vpfedge.push_back(edgeline2[e2]);
  }
  vpfnode.resize(8*max_n_pf_);
  //vpfedge.resize(2*max_n_pf_);
  input_0.toServer(pfnode);
  input_1.toServer(pfedge);

}


void SSLPuppiProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {
   
  mass_truth_gnn_evt.clear(); pt_truth_gnn_evt.clear();mass_truth_puppi_evt.clear();mass_truth_chs_evt.clear();pt_truth_puppi_evt.clear();mass_truth_pf_evt.clear(); pt_truth_pf_evt.clear();pt_truth_chs_evt.clear();
  mass_reso_gnn_evt.clear(); mass_reso_pf_evt.clear(); mass_reso_puppi_evt.clear(); mass_reso_chs_evt.clear();pt_reso_gnn_evt.clear(); pt_reso_pf_evt.clear(); pt_reso_puppi_evt.clear();pt_reso_chs_evt.clear();
   
   const auto& output0 = iOutput.at("OUTPUT0");
   const auto& outputs = output0.fromServer<float>(); 
   //std::auto_ptr<std::vector<float> > SSLscore( new std::vector<float> );
   auto  SSLscore = std::make_unique<std::vector<float>>();
   auto  pf_eta = std::make_unique<std::vector<float>>();
   auto  pf_phi = std::make_unique<std::vector<float>>();
   auto  pf_puppipt = std::make_unique<std::vector<float>>();
   auto  gen_eta = std::make_unique<std::vector<float>>();
   auto  gen_phi = std::make_unique<std::vector<float>>();
   auto  gen_pt = std::make_unique<std::vector<float>>();
   auto  mass_diff = std::make_unique<std::vector<float>>();
   unsigned int i=0;
   auto const& pfs = iEvent.get(pf_token_);
   edm::Handle<std::vector<pat::PackedGenParticle>> genParticles;
   iEvent.getByToken(genParticleSrc_, genParticles);
   edm::Handle<std::vector<pat::Jet> > ak8puppijets;
   iEvent.getByToken(ak8JetsToken_,ak8puppijets);
   edm::Handle<std::vector<reco::GenJet> > ak8GENjets;
   iEvent.getByToken(GENJetsToken_,ak8GENjets);
   /*edm::Handle<std::vector<pat::Jet> > ak8SDjets;
   iEvent.getByToken(ak8SDJetsToken_,ak8SDjets);
   edm::Handle<std::vector<reco::GenJet> > ak8GENSDjets;
   iEvent.getByToken(GENSDJetsToken_,ak8GENSDjets);*/
   std::vector<fastjet::PseudoJet> pfJetInputs;
   std::vector<fastjet::PseudoJet> puppiJetInputs;
   std::vector<fastjet::PseudoJet> gnnJetInputs;
   std::vector<fastjet::PseudoJet> chsJetInputs;
   //bool keepIt=true;
   int LV_num,PU_num;
   LV_num = 0; PU_num = 0;
   for (const auto& pf_count : pfs){
    if ((abs(pf_count.eta())>2.5)){
     SSLscore->push_back(-1);
     continue;
    } 
    else SSLscore->push_back(outputs[0][i]);
    pf_eta->push_back(pf_count.eta());
    pf_phi->push_back(pf_count.phi());
    pf_puppipt->push_back(pf_count.pt()*outputs[0][i]);
    if ((pf_count.charge()!=0) && (pf_count.pt() > 0.5) && (pf_count.fromPV()>2)&& (pf_count.puppiWeight()>0.99)) LV_num++;
    if ((pf_count.charge()!=0) && (pf_count.pt() > 0.5) && (pf_count.fromPV()<1)&& (pf_count.puppiWeight()<0.01)) PU_num++;
    
    if (pf_count.pt()> 0.5) {
	    TLorentzVector pf_,puppi_,gnn_,chs_;
	    pf_.SetPtEtaPhiM(pf_count.pt(),pf_count.eta(),pf_count.phi(),0);
      //if (pf_count.charge()!=0) std::cout<<"pf_pt: "<<pf_count.pt()<<"puppi: "<<pf_count.puppiWeight()<<"SSLscore:"<<outputs[0][i]<<std::endl;
      if(pf_count.charge()!=0){
        if(pf_count.puppiWeight()>0.99) h_LV_score->Fill(outputs[0][i]);
        if(pf_count.puppiWeight()<0.01) h_PU_score->Fill(outputs[0][i]);
      }
      if(pf_count.charge()==0) h_NE_score->Fill(outputs[0][i]);
      if (pf_count.charge()==0) gnn_.SetPtEtaPhiM(pf_count.pt()*outputs[0][i],pf_count.eta(),pf_count.phi(),0);
      else gnn_.SetPtEtaPhiM(pf_count.pt()*pf_count.puppiWeight(),pf_count.eta(),pf_count.phi(),0);
      if (pf_count.charge()==0) chs_.SetPtEtaPhiM(pf_count.pt(),pf_count.eta(),pf_count.phi(),0);
      else chs_.SetPtEtaPhiM(pf_count.pt()*pf_count.puppiWeight(),pf_count.eta(),pf_count.phi(),0);
      puppi_.SetPtEtaPhiM(pf_count.pt()*pf_count.puppiWeight(),pf_count.eta(),pf_count.phi(),0);
   	  pfJetInputs.emplace_back(pf_.Px(), pf_.Py(), pf_.Pz(), pf_.E());
      gnnJetInputs.emplace_back(gnn_.Px(), gnn_.Py(), gnn_.Pz(), gnn_.E());
      puppiJetInputs.emplace_back(puppi_.Px(), puppi_.Py(), puppi_.Pz(), puppi_.E());
      chsJetInputs.emplace_back(chs_.Px(), chs_.Py(), chs_.Pz(), chs_.E());
        }
    i++;
   }
   
   std::vector<fastjet::PseudoJet> GenJetInputs;
   for(const auto& particle : *genParticles){
	   if(particle.status()!=1) continue;
     if(abs(particle.pdgId())==12) continue;
     if(abs(particle.pdgId())==14) continue;
     if(abs(particle.pdgId())==16) continue;
     if(abs(particle.eta())>2.5) continue;
    gen_eta->push_back(particle.eta());
    gen_pt->push_back(particle.pt());
    gen_phi->push_back(particle.phi());    
    if (particle.pt() > 0.5) {
            TLorentzVector pfgen_;
            pfgen_.SetPtEtaPhiM(particle.pt(),particle.eta(),particle.phi(),0);
            GenJetInputs.emplace_back(pfgen_.Px(), pfgen_.Py(), pfgen_.Pz(), pfgen_.E());
        }

   }
   fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetRadius_);
   fastjet::ClusterSequence csPf(pfJetInputs, jetDef);
   fastjet::ClusterSequence csGNN(gnnJetInputs, jetDef);
   fastjet::ClusterSequence csPUPPI(puppiJetInputs, jetDef);
   fastjet::ClusterSequence csCHS(chsJetInputs, jetDef);
   fastjet::ClusterSequence csGen(GenJetInputs, jetDef);
   float JetPtmin=150;
   std::vector<fastjet::PseudoJet> jetsPf = sorted_by_pt(csPf.inclusive_jets(JetPtmin));
   std::vector<fastjet::PseudoJet> jetsPUPPI = sorted_by_pt(csPUPPI.inclusive_jets(JetPtmin));
   std::vector<fastjet::PseudoJet> jetsGNN = sorted_by_pt(csGNN.inclusive_jets(JetPtmin));
   std::vector<fastjet::PseudoJet> jetsGen = sorted_by_pt(csGen.inclusive_jets(JetPtmin));
   std::vector<fastjet::PseudoJet> jetsCHS = sorted_by_pt(csCHS.inclusive_jets(JetPtmin));
   std::vector<float> massdiff;
   
   //if((LV_num<5)||(PU_num<50)) keepIt = false;
   //if(keepIt){
   //Matching & calculate invmass
   for (const auto& jetGen : jetsGen) {
    if((LV_num<5)||(PU_num<50)) continue;
      TLorentzVector genp4;
      genp4.SetPxPyPzE(jetGen.px(), jetGen.py(), jetGen.pz(), jetGen.e());
	    for (const auto& jetpf : jetsPf){
       TLorentzVector pfp4; 		   
		   pfp4.SetPxPyPzE(jetpf.px(), jetpf.py(), jetpf.pz(), jetpf.e());
       	  
		   if (pfp4.DeltaR(genp4)<0.1){
			   h_pf->Fill((pfp4.M()-genp4.M())/genp4.M());
         h_pf_pt_diff->Fill((pfp4.Pt()-genp4.Pt())/genp4.Pt());
         h_pf_mass->Fill(pfp4.M());
         h_pf_pt->Fill(pfp4.Pt());
         mass_reso_pf.push_back((pfp4.M()-genp4.M())/genp4.M());
         pt_reso_pf.push_back((pfp4.Pt()-genp4.Pt())/genp4.Pt());
         mass_reso_pf_evt.push_back((pfp4.M()-genp4.M())/genp4.M());
         pt_reso_pf_evt.push_back((pfp4.Pt()-genp4.Pt())/genp4.Pt());
         h2pf->Fill(genp4.M(), (pfp4.M()-genp4.M())/genp4.M());
         mass_truth_pf_evt.push_back(genp4.M()); pt_truth_pf_evt.push_back(genp4.Pt());
         
		   }
	   }
     for (const auto& jetpuppi : jetsPUPPI){
       TLorentzVector puppip4; 		   
		   puppip4.SetPxPyPzE(jetpuppi.px(), jetpuppi.py(), jetpuppi.pz(), jetpuppi.e());	
		   if (puppip4.DeltaR(genp4)<0.1){
			   h_puppi->Fill((puppip4.M()-genp4.M())/genp4.M());
         h_puppi_pt_diff->Fill((puppip4.Pt()-genp4.Pt())/genp4.Pt());
         h_puppi_mass->Fill(puppip4.M());h_puppi_pt->Fill(puppip4.Pt());
         mass_reso_puppi.push_back((puppip4.M()-genp4.M())/genp4.M());
         pt_reso_puppi.push_back((puppip4.Pt()-genp4.Pt())/genp4.Pt());
         mass_reso_puppi_evt.push_back((puppip4.M()-genp4.M())/genp4.M());
         pt_reso_puppi_evt.push_back((puppip4.Pt()-genp4.Pt())/genp4.Pt());
         h2puppi->Fill(genp4.M(), (puppip4.M()-genp4.M())/genp4.M());
         mass_truth_puppi_evt.push_back(genp4.M()); pt_truth_puppi_evt.push_back(genp4.Pt());
		   }
	   }

     for (const auto& jetCHS : jetsCHS){
      TLorentzVector chsp4; 		   
      chsp4.SetPxPyPzE(jetCHS.px(), jetCHS.py(), jetCHS.pz(), jetCHS.e());	
      if (chsp4.DeltaR(genp4)<0.1){
        h_chs->Fill((chsp4.M()-genp4.M())/genp4.M());
        h_chs_pt_diff->Fill((chsp4.Pt()-genp4.Pt())/genp4.Pt());
        h_chs_mass->Fill(chsp4.M());h_puppi_pt->Fill(chsp4.Pt());
        mass_reso_chs.push_back((chsp4.M()-genp4.M())/genp4.M());
        pt_reso_chs.push_back((chsp4.Pt()-genp4.Pt())/genp4.Pt());
        mass_reso_chs_evt.push_back((chsp4.M()-genp4.M())/genp4.M());
        pt_reso_chs_evt.push_back((chsp4.Pt()-genp4.Pt())/genp4.Pt());
        h2chs->Fill(genp4.M(), (chsp4.M()-genp4.M())/genp4.M());
        mass_truth_chs_evt.push_back(genp4.M()); pt_truth_chs_evt.push_back(genp4.Pt());
      }
    }
     for (const auto& jetGNN : jetsGNN){
       TLorentzVector GNNp4; 		   
		   GNNp4.SetPxPyPzE(jetGNN.px(), jetGNN.py(), jetGNN.pz(), jetGNN.e());		
		   if (GNNp4.DeltaR(genp4)<0.1){
			   massdiff.push_back((GNNp4.M()-genp4.M())/genp4.M());
			   h_gnn->Fill((GNNp4.M()-genp4.M())/genp4.M());
         h_gnn_pt_diff->Fill((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         h_gnn_mass->Fill(GNNp4.M());h_gnn_pt->Fill(GNNp4.Pt());
         mass_reso_gnn.push_back((GNNp4.M()-genp4.M())/genp4.M());
         pt_reso_gnn.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         mass_reso_gnn_evt.push_back((GNNp4.M()-genp4.M())/genp4.M());
         pt_reso_gnn_evt.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         h2gnn->Fill(genp4.M(), (GNNp4.M()-genp4.M())/genp4.M());
         mass_truth_gnn_evt.push_back(genp4.M()); pt_truth_gnn_evt.push_back(genp4.Pt());
		   }
	   }

     

     /*for (unsigned int ak8puppi = 0; ak8puppi< ak8puppijets->size(); ++ak8puppi){
       TLorentzVector ak8puppip4; 
       const pat::Jet & ak8puppip4jet = ak8puppijets->at(ak8puppi);		   
		   ak8puppip4.SetPxPyPzE(ak8puppip4jet.px(), ak8puppip4jet.py(), ak8puppip4jet.pz(), ak8puppip4jet.energy());	
       if (ak8puppip4.Pt()<JetPtmin) continue;
		   if (ak8puppip4.DeltaR(genp4)<0.1){
			   //massdiff.push_back((GNNp4.M()-genp4.M())/genp4.M());
			   h_puppi_ak8->Fill((ak8puppip4.M()-genp4.M())/genp4.M());
         h_puppi_pt_diff_ak8->Fill((ak8puppip4.Pt()-genp4.Pt())/genp4.Pt());
         //h_gnn_mass->Fill(GNNp4.M());h_gnn_pt->Fill(GNNp4.Pt());
         //mass_reso_gnn.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //mass_reso_gnn_evt.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn_evt.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //h2gnn->Fill(genp4.M(), (GNNp4.M()-genp4.M())/genp4.M());
         //mass_truth_gnn_evt.push_back(genp4.M()); pt_truth_gnn_evt.push_back(genp4.Pt());
		   }
	   }*/
    }

   for (unsigned int ak8puppi = 0; ak8puppi < ak8puppijets->size(); ++ak8puppi){
       TLorentzVector ak8puppip4; 
       const pat::Jet & ak8puppip4jet = ak8puppijets->at(ak8puppi);		   
		   ak8puppip4.SetPxPyPzE(ak8puppip4jet.px(), ak8puppip4jet.py(), ak8puppip4jet.pz(), ak8puppip4jet.energy());
       if (ak8puppip4.Pt()<JetPtmin) continue;
       for (unsigned int ak8GEN = 0; ak8GEN < ak8GENjets->size(); ++ak8GEN){
        	TLorentzVector ak8GENp4; 
       const reco::GenJet & ak8GENp4jet = ak8GENjets->at(ak8GEN);		   
		   ak8GENp4.SetPxPyPzE(ak8GENp4jet.px(), ak8GENp4jet.py(), ak8GENp4jet.pz(), ak8GENp4jet.energy());	
		   if (ak8puppip4.DeltaR(ak8GENp4)<0.1){
			   //massdiff.push_back((GNNp4.M()-genp4.M())/genp4.M());
			   h_puppi_ak8->Fill((ak8puppip4.M()-ak8GENp4.M())/ak8GENp4.M());
         h_puppi_pt_diff_ak8->Fill((ak8puppip4.Pt()-ak8GENp4.Pt())/ak8GENp4.Pt());
         h_puppi_mass_ak8->Fill(ak8puppip4.M());
         //h_gnn_mass->Fill(GNNp4.M());h_gnn_pt->Fill(GNNp4.Pt());
         //mass_reso_gnn.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //mass_reso_gnn_evt.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn_evt.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //h2gnn->Fill(genp4.M(), (GNNp4.M()-genp4.M())/genp4.M());
         //mass_truth_gnn_evt.push_back(genp4.M()); pt_truth_gnn_evt.push_back(genp4.Pt());
		   }
      }
	   }

    /*for (unsigned int ak8sd = 0; ak8sd < ak8SDjets->size(); ++ak8sd){
       TLorentzVector ak8sdp4; 
       const pat::Jet & ak8sdp4jet = ak8SDjets->at(ak8sd);		   
		   ak8sdp4.SetPxPyPzE(ak8sdp4jet.px(), ak8sdp4jet.py(), ak8sdp4jet.pz(), ak8sdp4jet.energy());
       if (ak8sdp4.Pt()<JetPtmin) continue;
       for (unsigned int ak8sdGEN = 0; ak8sdGEN < ak8GENSDjets->size(); ++ak8sdGEN){
        	TLorentzVector ak8sdGENp4; 
       const reco::GenJet & ak8sdGENp4jet = ak8GENSDjets->at(ak8sdGEN);		   
		   ak8sdGENp4.SetPxPyPzE(ak8sdGENp4jet.px(), ak8sdGENp4jet.py(), ak8sdGENp4jet.pz(), ak8sdGENp4jet.energy());	
		   if (ak8sdp4.DeltaR(ak8sdGENp4)<0.1){
			   //massdiff.push_back((GNNp4.M()-genp4.M())/genp4.M());
			   h_sd_ak8->Fill((ak8sdp4.M()-ak8sdGENp4.M())/ak8sdGENp4.M());
         h_sd_pt_diff_ak8->Fill((ak8sdp4.Pt()-ak8sdGENp4.Pt())/ak8sdGENp4.Pt());
         h_sd_mass_ak8->Fill(ak8sdp4.M());
         //h_gnn_mass->Fill(GNNp4.M());h_gnn_pt->Fill(GNNp4.Pt());
         //mass_reso_gnn.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //mass_reso_gnn_evt.push_back((GNNp4.M()-genp4.M())/genp4.M());
         //pt_reso_gnn_evt.push_back((GNNp4.Pt()-genp4.Pt())/genp4.Pt());
         //h2gnn->Fill(genp4.M(), (GNNp4.M()-genp4.M())/genp4.M());
         //mass_truth_gnn_evt.push_back(genp4.M()); pt_truth_gnn_evt.push_back(genp4.Pt());
		   }
      }
	   }*/
   
   /*iEvent.put(std::move(SSLscore),"SSLscore");
   iEvent.put(std::move(pf_eta),"pfeta");
   iEvent.put(std::move(pf_phi),"pfphi");
   iEvent.put(std::move(pf_puppipt),"pfpuppipt");
   iEvent.put(std::move(gen_eta),"geneta");
   iEvent.put(std::move(gen_phi),"genphi");
   iEvent.put(std::move(gen_pt),"genpt");*/
   iEvent.put(std::move(mass_diff),"massdiff");
   //tree_->Fill();
   
   //for(int i=0; i<npf; i++) {SSLscore->push_back(outputs[0][i]);}
   //for(int i=0; i<npf; i++) {std::cout<<outputs[0][i]<<std::endl;} 

   
  }



// ------------ method called when starting to processes a run  ------------
/*
void
TestProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TestProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TestProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TestProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SSLPuppiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("genParticleSrc", edm::InputTag("packedGenParticles"));
  desc.add<edm::InputTag>("pf_src", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("JetsAK8", edm::InputTag("slimmedJetsAK8"));
  desc.add<edm::InputTag>("GENJetsAK8", edm::InputTag("slimmedGenJetsAK8"));
  //desc.add<edm::InputTag>("JetsAK8SD", edm::InputTag("slimmedJetsAK8PFPuppiSoftDropPacked"));
  //desc.add<edm::InputTag>("GENSDJetsAK8", edm::InputTag("slimmedGenJetsAK8SoftDropSubJets"));
  desc.add<unsigned int>("max_n_pf", 4500);
  descriptions.add("SSLPuppiProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SSLPuppiProducer);
