#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

int debug = 0;
double Wjet_Normalization_FromBkgEstimation = 44.4878;

int VarBins = 0; 
double bins[9] = {600, 700, 800, 900, 1000, 1200, 1500, 2000, 2500};
int NBINS = 8;

// To get above normalization using below command with proper file name and path
//
// grep -A 10 "_WJets01_xww+++" /eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/cards_em_HP/other_wwlvj_Signal_aQGC600_em_HP.txt | grep "Events Number in sideband_low from fitting\|Events Number in sideband_high from fitting" | awk -F ":" '{print $2}' | awk '{print $1}'| awk '{total += $0} END{print "sum="total}'

//	FUNCTION TO CONVERT TF1 TO TH1F
TH1F* convertTF1toTH1F(TF1* f, TH1F* hIn){

  TH1F* hOut = (TH1F*)hIn->Clone();
  hOut->Reset();

  for (unsigned int ibin=1; ibin<hOut->GetNbinsX()+1; ++ibin){
    double xlow = hOut->GetBinLowEdge(ibin);
    double xhigh = hOut->GetBinLowEdge(ibin+1);
	hOut->SetBinContent(ibin,f->Integral(xlow,xhigh)/(xhigh-xlow));
  }

  return hOut;
}

//	FUNCTION TO GET WJET INTO SIGNAL REGION FROM SIDEBAND USING ALPHA
TH1F* getWjetSignalRegion_usingAlpha(TH1F* wjet, TH1F* alpha){
   
   TH1F* hOut = (TH1F*)wjet->Clone();
   hOut->Reset();

   int AvgAlphaInInterestedRegion_CountFreq = 0;
   double AvgAlphaInInterestedRegion = 0.0;
   // Take the appropriate choice of average value of alpha
   for (unsigned int ibin=1; ibin<alpha->GetNbinsX()+1; ++ibin){
      if ( (alpha->GetBinLowEdge(ibin) > 2500) && (alpha->GetBinContent(ibin) > 0))
      {
      	AvgAlphaInInterestedRegion += alpha->GetBinContent(ibin);
	//cout<<"===> "<<alpha->GetBinContent(ibin)<<endl;
	AvgAlphaInInterestedRegion_CountFreq++;
      }
   }
   cout<<"\n\n AvgAlphaInInterestedRegion_CountFreq = "<<AvgAlphaInInterestedRegion_CountFreq <<"\n\n"<<endl;
   if (AvgAlphaInInterestedRegion_CountFreq == 0)
   {
     for (unsigned int ibin=1; ibin<alpha->GetNbinsX()+1; ++ibin){
      if ( (alpha->GetBinLowEdge(ibin) > 1500) && (alpha->GetBinContent(ibin) > 0))
      {
      	AvgAlphaInInterestedRegion += alpha->GetBinContent(ibin);
	//cout<<"iii===> "<<alpha->GetBinContent(ibin)<<endl;
	AvgAlphaInInterestedRegion_CountFreq++;
      }
     }
   }
   cout<<"\n\n AvgAlphaInInterestedRegion_CountFreq = "<<AvgAlphaInInterestedRegion_CountFreq <<"\n\n"<<endl;

   AvgAlphaInInterestedRegion = AvgAlphaInInterestedRegion/((double)AvgAlphaInInterestedRegion_CountFreq);

   for (unsigned int ibin=1; ibin<hOut->GetNbinsX()+1; ++ibin){
      if (alpha->GetBinContent(ibin)>0) hOut->SetBinContent(ibin, wjet->GetBinContent(ibin)*alpha->GetBinContent(ibin));
      else hOut->SetBinContent(ibin, wjet->GetBinContent(ibin)*AvgAlphaInInterestedRegion);
      hOut->SetBinError(ibin, 0.0);
   }

   return hOut;
}

//	RESET A HISTOGRAM TO THE DIFFERENT NUMBER OF BINS...
TH1F *ResetTo4bins(TH1F *wjet, int nbins, double xmin, double xmax) {
  if (debug) cout << "For " << wjet->GetName() << "\n\n" << endl;

  TH1F *h;
  if (VarBins)
  	h = new TH1F("hOut", ";M_{ZV} (Vjet Signal Region);Events", NBINS, bins);
  else
  	h = new TH1F("hOut", ";M_{ZV} (Vjet Signal Region);Events", nbins, xmin, xmax);
  if (wjet) for (Int_t i = 1; i <= wjet->GetNbinsX(); i++){
              h->Fill(wjet->GetBinCenter(i), wjet->GetBinContent(i));
	      if (debug) cout << i << "\t" << wjet->GetBinContent(i) << "\t" << wjet->GetBinLowEdge(i) << " , " << wjet->GetBinLowEdge(i+1) << endl;
	}
  if (debug){
  cout << h->GetBinContent(1) << endl;
  cout << h->GetBinContent(2) << endl;
  cout << h->GetBinContent(3) << endl;
  cout << h->GetBinContent(4) << endl;
  cout << h->GetBinContent(5) << endl;
  }
  // Add overflow bin
  h->SetBinContent(h->GetNbinsX() , h->GetBinContent(h->GetNbinsX()) + h->GetBinContent(h->GetNbinsX()+1) );
  return h;
}


void GetAll_Systematic_Shape_ZV() {


   string functionForFit = "pol1";
   char* fitFormula;	fitFormula = new char[functionForFit.size()+1];	strcpy(fitFormula, functionForFit.c_str());

   int nbins = 4;
   double xmin_fit = 600.0;
   double xmax_fit = 3000.0;
   
  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 3;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  //example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
  //  example_plot( iPeriod, 11 );  // left-aligned
  //  example_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  example_plot( iPeriod, 11 );  // default: left-aligned
  //  example_plot( iPeriod, 22 );  // centered
  //  example_plot( iPeriod, 33 );  // right-aligned  


   
   cout<< "\n\n===============\n\n \t TO GET CORRECT NORMALIZATION USE (with proper path of file): \n\ngrep -A 10 \"_WJets01_xww+++\" /eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/cards_em_HP/other_wwlvj_Signal_aQGC600_em_HP.txt | grep \"Events Number in sideband_low from fitting\\|Events Number in sideband_high from fitting\" | awk -F \":\" '{print $2}' | awk '{print $1}'| awk '{total += $0} END{print \"sum=\"total}' \n\n===============\n\n" << endl;
 
   // Open the file containing the tree.
   TFile *myFile = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/SecondStep/WWTree_After_CWR/2019_03_28_16h05/HaddedFiles/Hadds_for_BkgEstimation/WWTree_VJets.root","READ");

   // Open all necessary file that we get after background estimation:
   TFile *bkgEstFile = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto.root","READ");
   TFile *bkgEstFile_Up0 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Up_0.root","READ");
   TFile *bkgEstFile_Up1 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Up_2.root","READ");
   TFile *bkgEstFile_Down0 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Down_0.root","READ");
   TFile *bkgEstFile_Down1 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Down_2.root","READ");
   TFile *bkgEstFile_alternate = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_After_CWR_ZV/2019_03_31_09h49/wjetmodel_Ex__WJets01_xww__sb_lo_Exp_auto.root","READ");

   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in. (otree is the name of tree)
   TTreeReader fReader("otree", myFile);


   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<Float_t> LHEWeight = {fReader, "LHEWeight"};
   TTreeReaderValue<Int_t> event = {fReader, "event"};
   TTreeReaderValue<Int_t> nTotEvents = {fReader, "nTotEvents"};
   TTreeReaderValue<Int_t> nTotNegEvents = {fReader, "nTotNegEvents"};
   TTreeReaderValue<Int_t> nEvents = {fReader, "nEvents"};
   TTreeReaderValue<Int_t> nNegEvents = {fReader, "nNegEvents"};
   TTreeReaderValue<Int_t> nPV = {fReader, "nPV"};
   TTreeReaderValue<Float_t> wSampleWeight = {fReader, "wSampleWeight"};
   TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
   TTreeReaderValue<Float_t> trig_eff_Weight = {fReader, "trig_eff_Weight"};
   TTreeReaderValue<Float_t> trig_eff_Weight2 = {fReader, "trig_eff_Weight2"};
   TTreeReaderValue<Float_t> id_eff_Weight = {fReader, "id_eff_Weight"};
   TTreeReaderValue<Float_t> id_eff_Weight2 = {fReader, "id_eff_Weight2"};
   TTreeReaderValue<Float_t> totalEventWeight = {fReader, "totalEventWeight"};
   TTreeReaderValue<Float_t> pu_Weight = {fReader, "pu_Weight"};
   TTreeReaderValue<Float_t> totalEventWeight_2 = {fReader, "totalEventWeight_2"};
   TTreeReaderValue<Float_t> pu_Weight_up = {fReader, "pu_Weight_up"};
   TTreeReaderValue<Float_t> totalEventWeight_3 = {fReader, "totalEventWeight_3"};
   TTreeReaderValue<Float_t> totalEventWeight_2Lep = {fReader, "totalEventWeight_2Lep"};
   TTreeReaderValue<Float_t> pu_Weight_down = {fReader, "pu_Weight_down"};
   TTreeReaderValue<Float_t> pfMET = {fReader, "pfMET"};
   TTreeReaderValue<Float_t> pfMET_Phi = {fReader, "pfMET_Phi"};
   TTreeReaderValue<Float_t> pfMET_Corr = {fReader, "pfMET_Corr"};
   TTreeReaderValue<Float_t> pfMET_Corr_phi = {fReader, "pfMET_Corr_phi"};
   TTreeReaderValue<Float_t> nu_pz_type0 = {fReader, "nu_pz_type0"};
   TTreeReaderValue<Int_t> nu_pz_isre = {fReader, "nu_pz_isre"};
   TTreeReaderValue<Int_t> type = {fReader, "type"};
   TTreeReaderValue<Float_t> l_pt1 = {fReader, "l_pt1"};
   TTreeReaderValue<Float_t> l_eta1 = {fReader, "l_eta1"};
   TTreeReaderValue<Float_t> l_phi1 = {fReader, "l_phi1"};
   TTreeReaderValue<Float_t> l_e1 = {fReader, "l_e1"};
   TTreeReaderValue<Float_t> l_iso1 = {fReader, "l_iso1"};
   TTreeReaderValue<Float_t> l_pt2 = {fReader, "l_pt2"};
   TTreeReaderValue<Float_t> l_eta2 = {fReader, "l_eta2"};
   TTreeReaderValue<Float_t> l_phi2 = {fReader, "l_phi2"};
   TTreeReaderValue<Float_t> l_e2 = {fReader, "l_e2"};
   TTreeReaderValue<Float_t> dilep_pt = {fReader, "dilep_pt"};
   TTreeReaderValue<Float_t> dilep_eta = {fReader, "dilep_eta"};
   TTreeReaderValue<Float_t> dilep_phi = {fReader, "dilep_phi"};
   TTreeReaderValue<Float_t> dilep_m = {fReader, "dilep_m"};
   TTreeReaderValue<Float_t> l_iso2 = {fReader, "l_iso2"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt = {fReader, "ungroomed_AK8jet_pt"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jer = {fReader, "ungroomed_AK8jet_pt_jer"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_eta = {fReader, "ungroomed_AK8jet_eta"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_phi = {fReader, "ungroomed_AK8jet_phi"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_e = {fReader, "ungroomed_AK8jet_e"};
   TTreeReaderValue<Float_t> AK8jet_mass = {fReader, "AK8jet_mass"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr = {fReader, "AK8jet_mass_pr"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jer = {fReader, "AK8jet_mass_pr_jer"};
   TTreeReaderValue<Float_t> AK8jet_mass_so = {fReader, "AK8jet_mass_so"};
   TTreeReaderValue<Float_t> AK8jet_pt_so = {fReader, "AK8jet_pt_so"};
   TTreeReaderValue<Float_t> AK8jet_mass_tr = {fReader, "AK8jet_mass_tr"};
   TTreeReaderValue<Float_t> AK8jet_mass_fi = {fReader, "AK8jet_mass_fi"};
   TTreeReaderValue<Float_t> AK8jet_tau2tau1 = {fReader, "AK8jet_tau2tau1"};
   TTreeReaderValue<Float_t> AK8_jetID_loose = {fReader, "AK8_jetID_loose"};
   TTreeReaderValue<Float_t> AK8jet_qjet = {fReader, "AK8jet_qjet"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt = {fReader, "ungroomed_PuppiAK8_jet_pt"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jer = {fReader, "ungroomed_PuppiAK8_jet_pt_jer"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_eta = {fReader, "ungroomed_PuppiAK8_jet_eta"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_phi = {fReader, "ungroomed_PuppiAK8_jet_phi"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_e = {fReader, "ungroomed_PuppiAK8_jet_e"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass = {fReader, "PuppiAK8_jet_mass"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr = {fReader, "PuppiAK8_jet_mass_pr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jer = {fReader, "PuppiAK8_jet_mass_pr_jer"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_so = {fReader, "PuppiAK8_jet_mass_so"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_so_corr = {fReader, "PuppiAK8_jet_mass_so_corr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_pt_so = {fReader, "PuppiAK8_jet_pt_so"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_tr = {fReader, "PuppiAK8_jet_mass_tr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_fi = {fReader, "PuppiAK8_jet_mass_fi"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_tau2tau1 = {fReader, "PuppiAK8_jet_tau2tau1"};
   TTreeReaderValue<Float_t> PuppiAK8_jetID_loose = {fReader, "PuppiAK8_jetID_loose"};
   TTreeReaderValue<Float_t> PuppiAK8jet_qjet = {fReader, "PuppiAK8jet_qjet"};
   TTreeReaderValue<Float_t> isGen = {fReader, "isGen"};
   TTreeReaderValue<Float_t> v_pt_type0 = {fReader, "v_pt_type0"};
   TTreeReaderValue<Float_t> v_eta_type0 = {fReader, "v_eta_type0"};
   TTreeReaderValue<Float_t> v_phi = {fReader, "v_phi"};
   TTreeReaderValue<Float_t> v_mt_type0 = {fReader, "v_mt_type0"};
   TTreeReaderValue<Float_t> v_mass_type0 = {fReader, "v_mass_type0"};
   TTreeReaderValue<Float_t> mass_lvj_type0 = {fReader, "mass_lvj_type0"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jer = {fReader, "mass_lvj_type0_met_jer"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8 = {fReader, "mass_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> mt_lvj_type0_PuppiAK8 = {fReader, "mt_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> pt_lvj_type0_PuppiAK8 = {fReader, "pt_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> eta_lvj_type0_PuppiAK8 = {fReader, "eta_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> rapidity_lvj_type0_PuppiAK8 = {fReader, "rapidity_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> phi_lvj_type0_PuppiAK8 = {fReader, "phi_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> energy_lvj_type0_PuppiAK8 = {fReader, "energy_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> mass_llj_PuppiAK8 = {fReader, "mass_llj_PuppiAK8"};
   TTreeReaderValue<Float_t> pt_llj_PuppiAK8 = {fReader, "pt_llj_PuppiAK8"};
   TTreeReaderValue<Float_t> eta_llj_PuppiAK8 = {fReader, "eta_llj_PuppiAK8"};
   TTreeReaderValue<Float_t> phi_llj_PuppiAK8 = {fReader, "phi_llj_PuppiAK8"};
   TTreeReaderValue<Int_t> njets = {fReader, "njets"};
   TTreeReaderValue<Int_t> njetsPuppi = {fReader, "njetsPuppi"};
   TTreeReaderValue<Int_t> nGoodAK8jets = {fReader, "nGoodAK8jets"};
   TTreeReaderValue<Int_t> nGoodPuppiAK8jets = {fReader, "nGoodPuppiAK8jets"};
   TTreeReaderValue<Int_t> njets_unmerged = {fReader, "njets_unmerged"};
   TTreeReaderValue<Int_t> njetsPuppi_unmerged = {fReader, "njetsPuppi_unmerged"};
   TTreeReaderValue<Int_t> nBTagJet_loose = {fReader, "nBTagJet_loose"};
   TTreeReaderValue<Int_t> nBTagJet_medium = {fReader, "nBTagJet_medium"};
   TTreeReaderValue<Int_t> nBTagJet_tight = {fReader, "nBTagJet_tight"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_loose = {fReader, "nBTagJetPuppi_loose"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_medium = {fReader, "nBTagJetPuppi_medium"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_tight = {fReader, "nBTagJetPuppi_tight"};
   TTreeReaderValue<Int_t> nBTagJet_loose_unmerged = {fReader, "nBTagJet_loose_unmerged"};
   TTreeReaderValue<Int_t> nBTagJet_medium_unmerged = {fReader, "nBTagJet_medium_unmerged"};
   TTreeReaderValue<Int_t> nBTagJet_tight_unmerged = {fReader, "nBTagJet_tight_unmerged"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_loose_unmerged = {fReader, "nBTagJetPuppi_loose_unmerged"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_medium_unmerged = {fReader, "nBTagJetPuppi_medium_unmerged"};
   TTreeReaderValue<Int_t> nBTagJetPuppi_tight_unmerged = {fReader, "nBTagJetPuppi_tight_unmerged"};
   TTreeReaderValue<Float_t> L1_Prefweight = {fReader, "L1_Prefweight"};
   TTreeReaderValue<Float_t> btag0Wgt = {fReader, "btag0Wgt"};
   TTreeReaderValue<Float_t> btag1Wgt = {fReader, "btag1Wgt"};
   TTreeReaderValue<Float_t> btag2Wgt = {fReader, "btag2Wgt"};
   TTreeReaderValue<Float_t> btag0WgtUpHF = {fReader, "btag0WgtUpHF"};
   TTreeReaderValue<Float_t> btag0WgtDownHF = {fReader, "btag0WgtDownHF"};
   TTreeReaderValue<Float_t> btag0WgtUpLF = {fReader, "btag0WgtUpLF"};
   TTreeReaderValue<Float_t> btag0WgtDownLF = {fReader, "btag0WgtDownLF"};
   TTreeReaderValue<Float_t> btag1WgtUpHF = {fReader, "btag1WgtUpHF"};
   TTreeReaderValue<Float_t> btag1WgtDownHF = {fReader, "btag1WgtDownHF"};
   TTreeReaderValue<Float_t> btag1WgtUpLF = {fReader, "btag1WgtUpLF"};
   TTreeReaderValue<Float_t> btag1WgtDownLF = {fReader, "btag1WgtDownLF"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt = {fReader, "vbf_maxpt_j1_pt"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jer = {fReader, "vbf_maxpt_j1_pt_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta = {fReader, "vbf_maxpt_j1_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jer = {fReader, "vbf_maxpt_j1_eta_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_phi = {fReader, "vbf_maxpt_j1_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_e = {fReader, "vbf_maxpt_j1_e"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_mass = {fReader, "vbf_maxpt_j1_mass"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_bDiscriminatorCSV = {fReader, "vbf_maxpt_j1_bDiscriminatorCSV"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt = {fReader, "vbf_maxpt_j2_pt"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jer = {fReader, "vbf_maxpt_j2_pt_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta = {fReader, "vbf_maxpt_j2_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jer = {fReader, "vbf_maxpt_j2_eta_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_phi = {fReader, "vbf_maxpt_j2_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_e = {fReader, "vbf_maxpt_j2_e"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_mass = {fReader, "vbf_maxpt_j2_mass"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_bDiscriminatorCSV = {fReader, "vbf_maxpt_j2_bDiscriminatorCSV"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_pt = {fReader, "vbf_maxpt_jj_pt"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_eta = {fReader, "vbf_maxpt_jj_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_phi = {fReader, "vbf_maxpt_jj_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_m = {fReader, "vbf_maxpt_jj_m"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_Deta = {fReader, "vbf_maxpt_jj_Deta"};
   TTreeReaderValue<Float_t> LepWEta = {fReader, "LepWEta"};
   TTreeReaderValue<Float_t> LepWRapidity = {fReader, "LepWRapidity"};
   TTreeReaderValue<Float_t> HadWEta = {fReader, "HadWEta"};
   TTreeReaderValue<Float_t> HadWRapidity = {fReader, "HadWRapidity"};
   TTreeReaderValue<Float_t> WWEta = {fReader, "WWEta"};
   TTreeReaderValue<Float_t> WWEta_PuppiAK8 = {fReader, "WWEta_PuppiAK8"};
   TTreeReaderValue<Float_t> WWRapidity = {fReader, "WWRapidity"};
   TTreeReaderValue<Float_t> WWRapidity_PuppiAK8 = {fReader, "WWRapidity_PuppiAK8"};
   TTreeReaderValue<Float_t> ZeppenfeldWH = {fReader, "ZeppenfeldWH"};
   TTreeReaderValue<Float_t> ZeppenfeldWHPuppi = {fReader, "ZeppenfeldWHPuppi"};
   TTreeReaderValue<Float_t> costheta1Puppi_type0 = {fReader, "costheta1Puppi_type0"};
   TTreeReaderValue<Float_t> costheta2Puppi_type0 = {fReader, "costheta2Puppi_type0"};
   TTreeReaderValue<Float_t> costhetastarPuppi_type0 = {fReader, "costhetastarPuppi_type0"};
   TTreeReaderValue<Float_t> phiPuppi_type0 = {fReader, "phiPuppi_type0"};
   TTreeReaderValue<Float_t> phi1Puppi_type0 = {fReader, "phi1Puppi_type0"};
   TTreeReaderValue<Float_t> VBSCentralityPuppi_type0 = {fReader, "VBSCentralityPuppi_type0"};
   TTreeReaderValue<Float_t> RpTPuppi_type0 = {fReader, "RpTPuppi_type0"};
   TTreeReaderValue<Float_t> ZeppenfeldWLPuppi_type0 = {fReader, "ZeppenfeldWLPuppi_type0"};
   TTreeReaderValue<Float_t> LeptonProjectionPuppi_type0 = {fReader, "LeptonProjectionPuppi_type0"};
   TTreeReaderValue<Float_t> PtBalancePuppi_type0 = {fReader, "PtBalancePuppi_type0"};
   TTreeReaderValue<Float_t> costheta1_type0 = {fReader, "costheta1_type0"};
   TTreeReaderValue<Float_t> costheta2_type0 = {fReader, "costheta2_type0"};
   TTreeReaderValue<Float_t> costhetastar_type0 = {fReader, "costhetastar_type0"};
   TTreeReaderValue<Float_t> phi_type0 = {fReader, "phi_type0"};
   TTreeReaderValue<Float_t> phi1_type0 = {fReader, "phi1_type0"};
   TTreeReaderValue<Float_t> VBSCentrality_type0 = {fReader, "VBSCentrality_type0"};
   TTreeReaderValue<Float_t> RpT_type0 = {fReader, "RpT_type0"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0 = {fReader, "ZeppenfeldWL_type0"};
   TTreeReaderValue<Float_t> LeptonProjection_type0 = {fReader, "LeptonProjection_type0"};
   TTreeReaderValue<Float_t> PtBalance_type0 = {fReader, "PtBalance_type0"};
   TTreeReaderValue<Float_t> BosonCentrality_type0 = {fReader, "BosonCentrality_type0"};
   TTreeReaderValue<Float_t> BosonCentralityPuppi_type0 = {fReader, "BosonCentralityPuppi_type0"};
   TTreeReaderValue<Float_t> PtBalance_2Lep = {fReader, "PtBalance_2Lep"};
   TTreeReaderValue<Float_t> BosonCentrality_2Lep = {fReader, "BosonCentrality_2Lep"};
   TTreeReaderValue<Float_t> costheta1_2Lep = {fReader, "costheta1_2Lep"};
   TTreeReaderValue<Float_t> costheta2_2Lep = {fReader, "costheta2_2Lep"};
   TTreeReaderValue<Float_t> costhetastar_2Lep = {fReader, "costhetastar_2Lep"};
   TTreeReaderValue<Float_t> phi_2Lep = {fReader, "phi_2Lep"};
   TTreeReaderValue<Float_t> phi1_2Lep = {fReader, "phi1_2Lep"};
   TTreeReaderValue<Float_t> VBSCentrality_2Lep = {fReader, "VBSCentrality_2Lep"};
   TTreeReaderValue<Float_t> RpT_2Lep = {fReader, "RpT_2Lep"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_2Lep = {fReader, "ZeppenfeldWL_2Lep"};
   TTreeReaderValue<Float_t> LeptonProjection_2Lep = {fReader, "LeptonProjection_2Lep"};


   TFile *f;
   if (VarBins)	f = new TFile("ZV_bkg_estimation_8bins_ChargedHiggs.root", "RECREATE");
   else f = new TFile("ZV_bkg_estimation_4Bins_LeadingLep50GeV.root", "RECREATE");

   // Create a histogram for the values we read.
   TH1F *hMC_Signal_4bin;
   if (VarBins)
   	hMC_Signal_4bin   = new TH1F("hMC_Signal_4bin",   "hMC_Signal_4bin;M_{ZV};Events", NBINS, bins);
   else
   	hMC_Signal_4bin   = new TH1F("hMC_Signal_4bin",   "hMC_Signal_4bin;M_{ZV};Events", 4, 600, 2500);
   TH1F *hMC_Signal_15bin  = new TH1F("hMC_Signal_15bin",  "hMC_Signal_15bin;M_{ZV};Events", 15, 600, 3000);
   TH1F *hMC_Signal_88bin  = new TH1F("hMC_Signal_88bin",  "hMC_Signal_88bin;M_{ZV};Events", 88, 600, 3000);

   //TH1F *hSideBand_4bin  = new TH1F("hSideBand_4bin",  "hSideBand_4bin;M_{ZV};Events", 4, 600, 3000);
   TH1F *hSideBand_15bin = new TH1F("hSideBand_15bin", "hSideBand_15bin;M_{ZV};Events", 15, 600, 3000);
   TH1F *hSideBand_88bin = new TH1F("hSideBand_88bin", "hSideBand_88bin;M_{ZV};Events", 88, 600, 3000);


   // Loop over all entries of the TTree.
   while (fReader.Next()) {
      // Just access the data as if variables were iterators (note the '*' in front of them):
      if(!(*l_pt1>50 && (((*type==0)&&(abs(*l_eta1)<2.4)) || ((*type==1)&&((abs(*l_eta1)<2.5)&&!(abs(*l_eta1)>1.4442 && abs(*l_eta1)<1.566)))))) continue;
      if(!(*l_pt2>30 && (((*type==0)&&(abs(*l_eta2)<2.4)) || ((*type==1)&&((abs(*l_eta2)<2.5)&&!(abs(*l_eta2)>1.4442 && abs(*l_eta2)<1.566)))))) continue;
      if(!((*ungroomed_PuppiAK8_jet_pt>200)&&(abs(*ungroomed_PuppiAK8_jet_eta)<2.4)&&(*PuppiAK8_jet_tau2tau1<0.55))) continue;
      //if(!((*PuppiAK8_jet_mass_so_corr>65) && (*PuppiAK8_jet_mass_so_corr<105))) continue;
      if(!(*nBTagJet_loose==0)) continue;
      if(!(*dilep_m>76 && *dilep_m<107)) continue;
      if(!(*vbf_maxpt_jj_m>800)) continue;
      if(!(abs(*vbf_maxpt_j2_eta - *vbf_maxpt_j1_eta)>4.0)) continue;
      if(!((*vbf_maxpt_j1_pt>30) && (*vbf_maxpt_j2_pt>30))) continue;
      if(!(*mass_llj_PuppiAK8>600)) continue;
      

      // Fill histogram for signal region
      if ((*PuppiAK8_jet_mass_so_corr>65) && (*PuppiAK8_jet_mass_so_corr<105))
      {
      	hMC_Signal_4bin->Fill(*mass_llj_PuppiAK8, ((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      	hMC_Signal_15bin->Fill(*mass_llj_PuppiAK8,((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      	hMC_Signal_88bin->Fill(*mass_llj_PuppiAK8,((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      }
      else if ((*PuppiAK8_jet_mass_so_corr>40) && (*PuppiAK8_jet_mass_so_corr<150))
      {
       //hSideBand_4bin->Fill(*mass_llj_PuppiAK8,((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      	hSideBand_15bin->Fill(*mass_llj_PuppiAK8,((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      	hSideBand_88bin->Fill(*mass_llj_PuppiAK8,((*wSampleWeight)*(35867.06)*(*L1_Prefweight)*(*btag0Wgt)*(*pu_Weight)*(*totalEventWeight_2Lep)));
      }
      //cout<<*mass_llj_PuppiAK8<<"\t"<<*BosonCentrality_type0<<endl;
   }

   // include overflow bin
    hMC_Signal_4bin->SetBinContent(hMC_Signal_4bin->GetNbinsX(), hMC_Signal_4bin->GetBinContent(hMC_Signal_4bin->GetNbinsX())+ hMC_Signal_4bin->GetBinContent(hMC_Signal_4bin->GetNbinsX()+1));
   hMC_Signal_15bin->SetBinContent(hMC_Signal_15bin->GetNbinsX(),hMC_Signal_15bin->GetBinContent(hMC_Signal_15bin->GetNbinsX())+hMC_Signal_15bin->GetBinContent(hMC_Signal_15bin->GetNbinsX()+1));
   hMC_Signal_88bin->SetBinContent(hMC_Signal_88bin->GetNbinsX(),hMC_Signal_88bin->GetBinContent(hMC_Signal_88bin->GetNbinsX())+hMC_Signal_88bin->GetBinContent(hMC_Signal_88bin->GetNbinsX()+1));
   //hSideBand_4bin->SetBinContent(hSideBand_4bin->GetNbinsX(),  hSideBand_4bin->GetBinContent(hSideBand_4bin->GetNbinsX())+  hSideBand_4bin->GetBinContent(hSideBand_4bin->GetNbinsX()+1));
    hSideBand_15bin->SetBinContent(hSideBand_15bin->GetNbinsX(), hSideBand_15bin->GetBinContent(hSideBand_15bin->GetNbinsX())+ hSideBand_15bin->GetBinContent(hSideBand_15bin->GetNbinsX()+1));
    hSideBand_88bin->SetBinContent(hSideBand_88bin->GetNbinsX(), hSideBand_88bin->GetBinContent(hSideBand_88bin->GetNbinsX())+ hSideBand_88bin->GetBinContent(hSideBand_88bin->GetNbinsX()+1));

    if (debug){
    	cout << "============== Cross-Check	===============\n\n" << endl;
    	for (int i=1; i<=hMC_Signal_15bin->GetNbinsX()+1; i++)
    		cout<< i << "\t" << hMC_Signal_15bin->GetBinContent(i) << endl;

    	cout << "============== Cross-Check	===============\n\n" << endl;
    	for (int i=1; i<=hSideBand_15bin->GetNbinsX()+1; i++)
    		cout<< i << "\t" << hSideBand_15bin->GetBinContent(i) << endl;
    }

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART- 1: Get alpha by taking ratio of Vjet signal region and sideband region distribution. 
   //
   //	Note that in alpha over-flow bin in already included as the both signal and side-band distribution for wjets includes overflow bin.
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////
   TH1F* alpha = (TH1F*)hMC_Signal_15bin->Clone();
   alpha->Divide(hSideBand_15bin);
   alpha->SetMinimum(-1.0);
   alpha->SetMaximum(4.0);
   alpha->SetName("alpha");
   alpha->GetXaxis()->SetTitle("M_{ZV} (GeV)");
   alpha->GetYaxis()->SetTitle("alpha = #frac{N^{MC,V+jets}_{signal}}{N^{MC,V+jets}_{side-band}}");

   // Fit alpha using polynomial of order 1
   gStyle->SetOptFit(1);
   TF1* f1 = new TF1("f1","pol1",600,3000);
   f1->SetLineColor(2);
   TFitResultPtr fitRes = alpha->Fit("f1","S");

   TCanvas* c1 = new TCanvas("Alpha_SigmaBand","canvas", 1000,750);

   // Get sigma band for alpha
   TH1F *hConf1s = (TH1F*)alpha->Clone("hConf1s");
   TH1F *hConf2s = (TH1F*)alpha->Clone("hConf2s");
   hConf1s->Reset();
   hConf2s->Reset();
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hConf1s, 0.68);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hConf2s, 0.95);
   hConf1s->SetName("hAlpha_1Sigma_Band");
   hConf2s->SetName("hAlpha_2Sigma_Band");
   hConf1s->SetStats(kFALSE);
   hConf2s->SetStats(kFALSE);
   hConf1s->SetMarkerSize(0);
   hConf2s->SetMarkerSize(0);
   hConf1s->SetFillColor(kGreen+1);
   hConf2s->SetFillColor(kOrange);
   hConf1s->SetLineColor(kGreen+1);
   hConf2s->SetLineColor(kOrange);
   hConf2s->SetTitle("");
   hConf2s->Draw("e3");
   hConf1s->Draw("e3 same");
   alpha->Draw("E1 same");
   f1->Draw("same");

   TLegend* leg = new TLegend(0.20,0.95,0.65,0.75);
   leg->SetNColumns(2);
   leg->AddEntry(alpha, "MC", "lep");
   leg->AddEntry(f1, "Fit", "l");
   leg->AddEntry(hConf1s, "1 #sigma band","f");
   leg->AddEntry(hConf2s, "2 #sigma band","f");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.65, 0.75, 0.95,0.95, "brNDC");
   char chi2[200], a[200],b[200];
   sprintf(chi2, "Chi2/ndf = %f / %d",f1->GetChisquare(),f1->GetNDF());
   sprintf(a, "a = %f +/- %f",f1->GetParameter(0),f1->GetParError(0));
   sprintf(b, "b = %f +/- %f",f1->GetParameter(1),f1->GetParError(1));
   pt->AddText(chi2);   pt->AddText(a);   pt->AddText(b);
   pt->Draw();   c1->Draw();   c1->Write();

   c1->Clear();
   leg->Clear();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART- 2: Get Alpha function variation by +/- 1 sigma
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
   c1->SetName("AlphaSyst_ParVariation");

   alpha->SetStats(0);
   alpha->SetTitle("");
   alpha->Draw("E1");
   f1->SetLineColor(1);
   f1->Draw("E1 same");
   leg->AddEntry(alpha, "MC", "lep");
   leg->AddEntry(f1, "Fit", "l");
   // Get Correlation matrix...
   //float initialPars[2],initialParErrors[2];
   std::vector<double> initialPars;
   std::vector<double> initialParErrors;

   for(int ipar=0; ipar<2; ipar++)
   {
	initialPars.push_back(f1->GetParameter(ipar));
	initialParErrors.push_back(f1->GetParError(ipar));
	cout << "par " << ipar << "  : " << f1->GetParameter(ipar) << " +/- " << f1->GetParError(ipar) << endl; 
   }

   TMatrixDSym cov = fitRes->GetCovarianceMatrix();
   cov.Print();
   int numFitParameter = cov.GetNrows();
   //int numFitParameter = f1->GetNpar();
   //Find the Eigenvectors and Eigen value of the covariance matrix
   TMatrixDSymEigen eigencov = cov;
   //eigencov.Print();
   //Get the matrix of Eigen vector
   const TMatrixD eigenvector_matrix = eigencov.GetEigenVectors();
   const TVectorD eigenvalues = eigencov.GetEigenValues();

   eigenvector_matrix.Print();
   eigenvalues.Print();
   
   std::vector<TVectorD*> eigenvectors;
   cout << "eigenvector_matrix.GetNrows() = " <<  eigenvector_matrix.GetNrows() << endl;

   for ( int i = 0; i < eigenvector_matrix.GetNrows(); ++i ) {
     TVectorD* eigenvector = new TVectorD(eigenvector_matrix.GetNrows());
     for ( int j = 0; j < eigenvector_matrix.GetNrows(); ++j ) {
       (*eigenvector)(j) =  eigenvector_matrix(j, i);
     }
     eigenvectors.push_back(eigenvector);
     std::cout << "EigenValue #" << i << " = " << eigenvalues(i) << ": EigenVector = { ";
     for ( int j = 0; j < eigenvector_matrix.GetNrows(); ++j ) {
       std::cout << (*eigenvector)(j);
       if ( j < (eigenvector_matrix.GetNrows() - 1) ) std::cout << ",";
     }
     std::cout << " }" << std::endl;
   }

  
   std::vector<TF1*> systFunctions;
   std::vector<TString> names;
   for( int i = 0; i < numFitParameter; ++i ) {
     //Get the i-th Eigenvector
     const TVectorD& eigenvector = (*eigenvectors[i]);
     //Compute norm of i-th Eigenvector
     double norm = eigenvector.Norm2Sqr();
     //Get the Eigenvalue associate to i-th Eigenvector
     double eigenvalue = eigenvalues[i];
     //Compute uncertainty on fit parameter
     double sigma = sqrt(TMath::Abs(eigenvalue));
     if (debug)
       cout << i << " norm " << norm << " " << eigenvalue << " +/- " << sigma << endl;
     //Compute unit vector in direction of i-th Eigenvector
     TVectorD eigenvector_unit = (1./norm)*eigenvector;
     if (debug)
       cout << i << " eigenvector " << eigenvector[i] << "; unit eigenvector " << eigenvector_unit[i] <<endl;
     
     //    Vary parameters up 
     std::vector<double> upPars;
     for (int j=0;j<initialPars.size();++j){
       double newParUp = initialPars[j] + sigma*eigenvector_unit[j];
       upPars.push_back(newParUp);
     }
     TString f12UpName = Form("Par%iUp", i);
     TF1* f12Up = new TF1(f12UpName.Data(),fitFormula,xmin_fit,xmax_fit);
     for (int j=0;j<initialPars.size();++j){
       f12Up->FixParameter(j,upPars[j]);
     }
     if (i==0) {
       f12Up->SetLineColor(kBlue);
     }else{
       f12Up->SetLineColor(kRed);
     }
     f12Up->Draw("same");
     leg->AddEntry(f12Up,f12UpName);
     systFunctions.push_back(f12Up);
     names.push_back(f12UpName);
     
     //    Vary parameters down
     std::vector<double> downPars;
     for (int j=0;j<initialPars.size();++j){
       double newParDown = initialPars[j] - sigma*eigenvector_unit[j];
       downPars.push_back(newParDown);
     }
     TString f12DownName = Form("Par%iDown", i);
     TF1* f12Down = new TF1(f12DownName.Data(),fitFormula,xmin_fit,xmax_fit);
     for (int j=0;j<initialPars.size();++j){
       f12Down->SetParameter(j,downPars[j]);
     }
     if (i==0){
       f12Down->SetLineColor(kBlue);
       f12Down->SetLineStyle(2);
     } else {
       f12Down->SetLineColor(kRed);
       f12Down->SetLineStyle(2);
     }
     f12Down->Draw("same");
     leg->AddEntry(f12Down,f12DownName);
     systFunctions.push_back(f12Down);
     names.push_back(f12DownName);
     
     if (debug) {
       for (int j=0;j<initialPars.size();++j){
         cout << " _________ " << j << endl;
         cout << "Initial " << initialPars[j] << endl;
         cout << "Up   " << upPars[j] << endl;
         cout << "Down " << downPars[j] << endl;      
       }
     }
   }

   leg->Draw();	pt->Draw();

   TLine *line =new TLine(c1->GetUxmin(),0.0,c1->GetUxmax(),0.0);
   line->SetLineColor(kViolet);
   line->Draw();

   c1->Draw();
   c1->Write();
   for (int k=0; k<systFunctions.size(); k++)
	systFunctions[k]->Write();

   if (debug) {
     cout << "Number of fit functions " << systFunctions.size() << endl;
     cout << "Order of functions is such: " << endl;
     cout << "[0] -- par1 Up" << endl;
     cout << "[1] -- par1 Down" << endl;
     cout << "[2] -- par2 Up" << endl;
     cout << "[3] -- par2 Down" << endl;
     cout << "Etc......" << endl;
   }

   //First convert nominal fit into histogram

   c1->Clear();
   leg->Clear();
   c1->SetName("Alpha_Hist_Systematic_fits");
   TH1F* hCentral = convertTF1toTH1F(f1, hMC_Signal_88bin);
   hCentral->SetName("AlphaSyst_Nominal");
   hCentral->SetTitle("");
   hCentral->SetLineColor(2);
   hCentral->SetMarkerSize(1.5);
   hCentral->SetMarkerStyle(5);

   alpha->Draw();
   hCentral->Draw("same");
   leg->AddEntry(alpha,"Alpha MC","lep");
   leg->AddEntry(hCentral,"Alpha Nominal","l");
   cout<<"Alphaa MC = "<< alpha->Integral() <<"\t nominal = "<< hCentral->Integral() << endl;

   //Now convert systematic fits into histos
   std::vector<TH1F*> histos;
   TString catname = "AlphaSyst_";
   for (int ifun=0; ifun<systFunctions.size(); ++ifun){
     TF1* thisFunc = systFunctions[ifun];
     TH1F* thisHisto = convertTF1toTH1F(thisFunc, hMC_Signal_88bin);
     TString thisHistoName = catname+names[ifun];
     thisHisto->SetName(thisHistoName);
     thisHisto->SetTitle(thisHistoName);
     thisHisto->SetLineColor(ifun+3);
     thisHisto->SetMarkerColor(ifun+3);
     thisHisto->Draw("same");
     histos.push_back(thisHisto);
     leg->AddEntry(thisHisto,thisHistoName);
     cout<<"\t DEGUB: 4 ==> "<<names[ifun]<<"\t"<< thisFunc->GetName()<< ", \t Integral = "<< thisHisto->Integral() << endl;
   }
   leg->Draw();	pt->Draw();
   line->Draw();
   c1->Draw();
   c1->Write();
   leg->Clear();
   c1->Clear();

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 3: Get side-band Nominal Main histogram which is corrected using data 
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////
   //File is opened at start of program so that TFile will not interfer with the output TFile

   TLegend* leg2 = new TLegend(0.45,0.95,0.95,0.75);
   leg2->SetNColumns(2);
   c1->SetName("Compare_Vjet_SideBand_MC_CorrShape");
   TH1F* Wjet_Corr_Hist = (TH1F*)bkgEstFile->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist->SetStats(0);
   Wjet_Corr_Hist->SetLineColor(2);
   Wjet_Corr_Hist->SetMarkerColor(2);
   Wjet_Corr_Hist->SetName("Vjet_SideBand_CorrShape_From_Data");
   Wjet_Corr_Hist->SetTitle("Comparison of Vjet MC and corrected shape from data");
   Wjet_Corr_Hist->Write();	// Need to write explicitly. Somehow its not written automatically to output root file.

   Wjet_Corr_Hist->Draw("hist");
   hSideBand_88bin->SetStats(0);
   hSideBand_88bin->Draw("E1 same");
   leg2->AddEntry(Wjet_Corr_Hist,"SideBand Region (Corr)");
   leg2->AddEntry(hSideBand_88bin,"SideBand Region (MC)","lep");

   leg2->Draw();
   c1->Draw();
   c1->Write();
   c1->SetName("Compare_Vjet_SideBand_MC_CorrShape_Log");
   c1->SetLogy(1);
   c1->Write();
   c1->SetLogy(0);
   c1->Clear();
   leg2->Clear();

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 4: Get signal histogram from alpha multiplication
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////

   c1->SetName("AlphaSyst_Vjet_CorrShape_SignalRegion");
   TH1F* hCorr_Signal_Central = getWjetSignalRegion_usingAlpha(Wjet_Corr_Hist,hCentral);
   hCorr_Signal_Central->SetName("AlphaSyst_Vjet_SR_Nominal");
   hCorr_Signal_Central->SetTitle("Vjets Signal Region Corrected From Data");
   hCorr_Signal_Central->SetStats(0);
   hCorr_Signal_Central->SetLineColor(2);
   hCorr_Signal_Central->SetMarkerColor(2);
   hCorr_Signal_Central->SetMarkerSize(1.5);
   hCorr_Signal_Central->SetMarkerStyle(5);

   hMC_Signal_88bin->SetStats(0);
   hMC_Signal_88bin->SetTitle("Vjets Signal Region Corrected From Data");
   hMC_Signal_88bin->Draw();
   hCorr_Signal_Central->Draw("same");

   leg2->AddEntry(hMC_Signal_88bin,"Vjet MC","lep");
   leg2->AddEntry(hCorr_Signal_Central,"Vjet Corr Nominal","l");

   std::vector<TH1F*> histos_Vjet_SR;
   histos_Vjet_SR.push_back(hCorr_Signal_Central);
   catname = "AlphaSyst_Vjet_SR_";
   for (int ifun=0; ifun<histos.size(); ++ifun){
     TH1F* thisFunc = histos[ifun];
     TH1F* thisHisto = getWjetSignalRegion_usingAlpha(Wjet_Corr_Hist, thisFunc);
     TString thisHistoName = catname+names[ifun];
     thisHisto->SetName(thisHistoName);
     thisHisto->SetTitle(thisHistoName);
     //thisHisto->Write();
     thisHisto->SetLineColor(ifun+3);
     thisHisto->SetMarkerColor(ifun+3);
     thisHisto->SetStats(0);
     histos_Vjet_SR.push_back(thisHisto);
     thisHisto->Draw("same");
     TString thisLegName = "Vjet Corr "+names[ifun];
     leg2->AddEntry(thisHisto,thisLegName);
     cout<<"\t DEGUB: 3 ==> "<<names[ifun]<<"\t"<< thisFunc->GetName()<< ", Integral = "<<thisHisto->Integral()<<endl;
   }
   leg2->Draw();
   c1->Draw();
   c1->Write();
   c1->SetName("AlphaSyst_Vjet_CorrShape_SignalRegion_Log");
   c1->SetLogy(1);
   c1->Write();
   leg2->Clear();
   c1->Clear();
   c1->SetLogy(0);


   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 5: Reset to 4 bins
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////
   
   c1->SetName("AlphaSyst_Vjet_CorrShape_SignalRegion_4Bins");
   catname = "AlphaSyst_Vjet_SR_4bins_";
   cout << "\n====================== \n\n\t histos for alpha \n\n ======================"<< endl;
   for (int ifun=0; ifun<histos_Vjet_SR.size(); ++ifun){
     TH1F* thisFunc = histos_Vjet_SR[ifun];
     TH1F* thisHisto = ResetTo4bins(thisFunc, 4, 600, 2500);
     TString thisHistoName;
     if (ifun == 0) thisHistoName = catname+"Nominal";
     else thisHistoName = catname+names[ifun-1];
     thisHisto->SetName(thisHistoName);
     //thisHisto->SetTitle(thisHistoName);
     thisHisto->SetTitle("");
     thisHisto->SetStats(0);
     thisHisto->SetLineColor(ifun+2);
     thisHisto->SetMarkerColor(ifun+2);
     if (ifun == 0) thisHisto->Draw();
     else thisHisto->Draw("same");
     TString thisLegName;
     if (ifun == 0) thisLegName = "Vjet Corr Nominal";
     else thisLegName = "Vjet Corr "+names[ifun-1];
     leg2->AddEntry(thisHisto,thisLegName);
     if (ifun == 0) cout<<"\t DEGUB: 2 ==> Nominal"<<endl;
     else cout<<"\t DEGUB: 2 ==> "<<names[ifun-1]<<"\t"<< thisFunc->GetName() << ", Integral = "<<thisHisto->Integral()<<endl;
   }
   hMC_Signal_4bin->Draw("same");
   leg2->AddEntry(hMC_Signal_4bin,"MC","lep");
   leg2->Draw();	
   c1->Draw();
   c1->Write();
   c1->SetName("AlphaSyst_Vjet_CorrShape_SignalRegion_4Bins_Log");
   c1->SetLogy(1);
   c1->Write();
   leg2->Clear();
   c1->Clear();
   c1->SetLogy(0);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 6: Get Wjet fit all ystematics shapes including alternate shape 
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////
   
   std::vector<TH1F*> histos_WjetSyst;
   //histos_WjetSyst.push_back(Wjet_Corr_Hist)
   c1->SetName("WjetFitSyst_Vjet_CorrShape_SideBandRegion");
   TH1F* Wjet_Corr_Hist_Up0 = (TH1F*)bkgEstFile_Up0->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist_Up0->SetName("WjetFitSyst_SideBandRegion_Corr_Hist_From_Data_Par0Up");
   Wjet_Corr_Hist_Up0->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist_Up0->SetStats(0);
   Wjet_Corr_Hist_Up0->SetLineColor(3);
   Wjet_Corr_Hist_Up0->SetMarkerColor(3);
   histos_WjetSyst.push_back(Wjet_Corr_Hist_Up0);
   Wjet_Corr_Hist_Up0->Write();

   TH1F* Wjet_Corr_Hist_Down0 = (TH1F*)bkgEstFile_Down0->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist_Down0->SetName("WjetFitSyst_SideBandRegion_Corr_Hist_From_Data_Par0Down");
   Wjet_Corr_Hist_Down0->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist_Down0->SetStats(0);
   Wjet_Corr_Hist_Down0->SetLineColor(5);
   Wjet_Corr_Hist_Down0->SetMarkerColor(5);
   histos_WjetSyst.push_back(Wjet_Corr_Hist_Down0);
   Wjet_Corr_Hist_Down0->Write();

   TH1F* Wjet_Corr_Hist_Up1 = (TH1F*)bkgEstFile_Up1->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist_Up1->SetName("WjetFitSyst_SideBandRegion_Corr_Hist_From_Data_Par1Up");
   Wjet_Corr_Hist_Up1->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist_Up1->SetStats(0);
   Wjet_Corr_Hist_Up1->SetLineColor(4);
   Wjet_Corr_Hist_Up1->SetMarkerColor(4);
   histos_WjetSyst.push_back(Wjet_Corr_Hist_Up1);
   Wjet_Corr_Hist_Up1->Write();

   TH1F* Wjet_Corr_Hist_Down1 = (TH1F*)bkgEstFile_Down1->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist_Down1->SetName("WjetFitSyst_SideBandRegion_Corr_Hist_From_Data_Par1Down");
   Wjet_Corr_Hist_Down1->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist_Down1->SetStats(0);
   Wjet_Corr_Hist_Down1->SetLineColor(6);
   Wjet_Corr_Hist_Down1->SetMarkerColor(6);
   histos_WjetSyst.push_back(Wjet_Corr_Hist_Down1);
   Wjet_Corr_Hist_Down1->Write();

   TH1F* Wjet_Corr_Hist_alter = (TH1F*)bkgEstFile_alternate->Get("rrv_mass_lvj__rrv_mass_lvj");
   Wjet_Corr_Hist_alter->SetName("WjetFitSyst_SideBandRegion_Corr_Hist_From_Data_AlterNameShape");
   Wjet_Corr_Hist_alter->Scale(Wjet_Normalization_FromBkgEstimation);
   Wjet_Corr_Hist_alter->SetStats(0);
   Wjet_Corr_Hist_alter->SetLineColor(7);
   Wjet_Corr_Hist_alter->SetMarkerColor(7);
   histos_WjetSyst.push_back(Wjet_Corr_Hist_alter);
   Wjet_Corr_Hist_alter->Write();
   names.push_back("AlternateShape_Down");

   hSideBand_88bin->SetTitle("Vjet SideBand Region Corrected From Data");
   hSideBand_88bin->Draw("E1");
   Wjet_Corr_Hist->Draw("hist same");
   Wjet_Corr_Hist_Up0->Draw("hist same");
   Wjet_Corr_Hist_Up1->Draw("hist same");
   Wjet_Corr_Hist_Down0->Draw("hist same");
   Wjet_Corr_Hist_Down1->Draw("hist same");
   Wjet_Corr_Hist_alter->Draw("hist same");


   leg2->AddEntry(hSideBand_88bin,"MC","lep");
   leg2->AddEntry(Wjet_Corr_Hist,"Vjet Corr Nominal");
   leg2->AddEntry(Wjet_Corr_Hist_Up0,"Vjet Corr Par0 Up");
   leg2->AddEntry(Wjet_Corr_Hist_Up1,"Vjet Corr Par1 Up");
   leg2->AddEntry(Wjet_Corr_Hist_Down0,"Vjet Corr Par0 Down");
   leg2->AddEntry(Wjet_Corr_Hist_Down1,"Vjet Corr Par1 Down");
   leg2->AddEntry(Wjet_Corr_Hist_alter,"Vjet Corr Alternate");
   leg2->Draw();

   c1->Draw();
   c1->Write();
   c1->SetName("WjetFitSyst_Vjet_CorrShape_SideBandRegion_Log");
   c1->SetLogy(1);
   c1->Write();
   leg2->Clear();
   c1->Clear();
   c1->SetLogy(0);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 7: Get Wjet fit systematics (in signal region by multiplying it with alpha)
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   c1->SetName("WjetFitSyst_Vjet_CorrShape_SignalRegion");

   hMC_Signal_88bin->SetTitle("Vjets Signal Region Corrected From Data");
   leg2->AddEntry(hMC_Signal_88bin,"MC","lep");
   hMC_Signal_88bin->Draw("E1");

   hCorr_Signal_Central->SetLineColor(2);
   hCorr_Signal_Central->SetMarkerColor(2);
   hCorr_Signal_Central->Draw("hist same");
   leg2->AddEntry(hCorr_Signal_Central,"Vjet Nominal");

   std::vector<TH1F*> histos_WjetSyst_SR;
   histos_WjetSyst_SR.push_back(hCorr_Signal_Central);
   catname = "WjetFitSyst_SignalRegion_Corr_Hist_From_Data_";
   cout<<"\n\n\n=======================\t RK\n\n\n"<<endl;
   for (int ifun=0; ifun<histos_WjetSyst.size(); ++ifun){
     TH1F* thisFunc = histos_WjetSyst[ifun];
     TH1F* thisHisto = getWjetSignalRegion_usingAlpha(thisFunc, hCentral);
     TString thisHistoName;
     thisHistoName = catname+names[ifun];
     thisHisto->SetName(thisHistoName);
     thisHisto->SetStats(0);
     thisHisto->SetLineColor(ifun+3);
     thisHisto->SetMarkerColor(ifun+3);
     histos_WjetSyst_SR.push_back(thisHisto);
     thisHisto->Draw("same");
     TString thisLegName = "Vjet "+names[ifun];
     leg2->AddEntry(thisHisto,thisLegName);
     cout<<"\t==> "<<names[ifun]<<"\t"<< thisFunc->GetName()<<endl;
   }
   cout<<"\n\n\n=======================\t RK\n\n\n"<<endl;
   leg2->Draw();

   c1->Draw();
   c1->Write();
   c1->SetName("WjetFitSyst_Vjet_CorrShape_SignalRegion_Log");
   c1->SetLogy(1);
   c1->Write();
   leg2->Clear();
   c1->Clear();
   c1->SetLogy(0);

   
   ////////////////////////////////////////////////////////////////////////////////////////////////
   //
   //	PART - 8: Get Wjet fit systematics (in signal region by multiplying it with alpha)
   //
   //	Change to 4 bins
   //
   ////////////////////////////////////////////////////////////////////////////////////////////////

   c1->SetName("WjetFitSyst_Vjet_CorrShape_SignalRegion_4bins");
   catname = "WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_";
   std::vector<TH1F*> histos_WjetSyst_SR_4bin;
   for (int ifun=0; ifun<histos_WjetSyst_SR.size(); ++ifun){
     TH1F* thisFunc = histos_WjetSyst_SR[ifun];
     TH1F* thisHisto = ResetTo4bins(thisFunc, 4, 600, 2500);
     TString thisHistoName;
     if (ifun == 0) thisHistoName = catname+"Nominal";
     else thisHistoName = catname+names[ifun-1];
     thisHisto->SetName(thisHistoName);
     thisHisto->SetTitle("Vjets Signal Region Corrected From Data");
     //thisHisto->SetTitle(thisHistoName);
     thisHisto->SetTitle("");
     thisHisto->SetStats(0);
     thisHisto->SetLineColor(ifun+2);
     thisHisto->SetMarkerColor(ifun+2);
     //thisHisto->SetBinContent(thisHisto->GetNbinsX(), thisHisto->GetBinContent(thisHisto->GetNbinsX()) + thisHisto->GetBinContent(thisHisto->GetNbinsX()+1) );
     histos_WjetSyst_SR_4bin.push_back(thisHisto);
     if (ifun == 0) thisHisto->Draw();
     else thisHisto->Draw("same");
     TString thisLegName;
     if (ifun == 0) thisLegName = "Vjet Nominal";
     else thisLegName = "Vjet "+names[ifun-1];
     leg2->AddEntry(thisHisto,thisLegName);
     if (ifun == 0) cout << "\t DEGUB: 1 ==> Nominal"<<endl;
     else cout<<"\t DEGUB: 1 ==> "<<names[ifun-1]<<"\t"<< thisFunc->GetName()<<endl;
   }
   //	Get Altername shape down :
   //	we have only one alternate shape so just get the difference of main and alternate shape bin by bin and make another alternate shape...
   TH1F* mainWjetHist = histos_WjetSyst_SR_4bin[0];
   TH1F* alterWjetHist = histos_WjetSyst_SR_4bin[histos_WjetSyst_SR.size()-1];
   cout<<" ==> "<< mainWjetHist->GetName() << "\t" << alterWjetHist->GetName() << endl;
   TH1F* hAlter_WjetHist_Up;
   if (VarBins)
   	hAlter_WjetHist_Up = new TH1F("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_Up",";M_{ZV};Events", NBINS, bins);
   else
   	hAlter_WjetHist_Up = new TH1F("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_Up",";M_{ZV};Events",4, 600, 2500);
   for (int ibin=1; ibin<mainWjetHist->GetNbinsX()+2; ibin++)
   {
      hAlter_WjetHist_Up->SetBinContent(ibin,mainWjetHist->GetBinContent(ibin) + (mainWjetHist->GetBinContent(ibin) - alterWjetHist->GetBinContent(ibin)));
      cout<< ibin << "\t" << mainWjetHist->GetBinContent(ibin) << "\t" << alterWjetHist->GetBinContent(ibin) << "\t" << hAlter_WjetHist_Up->GetBinContent(ibin) << endl;
   }
   //hAlter_WjetHist_Up->SetBinContent(hAlter_WjetHist_Up->GetNbinsX(), hAlter_WjetHist_Up->GetBinContent(hAlter_WjetHist_Up->GetNbinsX()) + hAlter_WjetHist_Up->GetBinContent(hAlter_WjetHist_Up->GetNbinsX()+1) );
   hAlter_WjetHist_Up->SetLineColor(histos_WjetSyst_SR.size()+2);
   hAlter_WjetHist_Up->SetMarkerColor(histos_WjetSyst_SR.size()+2);

   hAlter_WjetHist_Up->Draw("same");
   leg2->AddEntry(hAlter_WjetHist_Up,"Vjet AlternateShape Up");
   hMC_Signal_4bin->Draw("same");
   leg2->AddEntry(hMC_Signal_4bin,"MC","lep");
   leg2->Draw();	
   c1->Draw();
   c1->Write();
   c1->SetName("WjetFitSyst_Vjet_CorrShape_SignalRegion_4bins_Log");
   c1->SetLogy(1);
   c1->Write();
   leg2->Clear();
   c1->Clear();
   c1->SetLogy(0);


    f1->Write();   f->Write();   f->Close();
}
