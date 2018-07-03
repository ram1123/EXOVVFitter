#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

int debug = 0;
double Wjet_Normalization_FromBkgEstimation = 140.626122846+101.117465308;
// To get above normalization using below command with proper file name and path
//
// grep -A 10 "_WJets01_xww+++" WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/cards_em_HP/other_wwlvj_Signal_aQGC600_em_HP.txt | grep "Events Number in sideband_low from fitting\|Events Number in sideband_high from fitting" | awk -F ":" '{print $2}' | awk '{print $1}'| awk '{total += $0} END{print "sum="total}'

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

TH1F* getWjetSignalRegion_usingAlpha(TH1F* wjet, TH1F* alpha){
   
   TH1F* hOut = (TH1F*)wjet->Clone();
   hOut->Reset();

   for (unsigned int ibin=1; ibin<hOut->GetNbinsX()+1; ++ibin){
      hOut->SetBinContent(ibin, wjet->GetBinContent(ibin)*alpha->GetBinContent(ibin));
      hOut->SetBinError(ibin, 0.0);
   }

   return hOut;
}

TH1F *ResetTo4bins(TH1F *wjet, int nbins, double xmin, double xmax) {
  if (debug) cout << "For " << wjet->GetName() << "\n\n" << endl;
  TH1F *h = new TH1F("hOut", ";M_{WW} (Vjet Signal Region);Events", nbins, xmin, xmax);
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
  return h;
}

//TH1F* ResetTo4bins(TH1F* wjet, int nbins, double xmin, double xmax){
//
//   TH1F* hOut = new TH1F("hOut",";M_{WW} (Vjet Signal Region);",nbins, xmin, xmax);
//
//   double binWidth = (xmax - xmin)/nbins;
//
//   int binFill = 0;
//
//   double hOutBinContent[4] = {0.0, 0.0, 0.0, 0.0};
//   cout << "bins = " << (sizeof(hOutBinContent) / sizeof(double)) -1 << endl;
//
//   cout << " total bins = "<< wjet->GetNbinsX()+1 << endl;
//   for (unsigned int ibin=1; ibin<wjet->GetNbinsX()+1; ++ibin){
//     double xlow = wjet->GetBinLowEdge(ibin);
//     double xhigh = wjet->GetBinLowEdge(ibin+1);
//
//     double xlow_New = hOut->GetBinLowEdge(binFill);
//     double xhigh_New = hOut->GetBinLowEdge(binFill+1);
//
//     cout<< "binFill = " << binFill << endl;
//     if (xlow < xlow_New){	cout<<" Don't fill lower bins..." << endl;cout<< "Debug 1: "<< binFill << endl;	}
//     else{
//        if (xhigh <= xhigh_New)
//	{
//     	   hOutBinContent[binFill] += wjet->GetBinContent(ibin);
//	   cout<< "Debug 2: "<< binFill << endl;
//	}
//        else
//	{
//     	   //binFill += 1;
//	   if (binFill < ((sizeof(hOutBinContent) / sizeof(double)) - 1 )) binFill++;
//	   else { std::cout << "Ouch!" << std::endl; break; }
//	   hOutBinContent[binFill] += wjet->GetBinContent(ibin);
//	   cout<< "Debug 3: "<< binFill << endl;
//	}
//     }
//     
//     cout<< "binFill = " << binFill << endl;
//     cout<<ibin <<  "  xlow = " <<  xlow  <<  " xhigh = " << xhigh <<  "  xlow_New = " <<  xlow_New <<  "  xhigh_New = " <<  xhigh_New <<  "  hOutBinContent[" << binFill << "] = " << hOutBinContent[binFill] <<endl;
//
//
//   }
//
//   return hOut;
//
//}


void Alpha_calculation() {

   cout<< "\n\n===============\n\n \t TO GET CORRECT NORMALIZATION USE (with proper path of file): \n\ngrep -A 10 \"_WJets01_xww+++\" WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/cards_em_HP/other_wwlvj_Signal_aQGC600_em_HP.txt | grep \"Events Number in sideband_low from fitting\\|Events Number in sideband_high from fitting\" | awk -F \":\" '{print $2}' | awk '{print $1}'| awk '{total += $0} END{print \"sum=\"total}' \n\n===============\n\n" << endl;

   string functionForFit = "pol1";

   char* fitFormula;	fitFormula = new char[functionForFit.size()+1];	strcpy(fitFormula, functionForFit.c_str());
   double xmin_fit = 600.0;
   double xmax_fit = 5000.0;
   
 
   // Open the file containing the tree.
   TFile *myFile = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/SecondStep/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15/HaddedFiles/Hadds_for_BkgEstimation/WWTree_VJets.root","READ");

   // Open file that we get after background estimation:
   TFile *bkgEstFile = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/2018_06_22_07h37/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto.root","READ");
   TFile *bkgEstFile_Up0 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/2018_06_22_07h37/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Up_0.root","READ");
   TFile *bkgEstFile_Up1 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/2018_06_22_07h37/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Up_2.root","READ");
   TFile *bkgEstFile_Down0 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/2018_06_22_07h37/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Down_0.root","READ");
   TFile *bkgEstFile_Down1 = TFile::Open("root:://cmseos.fnal.gov//eos/uscms/store/user/rasharma/BackgroundEstimation/WWTree_CommonNtuple_For1and2Lepton_2018_05_15_04h15_WV_600_5TeV/2018_06_22_07h37/wjetmodel_Ex__WJets0_xww__sb_lo_ExpTail_auto_Down_2.root","READ");

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
   TTreeReaderValue<Float_t> pfMET_jes_up = {fReader, "pfMET_jes_up"};
   TTreeReaderValue<Float_t> pfMET_jes_dn = {fReader, "pfMET_jes_dn"};
   TTreeReaderValue<Float_t> pfMET_Phi = {fReader, "pfMET_Phi"};
   TTreeReaderValue<Float_t> pfMET_Corr = {fReader, "pfMET_Corr"};
   TTreeReaderValue<Float_t> pfMET_Corr_phi = {fReader, "pfMET_Corr_phi"};
   TTreeReaderValue<Float_t> pfMET_Corr_Cov00 = {fReader, "pfMET_Corr_Cov00"};
   TTreeReaderValue<Float_t> pfMET_Corr_Cov01 = {fReader, "pfMET_Corr_Cov01"};
   TTreeReaderValue<Float_t> pfMET_Corr_Cov11 = {fReader, "pfMET_Corr_Cov11"};
   TTreeReaderValue<Float_t> pfMET_Corr_jerup = {fReader, "pfMET_Corr_jerup"};
   TTreeReaderValue<Float_t> pfMET_Corr_jerdn = {fReader, "pfMET_Corr_jerdn"};
   TTreeReaderValue<Float_t> pfMET_Corr_jenup = {fReader, "pfMET_Corr_jenup"};
   TTreeReaderValue<Float_t> pfMET_Corr_jendn = {fReader, "pfMET_Corr_jendn"};
   TTreeReaderValue<Float_t> pfMET_Corr_uncup = {fReader, "pfMET_Corr_uncup"};
   TTreeReaderValue<Float_t> pfMET_Corr_uncdn = {fReader, "pfMET_Corr_uncdn"};
   TTreeReaderValue<Float_t> pfMET_Corr_jrsup = {fReader, "pfMET_Corr_jrsup"};
   TTreeReaderValue<Float_t> pfMET_Corr_jrsdn = {fReader, "pfMET_Corr_jrsdn"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijerup = {fReader, "pfMET_Corr_phijerup"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijerdn = {fReader, "pfMET_Corr_phijerdn"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijenup = {fReader, "pfMET_Corr_phijenup"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijendn = {fReader, "pfMET_Corr_phijendn"};
   TTreeReaderValue<Float_t> pfMET_Corr_phiuncup = {fReader, "pfMET_Corr_phiuncup"};
   TTreeReaderValue<Float_t> pfMET_Corr_phiuncdn = {fReader, "pfMET_Corr_phiuncdn"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijrsup = {fReader, "pfMET_Corr_phijrsup"};
   TTreeReaderValue<Float_t> pfMET_Corr_phijrsdn = {fReader, "pfMET_Corr_phijrsdn"};
   TTreeReaderValue<Float_t> nu_pz_type0 = {fReader, "nu_pz_type0"};
   TTreeReaderValue<Float_t> nu_pz_type2 = {fReader, "nu_pz_type2"};
   TTreeReaderValue<Float_t> nu_pz_run2 = {fReader, "nu_pz_run2"};
   TTreeReaderValue<Float_t> nu_pz_run2_oth = {fReader, "nu_pz_run2_oth"};
   TTreeReaderValue<Int_t> nu_pz_run2_type = {fReader, "nu_pz_run2_type"};
   TTreeReaderValue<Int_t> nu_pz_isre = {fReader, "nu_pz_isre"};
   TTreeReaderValue<Int_t> type = {fReader, "type"};
   TTreeReaderValue<Float_t> l_pt1 = {fReader, "l_pt1"};
   TTreeReaderValue<Float_t> l_eta1 = {fReader, "l_eta1"};
   TTreeReaderValue<Float_t> l_phi1 = {fReader, "l_phi1"};
   TTreeReaderValue<Float_t> l_e1 = {fReader, "l_e1"};
   TTreeReaderValue<Float_t> l_charge1 = {fReader, "l_charge1"};
   TTreeReaderValue<Float_t> l_iso1 = {fReader, "l_iso1"};
   TTreeReaderValue<Float_t> l_pt2 = {fReader, "l_pt2"};
   TTreeReaderValue<Float_t> l_eta2 = {fReader, "l_eta2"};
   TTreeReaderValue<Float_t> l_phi2 = {fReader, "l_phi2"};
   TTreeReaderValue<Float_t> l_e2 = {fReader, "l_e2"};
   TTreeReaderValue<Float_t> dilep_pt = {fReader, "dilep_pt"};
   TTreeReaderValue<Float_t> dilep_eta = {fReader, "dilep_eta"};
   TTreeReaderValue<Float_t> dilep_phi = {fReader, "dilep_phi"};
   TTreeReaderValue<Float_t> dilep_m = {fReader, "dilep_m"};
   TTreeReaderValue<Float_t> l_charge2 = {fReader, "l_charge2"};
   TTreeReaderValue<Float_t> l_iso2 = {fReader, "l_iso2"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt = {fReader, "ungroomed_AK8jet_pt"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jes_up = {fReader, "ungroomed_AK8jet_pt_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jes_dn = {fReader, "ungroomed_AK8jet_pt_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jer = {fReader, "ungroomed_AK8jet_pt_jer"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jer_up = {fReader, "ungroomed_AK8jet_pt_jer_up"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_pt_jer_dn = {fReader, "ungroomed_AK8jet_pt_jer_dn"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_eta = {fReader, "ungroomed_AK8jet_eta"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_eta_jes_up = {fReader, "ungroomed_AK8jet_eta_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_eta_jes_dn = {fReader, "ungroomed_AK8jet_eta_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_phi = {fReader, "ungroomed_AK8jet_phi"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_phi_jes_up = {fReader, "ungroomed_AK8jet_phi_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_phi_jes_dn = {fReader, "ungroomed_AK8jet_phi_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_e = {fReader, "ungroomed_AK8jet_e"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_charge = {fReader, "ungroomed_AK8jet_charge"};
   TTreeReaderValue<Float_t> AK8jet_mass = {fReader, "AK8jet_mass"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_mass_jes_up = {fReader, "ungroomed_AK8jet_mass_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_AK8jet_mass_jes_dn = {fReader, "ungroomed_AK8jet_mass_jes_dn"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr = {fReader, "AK8jet_mass_pr"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jes_up = {fReader, "AK8jet_mass_pr_jes_up"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jes_dn = {fReader, "AK8jet_mass_pr_jes_dn"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jer = {fReader, "AK8jet_mass_pr_jer"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jer_up = {fReader, "AK8jet_mass_pr_jer_up"};
   TTreeReaderValue<Float_t> AK8jet_mass_pr_jer_dn = {fReader, "AK8jet_mass_pr_jer_dn"};
   TTreeReaderValue<Float_t> AK8jet_mass_so = {fReader, "AK8jet_mass_so"};
   TTreeReaderValue<Float_t> AK8jet_pt_so = {fReader, "AK8jet_pt_so"};
   TTreeReaderValue<Float_t> AK8jet_mass_tr = {fReader, "AK8jet_mass_tr"};
   TTreeReaderValue<Float_t> AK8jet_mass_fi = {fReader, "AK8jet_mass_fi"};
   TTreeReaderValue<Float_t> AK8jet_tau2tau1 = {fReader, "AK8jet_tau2tau1"};
   TTreeReaderValue<Float_t> AK8jet_sj1_pt = {fReader, "AK8jet_sj1_pt"};
   TTreeReaderValue<Float_t> AK8jet_sj1_eta = {fReader, "AK8jet_sj1_eta"};
   TTreeReaderValue<Float_t> AK8jet_sj1_phi = {fReader, "AK8jet_sj1_phi"};
   TTreeReaderValue<Float_t> AK8jet_sj1_m = {fReader, "AK8jet_sj1_m"};
   TTreeReaderValue<Float_t> AK8jet_sj1_q = {fReader, "AK8jet_sj1_q"};
   TTreeReaderValue<Float_t> AK8jet_sj2_pt = {fReader, "AK8jet_sj2_pt"};
   TTreeReaderValue<Float_t> AK8jet_sj2_eta = {fReader, "AK8jet_sj2_eta"};
   TTreeReaderValue<Float_t> AK8jet_sj2_phi = {fReader, "AK8jet_sj2_phi"};
   TTreeReaderValue<Float_t> AK8jet_sj2_m = {fReader, "AK8jet_sj2_m"};
   TTreeReaderValue<Float_t> AK8jet_sj2_q = {fReader, "AK8jet_sj2_q"};
   TTreeReaderValue<Float_t> AK8_jetID_loose = {fReader, "AK8_jetID_loose"};
   TTreeReaderValue<Float_t> AK8jet_e3_b1 = {fReader, "AK8jet_e3_b1"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_b1 = {fReader, "AK8jet_e3_v1_b1"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_b1 = {fReader, "AK8jet_e3_v2_b1"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_b1 = {fReader, "AK8jet_e4_v1_b1"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_b1 = {fReader, "AK8jet_e4_v2_b1"};
   TTreeReaderValue<Float_t> AK8jet_e3_b2 = {fReader, "AK8jet_e3_b2"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_b2 = {fReader, "AK8jet_e3_v1_b2"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_b2 = {fReader, "AK8jet_e3_v2_b2"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_b2 = {fReader, "AK8jet_e4_v1_b2"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_b2 = {fReader, "AK8jet_e4_v2_b2"};
   TTreeReaderValue<Float_t> AK8jet_e2_sdb1 = {fReader, "AK8jet_e2_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e3_sdb1 = {fReader, "AK8jet_e3_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_sdb1 = {fReader, "AK8jet_e3_v1_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_sdb1 = {fReader, "AK8jet_e3_v2_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_sdb1 = {fReader, "AK8jet_e4_v1_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_sdb1 = {fReader, "AK8jet_e4_v2_sdb1"};
   TTreeReaderValue<Float_t> AK8jet_e2_sdb2 = {fReader, "AK8jet_e2_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e3_sdb2 = {fReader, "AK8jet_e3_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_sdb2 = {fReader, "AK8jet_e3_v1_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_sdb2 = {fReader, "AK8jet_e3_v2_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_sdb2 = {fReader, "AK8jet_e4_v1_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_sdb2 = {fReader, "AK8jet_e4_v2_sdb2"};
   TTreeReaderValue<Float_t> AK8jet_e2_sdb4 = {fReader, "AK8jet_e2_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e3_sdb4 = {fReader, "AK8jet_e3_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_sdb4 = {fReader, "AK8jet_e3_v1_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_sdb4 = {fReader, "AK8jet_e3_v2_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_sdb4 = {fReader, "AK8jet_e4_v1_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_sdb4 = {fReader, "AK8jet_e4_v2_sdb4"};
   TTreeReaderValue<Float_t> AK8jet_e2_sdb05 = {fReader, "AK8jet_e2_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_e3_sdb05 = {fReader, "AK8jet_e3_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_e3_v1_sdb05 = {fReader, "AK8jet_e3_v1_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_e3_v2_sdb05 = {fReader, "AK8jet_e3_v2_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_e4_v1_sdb05 = {fReader, "AK8jet_e4_v1_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_e4_v2_sdb05 = {fReader, "AK8jet_e4_v2_sdb05"};
   TTreeReaderValue<Float_t> AK8jet_qjet = {fReader, "AK8jet_qjet"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt = {fReader, "ungroomed_PuppiAK8_jet_pt"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jes_up = {fReader, "ungroomed_PuppiAK8_jet_pt_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jes_dn = {fReader, "ungroomed_PuppiAK8_jet_pt_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jer = {fReader, "ungroomed_PuppiAK8_jet_pt_jer"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jer_up = {fReader, "ungroomed_PuppiAK8_jet_pt_jer_up"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_pt_jer_dn = {fReader, "ungroomed_PuppiAK8_jet_pt_jer_dn"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_eta = {fReader, "ungroomed_PuppiAK8_jet_eta"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_eta_jes_up = {fReader, "ungroomed_PuppiAK8_jet_eta_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_eta_jes_dn = {fReader, "ungroomed_PuppiAK8_jet_eta_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_phi = {fReader, "ungroomed_PuppiAK8_jet_phi"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_phi_jes_up = {fReader, "ungroomed_PuppiAK8_jet_phi_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_phi_jes_dn = {fReader, "ungroomed_PuppiAK8_jet_phi_jes_dn"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_e = {fReader, "ungroomed_PuppiAK8_jet_e"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_charge = {fReader, "ungroomed_PuppiAK8_jet_charge"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass = {fReader, "PuppiAK8_jet_mass"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_mass_jes_up = {fReader, "ungroomed_PuppiAK8_jet_mass_jes_up"};
   TTreeReaderValue<Float_t> ungroomed_PuppiAK8_jet_mass_jes_dn = {fReader, "ungroomed_PuppiAK8_jet_mass_jes_dn"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr = {fReader, "PuppiAK8_jet_mass_pr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jes_up = {fReader, "PuppiAK8_jet_mass_pr_jes_up"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jes_dn = {fReader, "PuppiAK8_jet_mass_pr_jes_dn"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jer = {fReader, "PuppiAK8_jet_mass_pr_jer"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jer_up = {fReader, "PuppiAK8_jet_mass_pr_jer_up"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_pr_jer_dn = {fReader, "PuppiAK8_jet_mass_pr_jer_dn"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_so = {fReader, "PuppiAK8_jet_mass_so"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_so_corr = {fReader, "PuppiAK8_jet_mass_so_corr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_pt_so = {fReader, "PuppiAK8_jet_pt_so"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_tr = {fReader, "PuppiAK8_jet_mass_tr"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_mass_fi = {fReader, "PuppiAK8_jet_mass_fi"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_tau2tau1 = {fReader, "PuppiAK8_jet_tau2tau1"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj1_pt = {fReader, "PuppiAK8_jet_sj1_pt"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj1_eta = {fReader, "PuppiAK8_jet_sj1_eta"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj1_phi = {fReader, "PuppiAK8_jet_sj1_phi"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj1_m = {fReader, "PuppiAK8_jet_sj1_m"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj1_q = {fReader, "PuppiAK8_jet_sj1_q"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj2_pt = {fReader, "PuppiAK8_jet_sj2_pt"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj2_eta = {fReader, "PuppiAK8_jet_sj2_eta"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj2_phi = {fReader, "PuppiAK8_jet_sj2_phi"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj2_m = {fReader, "PuppiAK8_jet_sj2_m"};
   TTreeReaderValue<Float_t> PuppiAK8_jet_sj2_q = {fReader, "PuppiAK8_jet_sj2_q"};
   TTreeReaderValue<Float_t> PuppiAK8_jetID_loose = {fReader, "PuppiAK8_jetID_loose"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_b1 = {fReader, "PuppiAK8jet_e3_b1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_b1 = {fReader, "PuppiAK8jet_e3_v1_b1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_b1 = {fReader, "PuppiAK8jet_e3_v2_b1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_b1 = {fReader, "PuppiAK8jet_e4_v1_b1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_b1 = {fReader, "PuppiAK8jet_e4_v2_b1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_b2 = {fReader, "PuppiAK8jet_e3_b2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_b2 = {fReader, "PuppiAK8jet_e3_v1_b2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_b2 = {fReader, "PuppiAK8jet_e3_v2_b2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_b2 = {fReader, "PuppiAK8jet_e4_v1_b2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_b2 = {fReader, "PuppiAK8jet_e4_v2_b2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e2_sdb1 = {fReader, "PuppiAK8jet_e2_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_sdb1 = {fReader, "PuppiAK8jet_e3_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_sdb1 = {fReader, "PuppiAK8jet_e3_v1_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_sdb1 = {fReader, "PuppiAK8jet_e3_v2_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_sdb1 = {fReader, "PuppiAK8jet_e4_v1_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_sdb1 = {fReader, "PuppiAK8jet_e4_v2_sdb1"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e2_sdb2 = {fReader, "PuppiAK8jet_e2_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_sdb2 = {fReader, "PuppiAK8jet_e3_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_sdb2 = {fReader, "PuppiAK8jet_e3_v1_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_sdb2 = {fReader, "PuppiAK8jet_e3_v2_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_sdb2 = {fReader, "PuppiAK8jet_e4_v1_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_sdb2 = {fReader, "PuppiAK8jet_e4_v2_sdb2"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e2_sdb4 = {fReader, "PuppiAK8jet_e2_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_sdb4 = {fReader, "PuppiAK8jet_e3_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_sdb4 = {fReader, "PuppiAK8jet_e3_v1_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_sdb4 = {fReader, "PuppiAK8jet_e3_v2_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_sdb4 = {fReader, "PuppiAK8jet_e4_v1_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_sdb4 = {fReader, "PuppiAK8jet_e4_v2_sdb4"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e2_sdb05 = {fReader, "PuppiAK8jet_e2_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_sdb05 = {fReader, "PuppiAK8jet_e3_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v1_sdb05 = {fReader, "PuppiAK8jet_e3_v1_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e3_v2_sdb05 = {fReader, "PuppiAK8jet_e3_v2_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v1_sdb05 = {fReader, "PuppiAK8jet_e4_v1_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_e4_v2_sdb05 = {fReader, "PuppiAK8jet_e4_v2_sdb05"};
   TTreeReaderValue<Float_t> PuppiAK8jet_qjet = {fReader, "PuppiAK8jet_qjet"};
   TTreeReaderValue<Float_t> ttb_ungroomed_jet_pt = {fReader, "ttb_ungroomed_jet_pt"};
   TTreeReaderValue<Float_t> ttb_ungroomed_jet_eta = {fReader, "ttb_ungroomed_jet_eta"};
   TTreeReaderValue<Float_t> ttb_ungroomed_jet_phi = {fReader, "ttb_ungroomed_jet_phi"};
   TTreeReaderValue<Float_t> ttb_ungroomed_jet_e = {fReader, "ttb_ungroomed_jet_e"};
   TTreeReaderValue<Float_t> ttb_jet_mass_pr = {fReader, "ttb_jet_mass_pr"};
   TTreeReaderValue<Float_t> ttb_jet_mass_so = {fReader, "ttb_jet_mass_so"};
   TTreeReaderValue<Float_t> ttb_jet_pt_so = {fReader, "ttb_jet_pt_so"};
   TTreeReaderValue<Float_t> ttb_jet_mass_tr = {fReader, "ttb_jet_mass_tr"};
   TTreeReaderValue<Float_t> ttb_jet_mass_fi = {fReader, "ttb_jet_mass_fi"};
   TTreeReaderValue<Float_t> ttb_jet_tau2tau1 = {fReader, "ttb_jet_tau2tau1"};
   TTreeReaderValue<Float_t> isGen = {fReader, "isGen"};
   TTreeReaderValue<Float_t> lep_pt_gen = {fReader, "lep_pt_gen"};
   TTreeReaderValue<Float_t> lep_eta_gen = {fReader, "lep_eta_gen"};
   TTreeReaderValue<Float_t> W_pt_gen = {fReader, "W_pt_gen"};
   TTreeReaderValue<Float_t> W_pz_gen = {fReader, "W_pz_gen"};
   TTreeReaderValue<Float_t> W_rap_gen = {fReader, "W_rap_gen"};
   TTreeReaderValue<Float_t> nu_pz_gen = {fReader, "nu_pz_gen"};
   TTreeReaderValue<Float_t> nu_pt_gen = {fReader, "nu_pt_gen"};
   TTreeReaderValue<Float_t> nu_phi_gen = {fReader, "nu_phi_gen"};
   TTreeReaderValue<Float_t> nu_eta_gen = {fReader, "nu_eta_gen"};
   TTreeReaderValue<Float_t> hadW_pt_gen = {fReader, "hadW_pt_gen"};
   TTreeReaderValue<Float_t> hadW_eta_gen = {fReader, "hadW_eta_gen"};
   TTreeReaderValue<Float_t> hadW_phi_gen = {fReader, "hadW_phi_gen"};
   TTreeReaderValue<Float_t> hadW_e_gen = {fReader, "hadW_e_gen"};
   TTreeReaderValue<Float_t> hadW_m_gen = {fReader, "hadW_m_gen"};
   TTreeReaderValue<Float_t> lepW_pt_gen = {fReader, "lepW_pt_gen"};
   TTreeReaderValue<Float_t> lepW_eta_gen = {fReader, "lepW_eta_gen"};
   TTreeReaderValue<Float_t> lepW_phi_gen = {fReader, "lepW_phi_gen"};
   TTreeReaderValue<Float_t> lepW_e_gen = {fReader, "lepW_e_gen"};
   TTreeReaderValue<Float_t> lepW_m_gen = {fReader, "lepW_m_gen"};
   TTreeReaderValue<Float_t> WW_mass_gen = {fReader, "WW_mass_gen"};
   TTreeReaderValue<Float_t> WW_mT_gen = {fReader, "WW_mT_gen"};
   TTreeReaderValue<Float_t> WW_pT_gen = {fReader, "WW_pT_gen"};
   TTreeReaderValue<Float_t> AK8_pt_gen = {fReader, "AK8_pt_gen"};
   TTreeReaderValue<Float_t> AK8_eta_gen = {fReader, "AK8_eta_gen"};
   TTreeReaderValue<Float_t> AK8_phi_gen = {fReader, "AK8_phi_gen"};
   TTreeReaderValue<Float_t> AK8_e_gen = {fReader, "AK8_e_gen"};
   TTreeReaderValue<Float_t> AK8_mass_gen = {fReader, "AK8_mass_gen"};
   TTreeReaderValue<Float_t> AK8_pruned_mass_gen = {fReader, "AK8_pruned_mass_gen"};
   TTreeReaderValue<Float_t> AK8_softdrop_mass_gen = {fReader, "AK8_softdrop_mass_gen"};
   TTreeReaderValue<Float_t> AK8_softdrop_pt_gen = {fReader, "AK8_softdrop_pt_gen"};
   TTreeReaderValue<Float_t> v_pt_type0 = {fReader, "v_pt_type0"};
   TTreeReaderValue<Float_t> v_pt_type2 = {fReader, "v_pt_type2"};
   TTreeReaderValue<Float_t> v_pt_run2 = {fReader, "v_pt_run2"};
   TTreeReaderValue<Float_t> v_eta_type2 = {fReader, "v_eta_type2"};
   TTreeReaderValue<Float_t> v_eta_type0 = {fReader, "v_eta_type0"};
   TTreeReaderValue<Float_t> v_eta_run2 = {fReader, "v_eta_run2"};
   TTreeReaderValue<Float_t> v_phi = {fReader, "v_phi"};
   TTreeReaderValue<Float_t> v_mt_type2 = {fReader, "v_mt_type2"};
   TTreeReaderValue<Float_t> v_mt_type0 = {fReader, "v_mt_type0"};
   TTreeReaderValue<Float_t> v_mt_run2 = {fReader, "v_mt_run2"};
   TTreeReaderValue<Float_t> v_mass_type2 = {fReader, "v_mass_type2"};
   TTreeReaderValue<Float_t> v_mass_type0 = {fReader, "v_mass_type0"};
   TTreeReaderValue<Float_t> v_mass_run2 = {fReader, "v_mass_run2"};
   TTreeReaderValue<Float_t> v_pt_type0_jer_up = {fReader, "v_pt_type0_jer_up"};
   TTreeReaderValue<Float_t> v_eta_type0_jer_up = {fReader, "v_eta_type0_jer_up"};
   TTreeReaderValue<Float_t> v_mt_type0_jer_up = {fReader, "v_mt_type0_jer_up"};
   TTreeReaderValue<Float_t> v_mass_type0_jer_up = {fReader, "v_mass_type0_jer_up"};
   TTreeReaderValue<Float_t> v_pt_type0_jer_dn = {fReader, "v_pt_type0_jer_dn"};
   TTreeReaderValue<Float_t> v_eta_type0_jer_dn = {fReader, "v_eta_type0_jer_dn"};
   TTreeReaderValue<Float_t> v_mt_type0_jer_dn = {fReader, "v_mt_type0_jer_dn"};
   TTreeReaderValue<Float_t> v_mass_type0_jer_dn = {fReader, "v_mass_type0_jer_dn"};
   TTreeReaderValue<Float_t> v_pt_type0_jes_up = {fReader, "v_pt_type0_jes_up"};
   TTreeReaderValue<Float_t> v_eta_type0_jes_up = {fReader, "v_eta_type0_jes_up"};
   TTreeReaderValue<Float_t> v_mt_type0_jes_up = {fReader, "v_mt_type0_jes_up"};
   TTreeReaderValue<Float_t> v_mass_type0_jes_up = {fReader, "v_mass_type0_jes_up"};
   TTreeReaderValue<Float_t> v_pt_type0_jes_dn = {fReader, "v_pt_type0_jes_dn"};
   TTreeReaderValue<Float_t> v_eta_type0_jes_dn = {fReader, "v_eta_type0_jes_dn"};
   TTreeReaderValue<Float_t> v_mt_type0_jes_dn = {fReader, "v_mt_type0_jes_dn"};
   TTreeReaderValue<Float_t> v_mass_type0_jes_dn = {fReader, "v_mass_type0_jes_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type0 = {fReader, "mass_lvj_type0"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jes_up = {fReader, "mass_lvj_type0_met_jes_up"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jes_dn = {fReader, "mass_lvj_type0_met_jes_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jer_up = {fReader, "mass_lvj_type0_met_jer_up"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jer_dn = {fReader, "mass_lvj_type0_met_jer_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_jer = {fReader, "mass_lvj_type0_met_jer"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_PuppiAK8_jes_up = {fReader, "mass_lvj_type0_met_PuppiAK8_jes_up"};
   TTreeReaderValue<Float_t> mass_lvj_type0_met_PuppiAK8_jes_dn = {fReader, "mass_lvj_type0_met_PuppiAK8_jes_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type2 = {fReader, "mass_lvj_type2"};
   TTreeReaderValue<Float_t> mass_lvj_run2 = {fReader, "mass_lvj_run2"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8 = {fReader, "mass_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8_jes_up = {fReader, "mass_lvj_type0_PuppiAK8_jes_up"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8_jes_dn = {fReader, "mass_lvj_type0_PuppiAK8_jes_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8_jer_up = {fReader, "mass_lvj_type0_PuppiAK8_jer_up"};
   TTreeReaderValue<Float_t> mass_lvj_type0_PuppiAK8_jer_dn = {fReader, "mass_lvj_type0_PuppiAK8_jer_dn"};
   TTreeReaderValue<Float_t> mass_lvj_type2_PuppiAK8 = {fReader, "mass_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> mass_lvj_run2_PuppiAK8 = {fReader, "mass_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> mt_lvj_type0_PuppiAK8 = {fReader, "mt_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> mt_lvj_type2_PuppiAK8 = {fReader, "mt_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> mt_lvj_run2_PuppiAK8 = {fReader, "mt_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> pt_lvj_type0_PuppiAK8 = {fReader, "pt_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> pt_lvj_type2_PuppiAK8 = {fReader, "pt_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> pt_lvj_run2_PuppiAK8 = {fReader, "pt_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> eta_lvj_type0_PuppiAK8 = {fReader, "eta_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> eta_lvj_type2_PuppiAK8 = {fReader, "eta_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> eta_lvj_run2_PuppiAK8 = {fReader, "eta_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> rapidity_lvj_type0_PuppiAK8 = {fReader, "rapidity_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> rapidity_lvj_type2_PuppiAK8 = {fReader, "rapidity_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> rapidity_lvj_run2_PuppiAK8 = {fReader, "rapidity_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> phi_lvj_type0_PuppiAK8 = {fReader, "phi_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> phi_lvj_type2_PuppiAK8 = {fReader, "phi_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> phi_lvj_run2_PuppiAK8 = {fReader, "phi_lvj_run2_PuppiAK8"};
   TTreeReaderValue<Float_t> energy_lvj_type0_PuppiAK8 = {fReader, "energy_lvj_type0_PuppiAK8"};
   TTreeReaderValue<Float_t> energy_lvj_type2_PuppiAK8 = {fReader, "energy_lvj_type2_PuppiAK8"};
   TTreeReaderValue<Float_t> energy_lvj_run2_PuppiAK8 = {fReader, "energy_lvj_run2_PuppiAK8"};
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
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jes_up = {fReader, "vbf_maxpt_j1_pt_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jes_dn = {fReader, "vbf_maxpt_j1_pt_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jer = {fReader, "vbf_maxpt_j1_pt_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jer_up = {fReader, "vbf_maxpt_j1_pt_jer_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_pt_jer_dn = {fReader, "vbf_maxpt_j1_pt_jer_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta = {fReader, "vbf_maxpt_j1_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jes_up = {fReader, "vbf_maxpt_j1_eta_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jes_dn = {fReader, "vbf_maxpt_j1_eta_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jer = {fReader, "vbf_maxpt_j1_eta_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jer_up = {fReader, "vbf_maxpt_j1_eta_jer_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_eta_jer_dn = {fReader, "vbf_maxpt_j1_eta_jer_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_phi = {fReader, "vbf_maxpt_j1_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_phi_jes_up = {fReader, "vbf_maxpt_j1_phi_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_phi_jes_dn = {fReader, "vbf_maxpt_j1_phi_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_e = {fReader, "vbf_maxpt_j1_e"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_mass = {fReader, "vbf_maxpt_j1_mass"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_mass_jes_up = {fReader, "vbf_maxpt_j1_mass_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_mass_jes_dn = {fReader, "vbf_maxpt_j1_mass_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_bDiscriminatorCSV = {fReader, "vbf_maxpt_j1_bDiscriminatorCSV"};
   TTreeReaderValue<Float_t> vbf_maxpt_j1_charge = {fReader, "vbf_maxpt_j1_charge"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt = {fReader, "vbf_maxpt_j2_pt"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jes_up = {fReader, "vbf_maxpt_j2_pt_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jes_dn = {fReader, "vbf_maxpt_j2_pt_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jer = {fReader, "vbf_maxpt_j2_pt_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jer_up = {fReader, "vbf_maxpt_j2_pt_jer_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_pt_jer_dn = {fReader, "vbf_maxpt_j2_pt_jer_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta = {fReader, "vbf_maxpt_j2_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jes_up = {fReader, "vbf_maxpt_j2_eta_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jes_dn = {fReader, "vbf_maxpt_j2_eta_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jer = {fReader, "vbf_maxpt_j2_eta_jer"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jer_up = {fReader, "vbf_maxpt_j2_eta_jer_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_eta_jer_dn = {fReader, "vbf_maxpt_j2_eta_jer_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_phi = {fReader, "vbf_maxpt_j2_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_phi_jes_up = {fReader, "vbf_maxpt_j2_phi_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_phi_jes_dn = {fReader, "vbf_maxpt_j2_phi_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_e = {fReader, "vbf_maxpt_j2_e"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_mass = {fReader, "vbf_maxpt_j2_mass"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_mass_jes_up = {fReader, "vbf_maxpt_j2_mass_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_mass_jes_dn = {fReader, "vbf_maxpt_j2_mass_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_bDiscriminatorCSV = {fReader, "vbf_maxpt_j2_bDiscriminatorCSV"};
   TTreeReaderValue<Float_t> vbf_maxpt_j2_charge = {fReader, "vbf_maxpt_j2_charge"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_pt = {fReader, "vbf_maxpt_jj_pt"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_pt_jes_up = {fReader, "vbf_maxpt_jj_pt_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_pt_jes_dn = {fReader, "vbf_maxpt_jj_pt_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_eta = {fReader, "vbf_maxpt_jj_eta"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_phi = {fReader, "vbf_maxpt_jj_phi"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_m = {fReader, "vbf_maxpt_jj_m"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_m_jes_up = {fReader, "vbf_maxpt_jj_m_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_m_jes_dn = {fReader, "vbf_maxpt_jj_m_jes_dn"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_Deta = {fReader, "vbf_maxpt_jj_Deta"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_Deta_jes_up = {fReader, "vbf_maxpt_jj_Deta_jes_up"};
   TTreeReaderValue<Float_t> vbf_maxpt_jj_Deta_jes_dn = {fReader, "vbf_maxpt_jj_Deta_jes_dn"};
   TTreeReaderValue<Float_t> LepWEta = {fReader, "LepWEta"};
   TTreeReaderValue<Float_t> LepWRapidity = {fReader, "LepWRapidity"};
   TTreeReaderValue<Float_t> HadWEta = {fReader, "HadWEta"};
   TTreeReaderValue<Float_t> HadWRapidity = {fReader, "HadWRapidity"};
   TTreeReaderValue<Float_t> WWEta = {fReader, "WWEta"};
   TTreeReaderValue<Float_t> WWEta_PuppiAK8 = {fReader, "WWEta_PuppiAK8"};
   TTreeReaderValue<Float_t> WWRapidity = {fReader, "WWRapidity"};
   TTreeReaderValue<Float_t> WWRapidity_PuppiAK8 = {fReader, "WWRapidity_PuppiAK8"};
   TTreeReaderValue<Float_t> ZeppenfeldWH = {fReader, "ZeppenfeldWH"};
   TTreeReaderValue<Float_t> ZeppenfeldWH_jes_up = {fReader, "ZeppenfeldWH_jes_up"};
   TTreeReaderValue<Float_t> ZeppenfeldWH_jes_dn = {fReader, "ZeppenfeldWH_jes_dn"};
   TTreeReaderValue<Float_t> ZeppenfeldWHPuppi = {fReader, "ZeppenfeldWHPuppi"};
   TTreeReaderValue<Float_t> costheta1Puppi_type0 = {fReader, "costheta1Puppi_type0"};
   TTreeReaderValue<Float_t> costheta2Puppi_type0 = {fReader, "costheta2Puppi_type0"};
   TTreeReaderValue<Float_t> costhetastarPuppi_type0 = {fReader, "costhetastarPuppi_type0"};
   TTreeReaderValue<Float_t> phiPuppi_type0 = {fReader, "phiPuppi_type0"};
   TTreeReaderValue<Float_t> phi1Puppi_type0 = {fReader, "phi1Puppi_type0"};
   TTreeReaderValue<Float_t> VBSCentralityPuppi_type0 = {fReader, "VBSCentralityPuppi_type0"};
   TTreeReaderValue<Float_t> costheta1Puppi_type2 = {fReader, "costheta1Puppi_type2"};
   TTreeReaderValue<Float_t> costheta2Puppi_type2 = {fReader, "costheta2Puppi_type2"};
   TTreeReaderValue<Float_t> costhetastarPuppi_type2 = {fReader, "costhetastarPuppi_type2"};
   TTreeReaderValue<Float_t> phiPuppi_type2 = {fReader, "phiPuppi_type2"};
   TTreeReaderValue<Float_t> phi1Puppi_type2 = {fReader, "phi1Puppi_type2"};
   TTreeReaderValue<Float_t> VBSCentralityPuppi_type2 = {fReader, "VBSCentralityPuppi_type2"};
   TTreeReaderValue<Float_t> costheta1Puppi_run2 = {fReader, "costheta1Puppi_run2"};
   TTreeReaderValue<Float_t> costheta2Puppi_run2 = {fReader, "costheta2Puppi_run2"};
   TTreeReaderValue<Float_t> costhetastarPuppi_run2 = {fReader, "costhetastarPuppi_run2"};
   TTreeReaderValue<Float_t> phiPuppi_run2 = {fReader, "phiPuppi_run2"};
   TTreeReaderValue<Float_t> phi1Puppi_run2 = {fReader, "phi1Puppi_run2"};
   TTreeReaderValue<Float_t> VBSCentralityPuppi_run2 = {fReader, "VBSCentralityPuppi_run2"};
   TTreeReaderValue<Float_t> RpTPuppi_type0 = {fReader, "RpTPuppi_type0"};
   TTreeReaderValue<Float_t> ZeppenfeldWLPuppi_type0 = {fReader, "ZeppenfeldWLPuppi_type0"};
   TTreeReaderValue<Float_t> LeptonProjectionPuppi_type0 = {fReader, "LeptonProjectionPuppi_type0"};
   TTreeReaderValue<Float_t> RpTPuppi_type2 = {fReader, "RpTPuppi_type2"};
   TTreeReaderValue<Float_t> ZeppenfeldWLPuppi_type2 = {fReader, "ZeppenfeldWLPuppi_type2"};
   TTreeReaderValue<Float_t> LeptonProjectionPuppi_type2 = {fReader, "LeptonProjectionPuppi_type2"};
   TTreeReaderValue<Float_t> RpTPuppi_run2 = {fReader, "RpTPuppi_run2"};
   TTreeReaderValue<Float_t> ZeppenfeldWLPuppi_run2 = {fReader, "ZeppenfeldWLPuppi_run2"};
   TTreeReaderValue<Float_t> LeptonProjectionPuppi_run2 = {fReader, "LeptonProjectionPuppi_run2"};
   TTreeReaderValue<Float_t> PtBalancePuppi_type0 = {fReader, "PtBalancePuppi_type0"};
   TTreeReaderValue<Float_t> PtBalancePuppi_type2 = {fReader, "PtBalancePuppi_type2"};
   TTreeReaderValue<Float_t> PtBalancePuppi_run2 = {fReader, "PtBalancePuppi_run2"};
   TTreeReaderValue<Float_t> costheta1_type0 = {fReader, "costheta1_type0"};
   TTreeReaderValue<Float_t> costheta2_type0 = {fReader, "costheta2_type0"};
   TTreeReaderValue<Float_t> costhetastar_type0 = {fReader, "costhetastar_type0"};
   TTreeReaderValue<Float_t> phi_type0 = {fReader, "phi_type0"};
   TTreeReaderValue<Float_t> phi1_type0 = {fReader, "phi1_type0"};
   TTreeReaderValue<Float_t> VBSCentrality_type0 = {fReader, "VBSCentrality_type0"};
   TTreeReaderValue<Float_t> costheta1_type2 = {fReader, "costheta1_type2"};
   TTreeReaderValue<Float_t> costheta2_type2 = {fReader, "costheta2_type2"};
   TTreeReaderValue<Float_t> costhetastar_type2 = {fReader, "costhetastar_type2"};
   TTreeReaderValue<Float_t> phi_type2 = {fReader, "phi_type2"};
   TTreeReaderValue<Float_t> phi1_type2 = {fReader, "phi1_type2"};
   TTreeReaderValue<Float_t> VBSCentrality_type2 = {fReader, "VBSCentrality_type2"};
   TTreeReaderValue<Float_t> costheta1_run2 = {fReader, "costheta1_run2"};
   TTreeReaderValue<Float_t> costheta2_run2 = {fReader, "costheta2_run2"};
   TTreeReaderValue<Float_t> costhetastar_run2 = {fReader, "costhetastar_run2"};
   TTreeReaderValue<Float_t> phi_run2 = {fReader, "phi_run2"};
   TTreeReaderValue<Float_t> phi1_run2 = {fReader, "phi1_run2"};
   TTreeReaderValue<Float_t> VBSCentrality_run2 = {fReader, "VBSCentrality_run2"};
   TTreeReaderValue<Float_t> RpT_type0 = {fReader, "RpT_type0"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0 = {fReader, "ZeppenfeldWL_type0"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0_jes_up = {fReader, "ZeppenfeldWL_type0_jes_up"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0_jes_dn = {fReader, "ZeppenfeldWL_type0_jes_dn"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0_jer_up = {fReader, "ZeppenfeldWL_type0_jer_up"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type0_jer_dn = {fReader, "ZeppenfeldWL_type0_jer_dn"};
   TTreeReaderValue<Float_t> LeptonProjection_type0 = {fReader, "LeptonProjection_type0"};
   TTreeReaderValue<Float_t> RpT_type2 = {fReader, "RpT_type2"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_type2 = {fReader, "ZeppenfeldWL_type2"};
   TTreeReaderValue<Float_t> LeptonProjection_type2 = {fReader, "LeptonProjection_type2"};
   TTreeReaderValue<Float_t> RpT_run2 = {fReader, "RpT_run2"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_run2 = {fReader, "ZeppenfeldWL_run2"};
   TTreeReaderValue<Float_t> LeptonProjection_run2 = {fReader, "LeptonProjection_run2"};
   TTreeReaderValue<Float_t> PtBalance_type0 = {fReader, "PtBalance_type0"};
   TTreeReaderValue<Float_t> PtBalance_type0_jes_up = {fReader, "PtBalance_type0_jes_up"};
   TTreeReaderValue<Float_t> PtBalance_type0_jes_dn = {fReader, "PtBalance_type0_jes_dn"};
   TTreeReaderValue<Float_t> PtBalance_type0_jer_up = {fReader, "PtBalance_type0_jer_up"};
   TTreeReaderValue<Float_t> PtBalance_type0_jer_dn = {fReader, "PtBalance_type0_jer_dn"};
   TTreeReaderValue<Float_t> PtBalance_type2 = {fReader, "PtBalance_type2"};
   TTreeReaderValue<Float_t> PtBalance_run2 = {fReader, "PtBalance_run2"};
   TTreeReaderValue<Float_t> BosonCentrality_type0 = {fReader, "BosonCentrality_type0"};
   TTreeReaderValue<Float_t> BosonCentrality_type0_jes_up = {fReader, "BosonCentrality_type0_jes_up"};
   TTreeReaderValue<Float_t> BosonCentrality_type0_jes_dn = {fReader, "BosonCentrality_type0_jes_dn"};
   TTreeReaderValue<Float_t> BosonCentrality_type0_jer_up = {fReader, "BosonCentrality_type0_jer_up"};
   TTreeReaderValue<Float_t> BosonCentrality_type0_jer_dn = {fReader, "BosonCentrality_type0_jer_dn"};
   TTreeReaderValue<Float_t> BosonCentrality_type2 = {fReader, "BosonCentrality_type2"};
   TTreeReaderValue<Float_t> BosonCentrality_run2 = {fReader, "BosonCentrality_run2"};
   TTreeReaderValue<Float_t> BosonCentralityPuppi_type0 = {fReader, "BosonCentralityPuppi_type0"};
   TTreeReaderValue<Float_t> BosonCentralityPuppi_type2 = {fReader, "BosonCentralityPuppi_type2"};
   TTreeReaderValue<Float_t> BosonCentralityPuppi_run2 = {fReader, "BosonCentralityPuppi_run2"};
   TTreeReaderValue<Float_t> PtBalance_2Lep = {fReader, "PtBalance_2Lep"};
   TTreeReaderValue<Float_t> BosonCentrality_2Lep = {fReader, "BosonCentrality_2Lep"};
   TTreeReaderValue<Float_t> BosonCentrality_2Lep_jes_up = {fReader, "BosonCentrality_2Lep_jes_up"};
   TTreeReaderValue<Float_t> BosonCentrality_2Lep_jes_dn = {fReader, "BosonCentrality_2Lep_jes_dn"};
   TTreeReaderValue<Float_t> costheta1_2Lep = {fReader, "costheta1_2Lep"};
   TTreeReaderValue<Float_t> costheta2_2Lep = {fReader, "costheta2_2Lep"};
   TTreeReaderValue<Float_t> costhetastar_2Lep = {fReader, "costhetastar_2Lep"};
   TTreeReaderValue<Float_t> phi_2Lep = {fReader, "phi_2Lep"};
   TTreeReaderValue<Float_t> phi1_2Lep = {fReader, "phi1_2Lep"};
   TTreeReaderValue<Float_t> VBSCentrality_2Lep = {fReader, "VBSCentrality_2Lep"};
   TTreeReaderValue<Float_t> RpT_2Lep = {fReader, "RpT_2Lep"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_2Lep = {fReader, "ZeppenfeldWL_2Lep"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_2Lep_jes_up = {fReader, "ZeppenfeldWL_2Lep_jes_up"};
   TTreeReaderValue<Float_t> ZeppenfeldWL_2Lep_jes_dn = {fReader, "ZeppenfeldWL_2Lep_jes_dn"};
   TTreeReaderValue<Float_t> LeptonProjection_2Lep = {fReader, "LeptonProjection_2Lep"};




   TFile f("bkg_estimation.root", "RECREATE");

   // Create a histogram for the values we read.
   TH1F *hMC_Signal_4bin   = new TH1F("hMC_Signal_4bin",   "hMC_Signal_4bin;M_{WW};Events", 4, 600, 5000);
   TH1F *hMC_Signal_15bin  = new TH1F("hMC_Signal_15bin",  "hMC_Signal_15bin;M_{WW};Events", 15, 600, 5000);
   TH1F *hMC_Signal_88bin  = new TH1F("hMC_Signal_88bin",  "hMC_Signal_88bin;M_{WW};Events", 88, 600, 5000);

   TH1F *hSideBand_4bin  = new TH1F("hSideBand_4bin",  "hSideBand_4bin;M_{WW};Events", 4, 600, 5000);
   TH1F *hSideBand_15bin = new TH1F("hSideBand_15bin", "hSideBand_15bin;M_{WW};Events", 15, 600, 5000);
   TH1F *hSideBand_88bin = new TH1F("hSideBand_88bin", "hSideBand_88bin;M_{WW};Events", 88, 600, 5000);


   // Loop over all entries of the TTree or TChain.
   while (fReader.Next()) {
      // Just access the data as if myPx and btag were iterators (note the '*'
      // in front of them):
      if(!(*l_pt2<0 && *l_pt1>30)) continue;
      if(!(((*type==0)&&(abs(*l_eta1)<2.4))||((*type==1)&&((abs(*l_eta1)<2.5)&&!(abs(*l_eta1)>1.4442 && abs(*l_eta1)<1.566))))) continue;
      if(!(((*type==0)&&(*pfMET_Corr>50)) || ((*type==1)&&(*pfMET_Corr>80)))) continue;
      if(!((*ungroomed_PuppiAK8_jet_pt>200)&&(abs(*ungroomed_PuppiAK8_jet_eta)<2.4)&&(*PuppiAK8_jet_tau2tau1<0.55))) continue;
      //if(!((*PuppiAK8_jet_mass_so_corr>65) && (*PuppiAK8_jet_mass_so_corr<105))) continue;
      if(!(*nBTagJet_loose==0)) continue;
      if(!(*vbf_maxpt_jj_m>800)) continue;
      if(!(abs(*vbf_maxpt_j2_eta-*vbf_maxpt_j1_eta)>4.0)) continue;
      if(!((*vbf_maxpt_j1_pt>30) && (*vbf_maxpt_j2_pt>30))) continue;
      if(!(*mass_lvj_type0_PuppiAK8>600)) continue;
      if(!(*BosonCentrality_type0>1.0)) continue;
      if(!((abs(*ZeppenfeldWL_type0)/abs(*vbf_maxpt_j2_eta-*vbf_maxpt_j1_eta))<0.3)) continue;
      if(!((abs(*ZeppenfeldWH)/abs(*vbf_maxpt_j2_eta-*vbf_maxpt_j1_eta))<0.3)) continue;
      

      // Fill histogram for signal region
      if ((*PuppiAK8_jet_mass_so_corr>65) && (*PuppiAK8_jet_mass_so_corr<105))
      {
      	hMC_Signal_4bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      	hMC_Signal_15bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      	hMC_Signal_88bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      }

      // Fill histogram for side-band region
      if ((*PuppiAK8_jet_mass_so_corr<65) || (*PuppiAK8_jet_mass_so_corr>105))
      {
      	hSideBand_4bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      	hSideBand_15bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      	hSideBand_88bin->Fill(*mass_lvj_type0_PuppiAK8,((*wSampleWeight)*(35867.06)*(*btag0Wgt)*(*genWeight)*(*trig_eff_Weight)*(*id_eff_Weight)*(*pu_Weight)));
      }
      //cout<<*mass_lvj_type0_PuppiAK8<<"\t"<<*BosonCentrality_type0<<endl;
   }

   // include overflow bin
   hMC_Signal_4bin->SetBinContent(4,hMC_Signal_4bin->GetBinContent(4)+hMC_Signal_4bin->GetBinContent(5));
   hMC_Signal_15bin->SetBinContent(4,hMC_Signal_15bin->GetBinContent(4)+hMC_Signal_15bin->GetBinContent(5));
   hMC_Signal_88bin->SetBinContent(4,hMC_Signal_88bin->GetBinContent(4)+hMC_Signal_88bin->GetBinContent(5));
   hSideBand_4bin->SetBinContent(4,hSideBand_4bin->GetBinContent(4)+hSideBand_4bin->GetBinContent(5));
   hSideBand_15bin->SetBinContent(4,hSideBand_15bin->GetBinContent(4)+hSideBand_15bin->GetBinContent(5));
   hSideBand_88bin->SetBinContent(4,hSideBand_88bin->GetBinContent(4)+hSideBand_88bin->GetBinContent(5));

   // Calculate alpha
   TH1F* alpha = (TH1F*)hMC_Signal_15bin->Clone();
   alpha->Divide(hSideBand_15bin);
   alpha->SetMinimum(-1.0);
   alpha->SetMaximum(4.0);
   alpha->SetName("alpha");
   alpha->GetXaxis()->SetTitle("M_{WW} (GeV)");
   alpha->GetYaxis()->SetTitle("alpha = #frac{N^{MC,V+jets}_{signal}}{N^{MC,V+jets}_{side-band}}");

   // Fit alpha using polynomial of order 1
   gStyle->SetOptFit(1);
   TF1* f1 = new TF1("f1","pol1",600,5000);
   f1->SetLineColor(2);
   TFitResultPtr fitRes = alpha->Fit("f1","S");

   TCanvas* c1 = new TCanvas("AlphaSigmaBand","canvas", 1000, 900);

   // Get sigma band for alpha
   TH1F *hConf1s = (TH1F*)alpha->Clone("hConf1s");
   TH1F *hConf2s = (TH1F*)alpha->Clone("hConf2s");
   hConf1s->Reset();
   hConf2s->Reset();
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hConf1s, 0.68);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hConf2s, 0.95);
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

   TLegend* leg = new TLegend(0.10,0.95,0.55,0.75);
   leg->SetNColumns(2);
   leg->AddEntry(alpha, "MC", "lep");
   leg->AddEntry(f1, "Fit", "l");
   leg->AddEntry(hConf1s, "1 #sigma band","f");
   leg->AddEntry(hConf2s, "2 #sigma band","f");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.55, 0.75, 0.90,0.95, "brNDC");
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
  //	PART- 2: Alpha function variation by +/- 1 sigma
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
   c1->SetName("AlphaParVariation");

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
  hCentral->SetName("Alpha_Central");
  hCentral->SetTitle("");

  hCentral->SetLineColor(1);
  hCentral->SetMarkerSize(1.5);
  hCentral->SetMarkerStyle(5);
  alpha->Draw();
  hCentral->Draw("same");
  leg->AddEntry(alpha,"Alpha MC","lep");
  leg->AddEntry(hCentral,"Alpha Nominal","l");
  //Now conver systematic fits into histos
  std::vector<TH1F*> histos;
  TString catname = "Alpha_Hist_Systematic_fits_";
  for (int ifun=0; ifun<systFunctions.size(); ++ifun){
    TF1* thisFunc = systFunctions[ifun];
    TH1F* thisHisto = convertTF1toTH1F(thisFunc, hMC_Signal_88bin);
    TString thisHistoName = catname+names[ifun];
    thisHisto->SetName(thisHistoName);
    thisHisto->SetTitle(thisHistoName);
    //thisHisto->Write();
    thisHisto->SetLineColor(ifun+1);
    thisHisto->Draw("same");
    histos.push_back(thisHisto);
    leg->AddEntry(thisHisto,thisHistoName);
  }
  leg->Draw();	pt->Draw();
  c1->Draw();
  c1->Write();
  leg->Clear();
  c1->Clear();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART - 3: Get signal histogram from alpha multiplication
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //File is opened at start of program so that TFile will not interfer with the output TFile

   TLegend* leg2 = new TLegend(0.35,0.90,0.90,0.70);
   leg2->SetNColumns(2);
  c1->SetName("Comp_Wjet_SideBand_MC_Corr");
   //TH1F* alpha = (TH1F*)hMC_Signal_15bin->Clone();
  TH1F* Wjet_Corr_Hist = (TH1F*)bkgEstFile->Get("rrv_mass_lvj__rrv_mass_lvj");
  Wjet_Corr_Hist->Scale(Wjet_Normalization_FromBkgEstimation);
  Wjet_Corr_Hist->SetStats(0);
  Wjet_Corr_Hist->SetLineColor(2);
  Wjet_Corr_Hist->SetMarkerColor(2);
  hSideBand_88bin->SetStats(0);

  Wjet_Corr_Hist->SetName("Wjet_SideBandRegion_Corr_Hist_From_Data");
  Wjet_Corr_Hist->Draw("hist");
  hSideBand_88bin->Draw("E1 same");
  leg2->AddEntry(Wjet_Corr_Hist,"SideBand Region (Corr)");
  leg2->AddEntry(hSideBand_88bin,"SideBand Region (MC)","lep");

  leg2->Draw();
  c1->Draw();
  c1->Write();
  c1->SetName("Comp_Wjet_SideBand_MC_Corr_Log");
  c1->SetLogy(1);
  c1->Write();
  c1->SetLogy(0);
  c1->Clear();
  leg2->Clear();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART - 3: Get signal histogram from alpha multiplication
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////

  c1->SetName("Vjet_Corr_SignalRegion_Hist_AlphaSystematic_fits");
  TH1F* hCorr_Signal_Central = getWjetSignalRegion_usingAlpha(Wjet_Corr_Hist,hCentral);
  hCorr_Signal_Central->SetName("Vjet_SRHist_Systematic_fits_Nominal");
  hCorr_Signal_Central->SetTitle("Corr Vjet Signal Region (Nominal)");

  hCorr_Signal_Central->SetLineColor(1);
  hCorr_Signal_Central->SetMarkerColor(1);
  hCorr_Signal_Central->SetMarkerSize(1.5);
  hCorr_Signal_Central->SetMarkerStyle(5);
  hMC_Signal_88bin->Draw();
  hCorr_Signal_Central->Draw("same");
  leg2->AddEntry(hMC_Signal_88bin,"Vjet MC","lep");
  leg2->AddEntry(hCorr_Signal_Central,"Vjet Corr SR Nominal","l");
  //Now conver systematic fits into histos
  std::vector<TH1F*> histos_Vjet_SR;
  histos_Vjet_SR.push_back(hCorr_Signal_Central);
  catname = "Vjet_SRHist_Systematic_fits_";
  for (int ifun=0; ifun<histos.size(); ++ifun){
    TH1F* thisFunc = histos[ifun];
    TH1F* thisHisto = getWjetSignalRegion_usingAlpha(Wjet_Corr_Hist, thisFunc);
    TString thisHistoName = catname+names[ifun];
    thisHisto->SetName(thisHistoName);
    thisHisto->SetTitle(thisHistoName);
    //thisHisto->Write();
    thisHisto->SetLineColor(ifun+1);
    thisHisto->SetMarkerColor(ifun+1);
    histos_Vjet_SR.push_back(thisHisto);
    thisHisto->Draw("same");
    leg2->AddEntry(thisHisto,thisHistoName);
  }
  leg2->Draw();
  c1->Draw();
  c1->Write();
  c1->SetName("Vjet_Corr_SignalRegion_Hist_AlphaSystematic_fits_Log");
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
  
  c1->SetName("Vjet_Corr_AlphaSystematics_SR_4bins");
  catname = "Vjet_Corr_AlphaSystematics_SR_4bins_";
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
    thisHisto->SetLineColor(ifun+1);
    if (ifun == 0) thisHisto->Draw();
    else thisHisto->Draw("same");
    leg2->AddEntry(thisHisto,thisHistoName);
  }
  leg2->Draw();	
  c1->Draw();
  c1->Write();
  c1->SetName("Vjet_Corr_AlphaSystematics_SR_4bins_Log");
  c1->SetLogy(1);
  c1->Write();
  leg2->Clear();
  c1->Clear();
  c1->SetLogy(0);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART - 6: Get Wjet fit systematics
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::vector<TH1F*> histos_WjetSyst;
  //histos_WjetSyst.push_back(Wjet_Corr_Hist)
  c1->SetName("Vjet_Corr_WjetFitSystematics_SideBandRegion");
  TH1F* Wjet_Corr_Hist_Up0 = (TH1F*)bkgEstFile_Up0->Get("rrv_mass_lvj__rrv_mass_lvj");
  Wjet_Corr_Hist_Up0->SetName("Wjet_SideBandRegion_Corr_Hist_From_Data_Par0Up");
  Wjet_Corr_Hist_Up0->Scale(Wjet_Normalization_FromBkgEstimation);
  Wjet_Corr_Hist_Up0->SetStats(0);
  Wjet_Corr_Hist_Up0->SetLineColor(3);
  Wjet_Corr_Hist_Up0->SetMarkerColor(3);
  histos_WjetSyst.push_back(Wjet_Corr_Hist_Up0);

  TH1F* Wjet_Corr_Hist_Up1 = (TH1F*)bkgEstFile_Up1->Get("rrv_mass_lvj__rrv_mass_lvj");
  Wjet_Corr_Hist_Up1->SetName("Wjet_SideBandRegion_Corr_Hist_From_Data_Par1Up");
  Wjet_Corr_Hist_Up1->Scale(Wjet_Normalization_FromBkgEstimation);
  Wjet_Corr_Hist_Up1->SetStats(0);
  Wjet_Corr_Hist_Up1->SetLineColor(4);
  Wjet_Corr_Hist_Up1->SetMarkerColor(4);
  histos_WjetSyst.push_back(Wjet_Corr_Hist_Up1);

  TH1F* Wjet_Corr_Hist_Down0 = (TH1F*)bkgEstFile_Down0->Get("rrv_mass_lvj__rrv_mass_lvj");
  Wjet_Corr_Hist_Down0->SetName("Wjet_SideBandRegion_Corr_Hist_From_Data_Par0Down");
  Wjet_Corr_Hist_Down0->Scale(Wjet_Normalization_FromBkgEstimation);
  Wjet_Corr_Hist_Down0->SetStats(0);
  Wjet_Corr_Hist_Down0->SetLineColor(5);
  Wjet_Corr_Hist_Down0->SetMarkerColor(5);
  histos_WjetSyst.push_back(Wjet_Corr_Hist_Down0);

  TH1F* Wjet_Corr_Hist_Down1 = (TH1F*)bkgEstFile_Down1->Get("rrv_mass_lvj__rrv_mass_lvj");
  Wjet_Corr_Hist_Down1->SetName("Wjet_SideBandRegion_Corr_Hist_From_Data_Par1Down");
  Wjet_Corr_Hist_Down1->Scale(Wjet_Normalization_FromBkgEstimation);
  Wjet_Corr_Hist_Down1->SetStats(0);
  Wjet_Corr_Hist_Down1->SetLineColor(6);
  Wjet_Corr_Hist_Down1->SetMarkerColor(6);
  histos_WjetSyst.push_back(Wjet_Corr_Hist_Down1);

  hSideBand_88bin->SetTitle("Vjet SideBand Corr (Wjet Fit systematics)");
  hSideBand_88bin->Draw("E1");
  Wjet_Corr_Hist->Draw("hist same");
  Wjet_Corr_Hist_Up0->Draw("hist same");
  Wjet_Corr_Hist_Up1->Draw("hist same");
  Wjet_Corr_Hist_Down0->Draw("hist same");
  Wjet_Corr_Hist_Down1->Draw("hist same");


  leg2->AddEntry(hSideBand_88bin,"MC");
  leg2->AddEntry(Wjet_Corr_Hist,"Vjet Corr Nominal");
  leg2->AddEntry(Wjet_Corr_Hist_Up0,"Vjet Corr Par0 Up");
  leg2->AddEntry(Wjet_Corr_Hist_Up1,"Vjet Corr Par1 Up");
  leg2->AddEntry(Wjet_Corr_Hist_Down0,"Vjet Corr Par0 Down");
  leg2->AddEntry(Wjet_Corr_Hist_Down1,"Vjet Corr Par1 Down");
  leg2->Draw();

  c1->Draw();
  c1->Write();
  c1->SetName("Vjet_Corr_WjetFitSystematics_SideBandRegion_Log");
  c1->SetLogy(1);
  c1->Write();
  leg2->Clear();
  c1->Clear();
  c1->SetLogy(0);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART - 6: Get Wjet fit systematics (in signal region by multiplying it with alpha)
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  c1->SetName("Vjet_Corr_WjetFitSystematics_SignalRegion");
  //TH1F* hCorr_Signal_Central = getWjetSignalRegion_usingAlpha(Wjet_Corr_Hist,hCentral);
  //hCorr_Signal_Central->SetName("Corr_Vjet_Signal_Region_Nominal");
  //hCorr_Signal_Central->SetTitle("Corr Vjet Signal Region (Nominal)");
  //hCorr_Signal_Central->SetLineColor(1);
  //hCorr_Signal_Central->SetMarkerSize(1.5);
  //hCorr_Signal_Central->SetMarkerStyle(5);
  hMC_Signal_88bin->Draw("E1");
  hCorr_Signal_Central->Draw("hist same");
  leg2->AddEntry(hMC_Signal_88bin,"MC");
  leg2->AddEntry(hCorr_Signal_Central,"Vjet Corr Nominal");


  std::vector<TH1F*> histos_WjetSyst_SR;
  histos_WjetSyst_SR.push_back(hCorr_Signal_Central);
  catname = "Vjet_Corr_WjetFitSystematics_SignalRegion_";
  for (int ifun=0; ifun<histos_WjetSyst.size(); ++ifun){
    TH1F* thisFunc = histos_WjetSyst[ifun];
    TH1F* thisHisto = getWjetSignalRegion_usingAlpha(thisFunc, hCentral);
    TString thisHistoName;
    thisHistoName = catname+names[ifun-1];
    thisHisto->SetName(thisHistoName);
    //thisHisto->SetTitle(thisHistoName);
    thisHisto->SetTitle("");
    thisHisto->SetStats(0);
    thisHisto->SetLineColor(ifun+2);
    histos_WjetSyst_SR.push_back(thisHisto);
    thisHisto->Draw("same");
    leg2->AddEntry(thisHisto,thisHistoName);
  }
  leg2->Draw();

  c1->Draw();
  c1->Write();
  c1->SetName("Vjet_Corr_WjetFitSystematics_SignalRegion_Log");
  c1->SetLogy(1);
  c1->Write();
  leg2->Clear();
  c1->Clear();
  c1->SetLogy(0);

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	PART - 6: Get Wjet fit systematics (in signal region by multiplying it with alpha)
  //	Change to 4 bins
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////

  c1->SetName("Vjet_Corr_WjetFitSystematics_SignalRegion_4bins");
  catname = "Vjet_Corr_WjetFitSystematics_SignalRegion_4bins_";
  for (int ifun=0; ifun<histos_WjetSyst_SR.size(); ++ifun){
    TH1F* thisFunc = histos_WjetSyst_SR[ifun];
    TH1F* thisHisto = ResetTo4bins(thisFunc, 4, 600, 2500);
    TString thisHistoName;
    if (ifun == 0) thisHistoName = catname+"Nominal";
    else thisHistoName = catname+names[ifun-1];
    thisHisto->SetName(thisHistoName);
    //thisHisto->SetTitle(thisHistoName);
    thisHisto->SetTitle("");
    thisHisto->SetStats(0);
    thisHisto->SetLineColor(ifun+1);
    if (ifun == 0) thisHisto->Draw();
    else thisHisto->Draw("same");
    leg2->AddEntry(thisHisto,thisHistoName);
  }
  leg2->Draw();	
  c1->Draw();
  c1->Write();
  c1->SetName("Vjet_Corr_WjetFitSystematics_SignalRegion_4bins_Log");
  c1->SetLogy(1);
  c1->Write();
  leg2->Clear();
  c1->Clear();
  c1->SetLogy(0);

   f1->Write();   f.Write();   f.Close();
}
