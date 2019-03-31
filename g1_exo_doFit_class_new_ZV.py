#python g1_exo_doFit_class.py -b -c mu > test.log
#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle
from array import array

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, TStopwatch

msgservice = ROOT.RooMsgService.instance()
msgservice.setGlobalKillBelow(RooFit.FATAL)

############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="em")
parser.add_option('--sample', action="store",type="string",dest="sample",default="Signal_aQGC")
parser.add_option('--mass', action="store",type="float",dest="mass",default="600")
parser.add_option('-s','--simple', action='store', dest='simple', default=False, help='pre-limit in simple mode')
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--type', action="store",type="string",dest="type",default="vbf")
parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="PuppiAK8_jet_mass_so_corr")
parser.add_option('--interpolate', action="store_true",dest="interpolate",default=False)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")

from ROOT import draw_error_band, draw_error_band2, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

fitresultsmj    = []
fitresultsmlvj  = []
fitresultsfinal = []
fit_DecoResults = []


class doFit_wj_and_wlvj:

    def __init__(self, in_channel,in_signal_sample, jetalgo, in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700, in_mj_min=30, in_mj_max=140, in_mlvj_min=400., in_mlvj_max=1400., fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", interpolate=False, input_workspace=None):

        tdrstyle.setTDRStyle()
        #TGaxis.SetMaxDigits(3)
	
        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ### set the channel type --> electron or muon
        self.channel=in_channel;
        self.leg = TLegend(); 
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;
	self.jetalgo = jetalgo
        self.interpolate = interpolate
                
        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuissance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj=5.;

        ### Set the binning for mlvj plots as a function of the model
        self.BinWidth_mlvj=50.;
            
        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor=1

        ## correct the binning of mj 
        self.BinWidth_mj=self.BinWidth_mj/self.narrow_factor;
        nbins_mj=int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max=in_mj_min+nbins_mj*self.BinWidth_mj;
                   
        ## correct the binning of mlvj 
        self.BinWidth_mlvj=self.BinWidth_mlvj/self.narrow_factor;
        nbins_mlvj=int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        in_mlvj_max=in_mlvj_min+nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
	varname = "Softdrop jet mass (GeV)"
	if self.jetalgo == "jet_mass_pr": varname = "Pruned jet mass (GeV)"
	if self.jetalgo == "PuppiAK8_jet_mass_pr": varname = "Pruned jet mass (GeV)"
	if self.jetalgo == "PuppiAK8_jet_mass_so_corr": varname = "Softdrop jet mass (GeV)"

        rrv_mass_j = RooRealVar("rrv_mass_j",varname,(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);

        ## define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","m_{ZV} (GeV)",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV");
        rrv_mass_lvj.setBins(nbins_mlvj);

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj=fit_model;
        self.MODEL_4_mlvj_alter=fit_model_alter;

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_");
        else:
            self.workspace4fit_ = input_workspace;
        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);

        #prepare workspace for unbin-Limit -> just fo the stuff on which running the limit 
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_");

        #define sidebands
        self.mj_sideband_lo_min = in_mj_min;
        self.mj_sideband_lo_max = 65
        self.mj_sideband_hi_min = 105
        self.mj_sideband_hi_max = in_mj_max;
        self.mj_signal_min = 65
        self.mj_signal_max = 105
        if (options.category).find('W') != -1:
         self.mj_signal_min = 65
         self.mj_signal_max = 85	        
	elif (options.category).find('Z') != -1:
         self.mj_signal_min = 85
         self.mj_signal_max = 105
	 	    
        ## zone definition in the jet mass 
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

        ## signal region definition in the mlvj variable in case of counting limit
        self.mlvj_signal_min = in_mlvj_signal_region_min
        self.mlvj_signal_max = in_mlvj_signal_region_max
        rrv_mass_lvj.setRange("signal_region",self.mlvj_signal_min,self.mlvj_signal_max);
        rrv_mass_lvj.setRange("high_mass",2500,in_mlvj_max);

        #prepare the data and mc files --> set the working directory and the files name
	self.file_Directory="/store/user/rasharma/SecondStep/WWTree_After_CWR/2019_03_28_16h05/HaddedFiles/Hadds_for_BkgEstimation/";
	#self.file_Directory="Ntuples2/";
                 
        #prepare background data and signal samples            
        self.signal_sample=in_signal_sample;

        self.file_data = ("WWTree_data_golden.root");#keep blind!!!!
#        self.file_data = ("WWTree_pseudodataS.root");#fake data
#        self.file_data = ("WWTree_pseudodata.root");#fake data
        self.file_signal     = ("WWTree_%s.root"%(self.signal_sample));
        self.file_WJets0_mc  = ("WWTree_VJets.root");
        #self.file_WJets0_mc  = ("WWTree_WJets.root");
        #self.file_VV_mc      = ("WWTree_VV.root");# WW+WZ
        self.file_VV_mc      = ("WWTree_VV_EWK_QCD.root");# WW+WZ
        self.file_TTbar_mc   = ("WWTree_TTbar.root");
        #self.file_TTbar_mc   = ("WWTree_VV.root");
        self.file_STop_mc    = ("WWTree_STop.root");

        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;

        self.categoryID=-1;
        if (self.channel=="el" or self.channel=="em") and self.wtagger_label.find("HP") != -1:
           self.categoryID=0;
           if self.wtagger_label.find("0") != -1: self.categoryID=6;
           elif self.wtagger_label.find("1") != -1: self.categoryID=8;
           elif self.wtagger_label.find("2") != -1: self.categoryID=10;
	if self.channel=="mu" and self.wtagger_label.find("HP") != -1:
	   self.categoryID=1;
           if self.wtagger_label.find("0") != -1: self.categoryID=7;
           elif self.wtagger_label.find("1") != -1: self.categoryID=9;
           elif self.wtagger_label.find("2") != -1: self.categoryID=11;
	if (self.channel=="el" or self.channel=="em") and self.wtagger_label.find("LP") != -1:
	   self.categoryID=2;
           if self.wtagger_label.find("0") != -1: self.categoryID=12;
           elif self.wtagger_label.find("1") != -1: self.categoryID=14;
           elif self.wtagger_label.find("2") != -1: self.categoryID=16;
	if self.channel=="mu" and self.wtagger_label.find("LP") != -1:
	   self.categoryID=3;
           if self.wtagger_label.find("0") != -1: self.categoryID=13;
           elif self.wtagger_label.find("1") != -1: self.categoryID=15;
           elif self.wtagger_label.find("2") != -1: self.categoryID=17;

        if self.channel=="mu" and self.wtagger_label.find("HP") != -1:
            self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.000);
            self.rrv_wtagger_eff_reweight_forT.setError(0.000);
            self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.000);
            self.rrv_wtagger_eff_reweight_forV.setError(0.000);
        if (self.channel=="el" or self.channel=="em") and self.wtagger_label.find("HP") != -1:
            self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.000);
            self.rrv_wtagger_eff_reweight_forT.setError(0.000);
            self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.000);
            self.rrv_wtagger_eff_reweight_forV.setError(0.000);
        if self.channel=="mu" and self.wtagger_label.find("LP") != -1:
            self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.000);
            self.rrv_wtagger_eff_reweight_forT.setError(0.000);
            self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.000);
            self.rrv_wtagger_eff_reweight_forV.setError(0.000);
        if (self.channel=="el" or self.channel=="em") and self.wtagger_label.find("LP") != -1:
            self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.000);
            self.rrv_wtagger_eff_reweight_forT.setError(0.000);
            self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.000);
            self.rrv_wtagger_eff_reweight_forV.setError(0.000);

        print "wtagger efficiency correction for Top sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forT.getVal(), self.rrv_wtagger_eff_reweight_forT.getError());
        print "wtagger efficiency correction for V sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forV.getVal(), self.rrv_wtagger_eff_reweight_forV.getError());
        
        self.mean_shift = 0.000
        self.sigma_scale=1.000
        
	self.plotsDir = 'plots_%s_%s' %(self.channel,self.wtagger_label)
        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file
        if not os.path.isdir("cards_%s_%s"%(self.channel,self.wtagger_label)):
            os.system("mkdir cards_%s_%s"%(self.channel,self.wtagger_label));
        self.rlt_DIR="cards_%s_%s/"%(self.channel,self.wtagger_label)

        ## extra text file
        self.file_rlt_txt = self.rlt_DIR+"other_wwlvj_%s_%s_%s.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"wwlvj_%s_%s_%s_workspace.root"%(self.signal_sample,self.channel,self.wtagger_label)
        ## datacard for the ubninned limit
        self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_%s_%s_unbin.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        ## workspace for the binned limit
        self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_%s_%s_counting.txt"%(self.signal_sample,self.channel,self.wtagger_label)
        
        self.file_out=open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        #self.file_out.close()
        self.file_out=open(self.file_rlt_txt,"a+");

        ## color palette 
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : ROOT.TColor.GetColor(222,90,106),
	    #'WJets' : 625,
            #'VV' : 607,
	    'VV' : ROOT.TColor.GetColor(250,202,255),
            #'STop' : 854,
            'STop' : 592,	# CWR comment
            'TTbar' : 592,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband=-1;
        self.datadriven_alpha_WJets_unbin=-1;
        self.datadriven_alpha_WJets_counting=-1;

        ### uncertainty for datacard
#        self.lumi_uncertainty      = 0.046;
        self.lumi_uncertainty      = 0.000; # updated feb 22 for Moriond
        #self.XS_STop_uncertainty = 0.30 ;
        #self.XS_VV_uncertainty   = 0.25 ;
        #self.XS_TTbar_NLO_uncertainty = 0.063 ;# from AN-12-368 table8
        self.XS_STop_NLO_uncertainty  = 0.00 ;# from AN-12-368 table8
        self.XS_VV_NLO_uncertainty    = 0.00 ;# from AN-12-368 table8
#        self.XS_sig_NNLO_uncertainty  = 0.00 ;

        #el and mu trigger and eff uncertainty, AN2012_368_v5 12.3
        self.lep_trigger_uncertainty = 0.00;
        self.lep_eff_uncertainty     = 0.00;
        #b tag scale uncertainty
        self.btag_scale_uncertainty  = 0.000;
        self.signal_btag_uncertainty = 0.000;
      
        if self.channel == "mu":
         self.signal_lepton_energy_scale_uncertainty = 0.000 ;
         self.signal_lepton_energy_res_uncertainty   = 0.000 ;
         self.signal_jet_energy_res_uncertainty      = 0.000 ;
        elif self.channel == "el":
         self.signal_lepton_energy_scale_uncertainty = 0.000 ;
         self.signal_lepton_energy_res_uncertainty   = 0.000 ;
         self.signal_jet_energy_res_uncertainty      = 0.000 ;
        elif self.channel == "em":
         self.signal_lepton_energy_scale_uncertainty = 0.000 ;
         self.signal_lepton_energy_res_uncertainty   = 0.000 ;
         self.signal_jet_energy_res_uncertainty      = 0.000 ;

        self.eff_vtag_model = 0. ;
        self.eff_vtag_mass_sf = 0.;

        self.uncertainties(self.signal_sample)
#	self.xs_rescale = 0.01 
	self.xs_rescale = 1.
        if (self.interpolate):
            self.xs_rescale = 1.
	                            
        #### sigma and mean signal systematic inflation
        self.mean_signal_uncertainty_jet_scale  = 0.000 ;
        self.mean_signal_uncertainty_lep_scale  = 0.000 ;
        self.sigma_signal_uncertainty_jet_scale = 0.000 ;
        self.sigma_signal_uncertainty_jet_res   = 0.000 ;
        self.sigma_signal_uncertainty_lep_scale = 0.000 ;

        #### Set systematic on the Wjets shape   and TTbar due to PS, fitting function etc..
        self.shape_para_error_WJets0 = 0.0;
        self.shape_para_error_alpha  = 0.0;
        self.shape_para_error_TTbar = 0.0;
        self.shape_para_error_VV    = 0.;
        self.shape_para_error_STop  = 0.;

    
                                                                
        # shape parameter uncertainty
        self.FloatingParams=RooArgList("floatpara_list");

    def uncertainties(self,label_tstring):
    
       self.signal_jet_energy_scale_uncertainty_up = -0.000
       self.signal_jet_energy_scale_uncertainty_down = 0.000
       infile = open('sys/JESsys_%s_%s.txt' %(self.channel,self.wtagger_label),'r')
       for l in infile:
        if l.find('mass') != -1: continue
        if l.split(' ')[0].find(label_tstring) != -1:
         self.signal_jet_energy_scale_uncertainty_up = float(l.split(' ')[1])/100. 
	 self.signal_jet_energy_scale_uncertainty_down = float(l.split(' ')[2])/100.
       infile.close()
       
       self.signal_jet_mass_scale_uncertainty_up = -0.000
       self.signal_jet_mass_scale_uncertainty_down = 0.00
       infile = open('sys/JMSsys_%s_%s.txt' %(self.channel,self.wtagger_label),'r')
       for l in infile:
        if l.find('mass') != -1: continue
        if l.split(' ')[0].find(label_tstring) != -1:
         self.signal_jet_mass_scale_uncertainty_up = float(l.split(' ')[1])/100. 
	 self.signal_jet_mass_scale_uncertainty_down = float(l.split(' ')[2])/100.
       infile.close()

       self.signal_jet_mass_res_uncertainty_up = 0.000
       self.signal_jet_mass_res_uncertainty_down = -0.000
       infile = open('sys/JMRsys_%s_%s.txt' %(self.channel,self.wtagger_label),'r')
       for l in infile:
        if l.find('mass') != -1: continue
        if l.split(' ')[0].find(label_tstring) != -1:
         self.signal_jet_mass_res_uncertainty_up = float(l.split(' ')[1])/100. 
	 self.signal_jet_mass_res_uncertainty_down = float(l.split(' ')[2])/100.
       infile.close()

       self.signal_xsec_NNLO_uncertainty = 0.00 # scale = 11%, PDF 13% at 1000 GeV
       infile = open('sys/signal_xsec_uncertainty.txt','r')
       for l in infile:
        if l.find('mass') != -1: continue
        if l.split(' ')[0].find(label_tstring) != -1:
         self.signal_xsec_NNLO_uncertainty = float(l.split(' ')[1])
       infile.close()
       
    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0.4, y_offset_low=0.1, x_offset_high =0.27, y_offset_high =0.1, TwoCoulum =1., isalpha=False, ismj=False):
        print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.37+x_offset_low, 0.50+y_offset_low, 0.72+x_offset_high, 0.82+y_offset_high, "", "NDC");            
            theLeg.SetName("theLegend");
            if ismj: theLeg = TLegend(0.3715365,0.505,0.8526448,0.845, "", "NDC"); 
            if TwoCoulum :
                theLeg.SetNColumns(2);
            if isalpha: theLeg = TLegend(0.3944724,0.4370629,0.7650754,0.8374126, "", "NDC");  
            
        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.05);
        theLeg.SetTextFont(42);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        

        #if   self.categoryID == 0 or self.categoryID == 2: legHeader="e#nu";
        #elif self.categoryID == 1 or self.categoryID == 3: legHeader="#mu#nu";

        if   self.channel == 'el': legHeader="W#rightarrowe#nu";
        elif self.channel == 'mu': legHeader="W#rightarrow#mu#nu";
        elif self.channel == 'em': legHeader="W#rightarrow#ell#nu";
        
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
	  #if objName.find("TPave") != -1: continue
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "Observed","PE");  objName_before=objName;                 
            #elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "Data "+legHeader,"PE");  objName_before=objName;                 
            else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        
                   
        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "V+jets","F");  objName_before=objName;                 
            else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
	    if objName.find("TPave") != -1: continue
            if objName == "errorband" : objName = "Bkg. Uncertainty";
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Stat. Uncertainty","F");
                else:
                    #if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "SingleTop","F");
                    if TString(objName).Data()=="STop" : print "Do not add Single Top to legend"; # Commendted CWR comments
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "Top quark","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "QCD+EW ZV","F");
                    elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
                    elif TString(objName).Data()=="WJets" : objName_before=objName; entryCnt = entryCnt+1; continue;
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Contains("Wprime") and TString(objName).Contains("M2000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "W' M_{W'}=2 TeV";
                    elif TString(objName).Contains("RS1") or TString(objName).Contains("Bulk"):
                       print "!!!! objName ",  objName
                       prefix = ""
                       if TString(objName).Contains("RS1"): prefix = "RS"
                       elif TString(objName).Contains("Bulk"): prefix = "Bulk"
                       elif TString(objName).Contains("Higgs"): prefix = "Higgs"
                       if TString(objName).Contains("_Higgs650"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = prefix+" M_{H}=0.65 TeV (#times100)";
                       if TString(objName).Contains("_Higgs750"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = prefix+" M_{H}=0.75 TeV (#times100)";
                       if TString(objName).Contains("_Higgs1000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = prefix+" M_{H}=1.0 TeV (#times100)";
                       if TString(objName).Contains("_BulkGraviton600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=0.6 TeV (#times100)";
                       if TString(objName).Contains("_BulkGraviton700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=0.7 TeV (#times100)";
                       if TString(objName).Contains("BulkGraviton750"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=750 GeV (#times20)";
                       if TString(objName).Contains("_BulkGraviton800"):
                           objName_signal_graviton = theObj ; 
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=0.8 TeV (#times100)";
                       if TString(objName).Contains("_BulkGraviton900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=0.9 TeV (#times100)";
                       if TString(objName).Contains("_BulkGraviton1000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.1 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.2 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.3 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.4 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.5 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.6 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.7 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.8 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M1900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=1.9 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2.1 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2.2 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2.3 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2.4 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M2500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=2.5 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M3000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=3 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_M3500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=3.5 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_4000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=4 TeV (#times100)";
                       if TString(objName).Contains("_WW_lvjj_4500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "G_{"+prefix+"}"+" M_{G}=4.5 TeV (#times100)";
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;

    def get_canvas(self,cname,isalpha=False):

       #tdrstyle.setTDRStyle()
       CMS_lumi.lumi_13TeV = "35.9 fb^{-1}"
       CMS_lumi.writeExtraText = 1
       CMS_lumi.extraText = "Preliminary"

       iPos = 11
       if( iPos==0 ): CMS_lumi.relPosX = 0.15

       H_ref = 600; 
       W_ref = 800; 
       W = W_ref
       H  = H_ref

       T = 0.08*H_ref
       B = 0.12*H_ref 
       L = 0.12*W_ref
       R = 0.06*W_ref

       canvas = ROOT.TCanvas(cname,cname,W,H)
       canvas.SetFillColor(0)
       canvas.SetBorderMode(0)
       canvas.SetFrameFillStyle(0)
       canvas.SetFrameBorderMode(0)
       canvas.SetLeftMargin( L/W )
       canvas.SetRightMargin( R/W )
       canvas.SetTopMargin( T/H )
       canvas.SetBottomMargin( B/H+0.03 )
       canvas.SetTickx()
       canvas.SetTicky()
       if isalpha:
        canvas.SetTicky(0)
       
       return canvas

    #### just drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0, isalpha=0):

        print "############### draw the canvas without pull ########################"
        cMassFit = self.get_canvas("cMassFit",isalpha)#TCanvas("cMassFit","cMassFit", 600,600);

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/200)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())
        #elif isalpha:    
        #    in_obj.GetYaxis().SetRangeUser(0.00001,0.28)

        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.04);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        self.leg.SetTextSize(0.031); 
        if isalpha: self.leg.SetTextSize(0.038)

        #banner = self.banner4Plot();
        #banner.Draw();

	cMassFit.Update()
	cMassFit.cd()
	CMS_lumi.CMS_lumi(cMassFit, 4, 11)	
	cMassFit.cd()
	cMassFit.Update()
	cMassFit.RedrawAxis()
	frame = cMassFit.GetFrame()
	frame.Draw()   
	cMassFit.cd()
	cMassFit.Update()
	        
        Directory=TString(in_directory); #+"_%02d_%02d/"%(options.cprime,options.BRnew));
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
        else:
            rlt_file.ReplaceAll(".root","");
            rlt_file = rlt_file.Append(".png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        if logy:
            #in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
            in_obj.GetYaxis().SetRangeUser(1e-3,1e3); 	# CWR Comment
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):
      print "############### draw the banner ########################"

      if iswithpull:
       if self.channel=="el":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       elif self.channel=="mu":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       elif self.channel=="em":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       banner.SetNDC(); banner.SetTextSize(0.041);
      else:
       if self.channel=="el":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       if self.channel=="mu":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       if self.channel=="em":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
        banner = TLatex(0.18,0.96,"CM                                          L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
       banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
      return banner;
    
    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, rrv_x, datahist, mplot, mplot_pull,ndof,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0,ismj=0,isPull=0):# mplot + pull

        print "############### draw the canvas with pull ########################" 
	#hist_ = datahist.createHistogram(rrv_x.GetName(),int(rrv_x.getBins()/self.narrow_factor))
        #chi2_ = self.calculate_chi2(datahist,rrv_x,mplot,ndof,ismj)
	mplot.GetXaxis().SetTitle("")
	#mplot.GetXaxis().SetMaxDigits(3)
	#mplot.GetXaxis().SetTitleOffset(1.1);
        #mplot.GetYaxis().SetTitleOffset(1.3);
        #mplot.GetXaxis().SetTitleSize(0.055);
        #mplot.GetYaxis().SetTitleSize(0.055);
        #mplot.GetXaxis().SetLabelSize(0.045);
        #mplot.GetYaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
	mplot.GetXaxis().SetLabelSize(0);
        #mplot_pull.GetXaxis().SetLabelSize(0.14);
        #mplot_pull.GetYaxis().SetLabelSize(0.14);
        #mplot_pull.GetYaxis().SetTitleSize(0.15);
        #mplot_pull.GetYaxis().SetNdivisions(205);
	
        cMassFit = self.get_canvas("cMassFit")#TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,1,0.30); #pad1 - pull
         pad2=TPad("pad2","pad2",0.,0.3,1.,1. ); #pad0

   	 pad2.SetRightMargin(0.1);
   	 pad2.SetTopMargin(0.1);
   	 pad2.SetBottomMargin(0.0001);
   	 pad1.SetRightMargin(0.1)
   	 pad1.SetTopMargin(0)
   	 pad1.SetBottomMargin(0.4)   
	 #pad1.SetRightMargin(0.05)
	 #pad2.SetRightMargin(0.05)
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
	
	'''	
	if ismj:
	 pt = ROOT.TPaveText(0.6243719,0.4080919,0.8756281,0.547952,"NDC")
	 pt.SetTextSize(0.03746254)
	else:
	 pt = ROOT.TPaveText(0.5175879,0.7152847,0.8027638,0.8551449,"NDC")
	 pt.SetTextSize(0.054)
	 
	pt.SetTextFont(62)	
	pt.SetTextAlign(12)
	pt.SetFillColor(0)
	pt.SetBorderSize(0)
	pt.SetFillStyle(0)
	text = pt.AddText("#chi^2/d.o.f = %.2f/%i = %.2f" %(chi2_[0],chi2_[1],chi2_[0]/chi2_[1]))
	text.SetTextFont(62)
	'''
	
        mplot.Draw();
	#pt.Draw()

        #banner = self.banner4Plot(1);
        #banner.Draw();

        pad1.cd();
        mplot_pull.Draw("AP");

        mplot.Print("v");
        if mplot.FindObject("errorband") != 0: 
          print " Pull  !!!! "       
          mplot_pull.Print("v");
#        TGraphAsymmErrors::errorband
#        npoints = errorband.GetN();
#        errorband_pull = errorband.Clone();
#        for ipoint in range(1,errorband_pull->GetN()):

        medianLine = TLine(mplot.GetXaxis().GetXmin(),0.,mplot.GetXaxis().GetXmax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
	medianLine.Draw()
	mplot_pull.Draw("Psame");
	
        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

	cMassFit.Update()
	pad2.cd()
	CMS_lumi.CMS_lumi(pad2, 4, 11)	
	pad2.cd()
	pad2.Update()
	pad2.RedrawAxis()
	frame = pad2.GetFrame()
	frame.Draw()   
	cMassFit.cd()
	cMassFit.Update()
			
        ## create the directory where store the plots
        Directory = TString(in_directory);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());
        
        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".root",".C");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".C","_log.png");
	pad2.SetLogy()
        mplot.GetYaxis().SetRangeUser(0.06,mplot.GetMaximum()*200);
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(in_file_name);
        if string_file_name.EndsWith(".root"):
            string_file_name.ReplaceAll(".root","_"+in_model_name);
        else:
            string_file_name.ReplaceAll(".root","");
            string_file_name.Append("_"+in_model_name);

        #if logy:
        #    mplot.GetYaxis().SetRangeUser(0.002,mplot.GetMaximum()*200);
	#    mplot.GetXaxis().SetTitle("m_{ZV} (GeV)")
        #    pad2.SetLogy() ;
        #    pad2.Update();
        #    cMassFit.Update();
        #    rlt_file.ReplaceAll(".root","_log.root");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".root",".pdf");
        #    cMassFit.SaveAs(rlt_file.Data());
        #    rlt_file.ReplaceAll(".pdf",".png");
        #    cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    #### draw canvas with plots with pull
    def draw_canvas_with_pull2(self, rrv_x, datahist, mplot, mplotP, mplot_pull,ndof,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0,ismj=0,isPull=0):# mplot + pull

        print "############### draw the canvas with pull ########################" 
	#hist_ = datahist.createHistogram(rrv_x.GetName(),int(rrv_x.getBins()/self.narrow_factor))
        chi2_ = self.calculate_chi2(datahist,rrv_x,mplot,ndof,ismj)
	mplot.GetXaxis().SetTitle("")
	#mplot.GetXaxis().SetTitleOffset(1.1);
        #mplot.GetYaxis().SetTitleOffset(1.3);
        #mplot.GetXaxis().SetTitleSize(0.055);
        #mplot.GetYaxis().SetTitleSize(0.055);
        #mplot.GetXaxis().SetLabelSize(0.045);
        #mplot.GetYaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetTitleSize(0.07)
        mplot.GetYaxis().SetTitleOffset(0.9)
        mplot.GetYaxis().SetLabelSize(0.06)
	mplot.GetXaxis().SetLabelSize(0);
        #mplot_pull.GetXaxis().SetLabelSize(0.14);
        #mplot_pull.GetYaxis().SetLabelSize(0.14);
        #mplot_pull.GetYaxis().SetTitleSize(0.15);
        #mplot_pull.GetYaxis().SetNdivisions(205);
	
        cMassFit = self.get_canvas("cMassFit")#TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,1,0.30); #pad1 - pull
         pad2=TPad("pad2","pad2",0.,0.3,1.,1. ); #pad0

   	 pad2.SetRightMargin(0.1);
   	 pad2.SetTopMargin(0.1);
   	 pad2.SetBottomMargin(0.0001);
   	 pad1.SetRightMargin(0.1)
   	 pad1.SetTopMargin(0)
   	 pad1.SetBottomMargin(0.4)   
	 #pad1.SetRightMargin(0.05)
	 #pad2.SetRightMargin(0.05)
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
	
		
	if ismj:
	 pt = ROOT.TPaveText(0.6243719,0.4080919,0.8756281,0.547952,"NDC")
	 pt.SetTextSize(0.03746254)
	else:
	 pt = ROOT.TPaveText(0.5175879,0.7152847,0.8027638,0.8551449,"NDC")
	 pt.SetTextSize(0.054)
	 
	pt.SetTextFont(62)	
	pt.SetTextAlign(12)
	pt.SetFillColor(0)
	pt.SetBorderSize(0)
	pt.SetFillStyle(0)
	text = pt.AddText("#chi^2/d.o.f = %.2f/%i = %.2f" %(chi2_[0],chi2_[1],chi2_[0]/chi2_[1]))
	text.SetTextFont(62)
	
	
        mplot.Draw();
	#pt.Draw()

        #banner = self.banner4Plot(1);
        #banner.Draw();

        pad1.cd();
        mplot_pull.Draw("AP");
        mplotP.Draw("same");
        medianLine = TLine(mplot.GetXaxis().GetXmin(),0.,mplot.GetXaxis().GetXmax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
	medianLine.Draw()
	mplot_pull.Draw("Psame");
	
        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

	cMassFit.Update()
	pad2.cd()
	CMS_lumi.CMS_lumi(pad2, 4, 11)	
	pad2.cd()
	pad2.Update()
	pad2.RedrawAxis()
	frame = pad2.GetFrame()
	frame.Draw()   
	cMassFit.cd()
	cMassFit.Update()
			
        ## create the directory where store the plots
        Directory = TString(in_directory);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
        if rlt_file.EndsWith(".root"):
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
        else:
            TString(in_model_name).ReplaceAll(".root","");
            rlt_file.ReplaceAll(".root","");
            rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());
        
        rlt_file.ReplaceAll(".pdf",".root");
        cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(in_file_name);
        if string_file_name.EndsWith(".root"):
            string_file_name.ReplaceAll(".root","_"+in_model_name);
        else:
            string_file_name.ReplaceAll(".root","");
            string_file_name.Append("_"+in_model_name);

        if logy:
            mplot.GetYaxis().SetRangeUser(0.002,mplot.GetMaximum()*200);
	    mplot.GetXaxis().SetTitle("m_{ZV} (GeV)")
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    def calculate_chi2(self,hist,rrv_x,mplot_orig,ndof,ismj):

        pulls = array('d',[])
        print "############### calculate chi2 (new) ########################"
        hpull = mplot_orig.pullHist();
	bins = 0
	bins_ = 0
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
	  hist.get(bins_)
	  hist.weightError(RooAbsData.SumW2)
	  print x,y,bins_,hist.get(bins_).getRealValue(rrv_x.GetName()),hist.weight(),hist.weightError(RooAbsData.SumW2)
	  if hist.weight() != 0: pulls.append(y)
	  #print x,y,hist.GetBinCenter(bins_),hist.GetBinContent(bins_)
          #if not(ismj) and y != 0 and TMath.Abs(y) < 4: pulls.append(y)
	  #elif ismj and y != 0 and (x < 65 or x > 135):
	  # pulls.append(y) 
          # bins+=1
	  #else: print "Bin %f is empty!" %x 
	  bins_+=1
	  
	chi2 = 0
	for p in pulls:
	 chi2+=(p*p)
	 
	#if ismj:
	# ndof = ndof - (hpull.GetN()-bins) 
	# print hpull.GetN(),bins,ndof
	print "Chi2/ndof = %f/%f = %f" %(chi2,ndof,chi2/ndof)
	return chi2,ndof
	       
    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
	  #print x,y
          if(y == 0):
           hpull.SetPoint(ipoint,x,10)
       
   	gt = ROOT.TH1F("gt","gt",int(rrv_x.getBins()/self.narrow_factor),rrv_x.getMin(),rrv_x.getMax());
   	gt.SetMinimum(-3.999);
   	gt.SetMaximum(3.999);
   	gt.SetDirectory(0);
   	gt.SetStats(0);
   	gt.SetLineStyle(0);
   	gt.SetMarkerStyle(20);
   	gt.GetXaxis().SetTitle(rrv_x.GetTitle());
   	gt.GetXaxis().SetLabelFont(42);
   	gt.GetXaxis().SetLabelOffset(0.02);
   	gt.GetXaxis().SetLabelSize(0.15);
   	gt.GetXaxis().SetTitleSize(0.15);
   	gt.GetXaxis().SetTitleOffset(1.2);
   	gt.GetXaxis().SetTitleFont(42);
   	gt.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
   	gt.GetYaxis().CenterTitle(True);
   	gt.GetYaxis().SetNdivisions(205);
   	gt.GetYaxis().SetLabelFont(42);
   	gt.GetYaxis().SetLabelOffset(0.007);
   	gt.GetYaxis().SetLabelSize(0.15);
   	gt.GetYaxis().SetTitleSize(0.15);
   	gt.GetYaxis().SetTitleOffset(0.35);
   	gt.GetYaxis().SetTitleFont(42);
   	#gt.GetXaxis().SetNdivisions(505)
	hpull.SetHistogram(gt)
	return hpull
        #mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        #medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        #mplot_pull.addObject(medianLine);
        #mplot_pull.addPlotable(hpull,"P");
        #mplot_pull.SetTitle("");
        #mplot_pull.GetXaxis().SetTitle("");
        #mplot_pull.GetYaxis().SetRangeUser(-5,5);
        #mplot_pull.GetYaxis().SetTitleSize(0.10);
        #mplot_pull.GetYaxis().SetLabelSize(0.10);
        #mplot_pull.GetXaxis().SetTitleSize(0.10);
        #mplot_pull.GetXaxis().SetLabelSize(0.10);
        #mplot_pull.GetYaxis().SetTitleOffset(0.40);
        #mplot_pull.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
        #mplot_pull.GetYaxis().CenterTitle();

        #return mplot_pull;

    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1);
        
        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
          N = data_plot.GetY()[iPoint];
          if N==0 : L = 0;
          else : L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
          U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
          data_plot.SetPointEYlow(iPoint, N-L);
          data_plot.SetPointEYhigh(iPoint,U-N);
          data_plot.SetPointEXlow(iPoint,0);        
          data_plot.SetPointEXhigh(iPoint,0);        
        
        mplot.addPlotable(data_plot,"PE");

    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region):
        return self.workspace4fit_.pdf("model"+label+mlvj_region+"_"+self.channel+"_mlvj");

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_signal_region"):
        print "########### Fixing a general mlvj model  ############"
        rdataset_General_mlvj = self.workspace4fit_.data("rdataset%s%s_%s_mlvj"%(label, mlvj_region,self.channel))
        model_General = self.get_mlvj_Model(label,mlvj_region);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_mlvj_Model(label,mlvj_region);

    ###### get TTbar model mlvj in a region 
    def get_TTbar_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing TTbar mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_TTbar_xww",mlvj_region);

    ###### get Single Top model mlvj in a region 
    def get_STop_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_STop_xww",mlvj_region);

    ###### get VV mlvj in a region 
    def get_VV_mlvj_Model(self, mlvj_region="_signal_region"):
        print "########### Fixing VV mlvj for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_VV_xww",mlvj_region);
	
    ### get an mj model from the workspace givin the label
    def get_mj_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj")

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label ):
        print "########### Fixing a general mj model  ############"
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_General = self.get_mj_Model(label);
        rdataset_General_mj.Print();
        model_General.Print();
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                param.Print();
	    param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_mj_Model(self,label="_TTbar_xww"):
        print "########### Fixing only the TTbar mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the stop component using the default label --> for extended pdf
    def get_STop_mj_Model(self,label="_STop_xww"):
        print "########### Fixing only the Stop mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_mj_Model(self,label="_VV_xww"):
        print "########### Fixing only the VV mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_mj_Model(self,label):
        print "########### Fixing only the WJets mj Shape --> just the printed parameters  ############"
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.channel))
        model_WJets = self.get_mj_Model(label);
        rdataset_WJets_mj.Print();
        model_WJets.Print();
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj);
        par=parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            if ( paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")):
             param.setConstant(kTRUE);
             param.Print();
            else:
             param.setConstant(0);
            param=par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.channel))

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj  ############"        
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        model = self.get_mlvj_Model(label,mlvj_region);
        model.Print();
        rdataset.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param=par.Next()

    ### fix a pdf in a different way --> for RooAbsPdf 
    def fix_Pdf(self,model_pdf,argset_notparameter):
        print "########### Fixing a RooAbsPdf for mlvj or mj  ############"        
        parameters = model_pdf.getParameters(argset_notparameter);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()


#############

    #### Method to make a RooAbsPdf setting parameters from external interpolation (available only for DoubleCB)    
    def make_Pdf_from_interpolation(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):
        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4fit_.var("rrv_mass_lvj");

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "DoubleCB_v1":
            label_tstring=TString(label);

            fileIn_name = TString("interpolationFiles/"+options.sample+"_"+self.channel+"_"+options.category+".root")#str(int(options.mass)));
            print "open file for interpolation of DoubleCB shape: ",fileIn_name

            fileIn = TFile(fileIn_name.Data());
            g_mean = fileIn.Get("mean");
            g_sigma = fileIn.Get("sigma");
            g_n1 = fileIn.Get("n1");
            g_n2 = fileIn.Get("n2");
            g_alpha1 = fileIn.Get("alpha1");
            g_alpha2 = fileIn.Get("alpha2");

            print "########### Double CB for Bulk graviton mlvj ############"

#            if label_tstring.Contains("RS1G_WW" or label_tstring.Contains("BulkGraviton") or label_tstring.Contains("Wprime_WZ")):
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, g_mean.Eval(int(options.mass)));
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, g_sigma.Eval(int(options.mass)));
            rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, g_n1.Eval(int(options.mass)));
            rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,g_alpha2.Eval(int(options.mass)));
            rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,g_n2.Eval(int(options.mass)));
            rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,g_alpha1.Eval(int(options.mass)));

            #set parameters of the CB to the interpolated value
#            rrv_mean_CB.setVal(g_mean.Eval(int(options.mass))); 
#            rrv_sigma_CB.setVal(g_sigma.Eval(int(options.mass))); 
#            rrv_n1_CB.setVal(g_n1.Eval(int(options.mass))); 
#            rrv_n2_CB.setVal(g_n2.Eval(int(options.mass))); 
#            rrv_alpha1_CB.setVal(g_alpha1.Eval(int(options.mass))); 
#            rrv_alpha2_CB.setVal(g_alpha2.Eval(int(options.mass))); 

            rrv_mean_CB.Print()
            rrv_sigma_CB.Print()
            rrv_n1_CB.Print()
            rrv_n2_CB.Print()
            rrv_alpha1_CB.Print()
            rrv_alpha2_CB.Print()
                         
            rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes_13TeV","CMS_sig_p1_jes_13TeV",0);
            rrv_mean_scale_p1.setConstant(kTRUE);
            if self.channel == "mu" :             
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_m_13TeV","CMS_sig_p1_scale_m_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.channel == "el" :
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_e_13TeV","CMS_sig_p1_scale_e_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.channel == "em":
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em_13TeV","CMS_sig_p1_scale_em_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
                
            rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_lep_scale));
            rrv_mean_scale_X1.setConstant(kTRUE);
            rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_jet_scale));
            rrv_mean_scale_X2.setConstant(kTRUE);

            rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label+"_"+self.channel,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));
            
            if self.channel == "mu":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_m_13TeV","CMS_sig_p2_scale_m_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.channel == "el":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_e_13TeV","CMS_sig_p2_scale_e_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.channel == "em":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em_13TeV","CMS_sig_p2_scale_em_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
                         
            rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer_13TeV","CMS_sig_p2_jer_13TeV",0);
            rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes_13TeV","CMS_sig_p2_jes_13TeV",0);
            rrv_sigma_scale_p2.setConstant(kTRUE);
            rrv_sigma_scale_p3.setConstant(kTRUE);

            rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_lep_scale));
            rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_scale));
            rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_res));
            rrv_mean_sigma_X1.setConstant(kTRUE);
            rrv_mean_sigma_X2.setConstant(kTRUE);
            rrv_mean_sigma_X3.setConstant(kTRUE);

            rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label+"_"+self.channel,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

            model_pdf = ROOT.RooDoubleCrystalBall("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB);

            
        ## return the pdf
        getattr(self.workspace4fit_,"import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)

###########

    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):
        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4fit_.var("rrv_mass_lvj");

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "DoubleCB_v1":
            label_tstring=TString(label);
            print "########### Double CB for Bulk graviton mlvj ############"
            print " %s %s " % (label,self.channel)
            print "########### Double CB for Bulk graviton mlvj ############"

            if label_tstring.Contains("600") and (label_tstring.Contains("RS1G_WW") or label_tstring.Contains("BulkGraviton") or label_tstring.Contains("Higgs") or label_tstring.Contains("Wprime_WZ")):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 600, 550, 650);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);
                    rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3.,0.5,6.);
                    rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,3,0.5,6.);

            elif label_tstring.Contains("4500") and (label_tstring.Contains("RS1G_WW") or label_tstring.Contains("BulkGraviton") or label_tstring.Contains("Higgs") or label_tstring.Contains("Wprime_WZ")):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,4520,4500,4540);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,225,160,260);
                    rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 10.,0.01,35);
                    rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);
                    rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,0.01,35);
                    rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1.5,0.5,3.);

            rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes_13TeV","CMS_sig_p1_jes_13TeV",0);
            rrv_mean_scale_p1.setConstant(kTRUE);
            if self.channel == "mu" :             
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_m_13TeV","CMS_sig_p1_scale_m_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.channel == "el" :
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_e_13TeV","CMS_sig_p1_scale_e_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.channel == "em":
             rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em_13TeV","CMS_sig_p1_scale_em_13TeV",0);
             rrv_mean_scale_p2.setConstant(kTRUE);
                
            rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_shift_scale_lep"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_lep_scale));
            rrv_mean_scale_X1.setConstant(kTRUE);
            rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_shift_scale_jes"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_jet_scale));
            rrv_mean_scale_X2.setConstant(kTRUE);

            rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label+"_"+self.channel,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));
            
            if self.channel == "mu":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_m_13TeV","CMS_sig_p2_scale_m_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.channel == "el":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_e_13TeV","CMS_sig_p2_scale_e_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.channel == "em":
             rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em_13TeV","CMS_sig_p2_scale_em_13TeV",0);
             rrv_sigma_scale_p1.setConstant(kTRUE);
                         
            rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer_13TeV","CMS_sig_p2_jer_13TeV",0);
            rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes_13TeV","CMS_sig_p2_jes_13TeV",0);
            rrv_sigma_scale_p2.setConstant(kTRUE);
            rrv_sigma_scale_p3.setConstant(kTRUE);

            rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_lep_scale));
            rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_scale));
            rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_shift_res"+label+"_"+self.channel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_res));
            rrv_mean_sigma_X1.setConstant(kTRUE);
            rrv_mean_sigma_X2.setConstant(kTRUE);
            rrv_mean_sigma_X3.setConstant(kTRUE);

            rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label+"_"+self.channel,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

            model_pdf = ROOT.RooDoubleCrystalBall("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB);

        ## ExpSlowFastFall pdf for W+jets bkg fit
        if in_model_name == "ExpSlowFastFall":

            print "########### ExpSlowFastFall  funtion for W+jets mlvj ############"
            rrv_c_ExpSlowFastFall = RooRealVar("rrv_c_ExpSlowFastFall"+label+"_"+self.channel,"rrv_c_ExpSlowFastFall"+label+"_"+self.channel,-2.28320e-03,-10.28320e-03,-0.28320e-05);
            rrv_n_ExpSlowFastFall = RooRealVar("rrv_n_ExpSlowFastFall"+label+"_"+self.channel,"rrv_n_ExpSlowFastFall"+label+"_"+self.channel,-2.28320e-03,-10.28320e-03,-0.28320e-05);
            #model_pdf = ROOT.RooExpSlowFastFall("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpSlowFastFall, rrv_n_ExpSlowFastFall);
            model_pdf_1 = ROOT.RooExpDan("model_pdf_1"+label+"_"+self.channel+mass_spectrum,"model_pdf_1"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpSlowFastFall, rrv_n_ExpSlowFastFall);

            rrv_c_Gaus = RooRealVar("rrv_c_Gaus"+label+"_"+self.channel,"rrv_c_Gaus"+label+"_"+self.channel,583,200, 700);
            rrv_n_Gaus = RooRealVar("rrv_n_Gaus"+label+"_"+self.channel,"rrv_n_Gaus"+label+"_"+self.channel,100,5,400);
            model_pdf_2 = RooGaussian("model_pdf_2"+label+"_"+self.channel+mass_spectrum,"model_pdf_2"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Gaus, rrv_n_Gaus);

            rrv_frac_Gaus = RooRealVar("rrv_frac_Gaus"+label+"_"+self.channel,"rrv_frac_Gaus"+label+"_"+self.channel,0.5,0.0,1.0);

            model_pdf = ROOT.RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(rrv_frac_Gaus));
            #model_pdf = ROOT.RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(frac),1);
                                                                                                             
        ## Landau pdf for W+jets bkg fit
        if in_model_name == "Landau":
            print "########### Landau  funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if label_tstring.Contains("signal_region"):
	    	print "==> initialize Signal region pars"
            	rrv_MPV_Landau = RooRealVar("rrv_MPV_Landau"+label+"_"+self.channel,"rrv_MPV_Landau"+label+"_"+self.channel,483,200,800);
            	rrv_Sigma_Landau = RooRealVar("rrv_Sigma_Landau"+label+"_"+self.channel,"rrv_Sigma_Landau"+label+"_"+self.channel,20.0,1.0,100.0);
            	rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-4.8438e-03,-59.5655e-01,-5.5655e-04);
            	rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, -985.40, -1285.4, -185.4);
            	rrv_frac_ExpN = RooRealVar("rrv_frac_ExpN"+label+"_"+self.channel,"rrv_frac_ExpN"+label+"_"+self.channel,0.5,0.0,0.8);
            elif label_tstring.Contains("sb_lo"):
	    	print "==> initialize side-band region pars"
            	rrv_MPV_Landau = RooRealVar("rrv_MPV_Landau"+label+"_"+self.channel,"rrv_MPV_Landau"+label+"_"+self.channel,500,100,1000);
            	rrv_Sigma_Landau = RooRealVar("rrv_Sigma_Landau"+label+"_"+self.channel,"rrv_Sigma_Landau"+label+"_"+self.channel,20.0,1.0,100.0);
            	rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-5.32e-3,-1e-1,-2e-6);
            	rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, -64.3, -1e4, 1e5);
            	rrv_frac_ExpN = RooRealVar("rrv_frac_ExpN"+label+"_"+self.channel,"rrv_frac_ExpN"+label+"_"+self.channel,0.5,0.0,0.8);
            
	    model_pdf = ROOT.RooLandauPlusGaus("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_frac_ExpN,rrv_MPV_Landau,rrv_Sigma_Landau, rrv_c_ExpN, rrv_n_ExpN)
                                                                                                             
	## Function for mWW fit from 160 to tail
        if in_model_name == "ErfPow2":
            rrv_c0     = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.channel,"rrv_c0_ErfPow2"+label+"_"+self.channel,20.696,-500.,500) 
            rrv_c1     = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.channel,"rrv_c1_ErfPow2"+label+"_"+self.channel,3.092,-100,100)
            rrv_c2     = RooRealVar("rrv_c2_ErfPow2"+label+"_"+self.channel,"rrv_c2_ErfPow2"+label+"_"+self.channel,3.092,-100,100)
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.channel,"rrv_offset_ErfPow2"+label+"_"+self.channel,-482,-900,100)
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.channel,"rrv_width_ErfPow2"+label+"_"+self.channel,-4665,-7000,100.0)
            
	    model_pdf = ROOT.RooErfPow3Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_c2,rrv_offset,rrv_width)

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":
            label_tstring=TString(label);

            print "########### ExpN funtion for W+jets mlvj ############"
            if label_tstring.Contains("signal_region"):
               rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-0.00347556,-1e-1,-1e-6);
               rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, -173.987, -1e4, 1e5);
            elif label_tstring.Contains("sb_lo"):
               rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.channel,"rrv_c_ExpN"+label+"_"+self.channel,-0.00359157,-1e-1,-1e-6);
               rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.channel,"rrv_n_ExpN"+label+"_"+self.channel, -4.47341, -1e4, 1e5);

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);
                                                                                                             
        ## levelled exp for W+jets bkg fit
        if in_model_name == "LandauExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.channel == "mu" or self.channel == "em":
             if ismc == 1 and label_tstring.Contains("sb_lo"):
               #rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 99,10,255);
               rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
               rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-1e-2,7.5e-2);                   
               #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-2e-2,7.5e-2);                   
             elif ismc == 1 and label_tstring.Contains("signal_region"):
              # rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 110,20,242);
               rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 105,20,500);
               #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1e-2,7.5e-2);
               rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1,7.5e-2);
             elif ismc == 0 :  
                 rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,40,280);
                 rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);    
      
      	    print "HERE I AM"
	    rrv_s_ExpTail.Print()     
	    rrv_a_ExpTail.Print()     
            model_pdf_1     = ROOT.RooExpTailPdf("model_pdf_1"+label+"_"+self.channel+mass_spectrum,"model_pdf_1"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

            rrv_MPV_Landau = RooRealVar("rrv_MPV_Landau"+label+"_"+self.channel,"rrv_MPV_Landau"+label+"_"+self.channel,500,100,2500);
            rrv_Sigma_Landau = RooRealVar("rrv_Sigma_Landau"+label+"_"+self.channel,"rrv_Sigma_Landau"+label+"_"+self.channel,150,100,500);
            model_pdf_2 = ROOT.RooLandau("model_pdf_2"+label+"_"+self.channel+mass_spectrum,"model_pdf_2"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_MPV_Landau, rrv_Sigma_Landau);

	    frac = RooRealVar("frac","frac",0.8,0.,1.);
            model_pdf = ROOT.RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(frac),1);

        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.wtagger_label.find("LP") != -1:
             rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
             rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e-2,1e6);
            else:
                if self.channel == "el" :
                 if ismc == 1 and label_tstring.Contains("sb_lo"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 139,0.,355);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2e-2,-1.e-2,5.5e-2);                     
                 elif ismc == 1 and label_tstring.Contains("signal_region"):
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 162,18,395);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1.6e-2,-1.e-2,5.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,70,240);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);
                           
                if self.channel == "mu" or self.channel == "em":
                 if ismc == 1 and label_tstring.Contains("sb_lo"):
                   #rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 99,10,255);
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-1e-2,7.5e-2);                   
                   #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 3e-2,-2e-2,7.5e-2);                   
                 elif ismc == 1 and label_tstring.Contains("signal_region"):
                  # rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 105,20,242);
                   rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 105,20,500);
                   #rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1e-2,7.5e-2);
                   rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 2.9e-2,-1,7.5e-2);
                 elif ismc == 0 :  
                     rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 161,40,280);
                     rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 8e-3,-1e-2,1.3e-1);    
      
      	    print "HERE I AM"
	    rrv_s_ExpTail.Print()     
	    rrv_a_ExpTail.Print()     
            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

        ## sum of two exponential 
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,-0.00358516,-10,10.);
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);

        ## Erf times for mj spectrum
        if in_model_name == "ErfExp" :
            print "########### Erf*Exp for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.0323819,-0.1,-1e-4);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,65.0,0.,200);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,34.71,0., 200.);
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v1" :
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.006,-0.1,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400.,550.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,70.,10,100.);
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v2" : 
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400.,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel,"rrv_residue_ErfExp"+label+"_"+self.channel,0.,0.,1.);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2. "%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v3" : #different init-value and range
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,450.,400,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.channel,"rrv_residue_ErfExp"+label+"_"+self.channel,0.,0.,1.);
            rrv_high_ErfExp    = RooRealVar("rrv_high_ErfExp"+label+"_"+self.channel,"rrv_high_ErfExp"+label+"_"+self.channel,1.,0.,400);
            rrv_high_ErfExp.setConstant(kTRUE);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)* TMath::Power( ((1+TMath::Erf((%s-%s)/%s))/2.), %s )"%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(),rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName(), rrv_high_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_high_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## User1 function 
        if in_model_name == "User1":
            print "########### User 1 Pdf  for mlvj fit ############"
            rrv_p0     = RooRealVar("rrv_p0_User1"+label+"_"+self.channel,"rrv_p0_User1"+label+"_"+self.channel, 31.5438, 10, 90);
            if self.wtagger_label=="HP": #change this!
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -4.93731, -9, -2);
            else:
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel,"rrv_p1_User1"+label+"_"+self.channel, -2, -4, 0.);
            model_pdf=RooUser1Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

        ## Exp+Gaus or mj spectrum
        if in_model_name == "ExpGaus":
            print "########### Exp + Gaus for mj  fit  ############"
            rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,0.05,-0.2,0.2);
            exp             = ROOT.RooExponential("exp"+label+"_"+self.channel,"exp"+label+"_"+self.channel,rrv_x,rrv_c_Exp);

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,78,88);
            rrv_sigma1_gaus = RooRealVar("rrv_smgma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,10);
            rrv_high        = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.5,0.,1.);
            gaus            = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            model_pdf       = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))

        ## Erf*Exp + 2Gaus for mj spectrum
        if in_model_name == "2Gaus_ErfExp":

            print "########### 2Gaus + Erf*Exp for mj fit  ############"
            mean1_tmp      = 8.01411e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#, deltamean_tmp, deltamean_tmp);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel,"rrv_frac_2gaus"+label+"_"+self.channel,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp     = -0.0461211 ; c0_tmp_err     = 6.83e-03;
            offset_tmp = 7.9350e+01  ; offset_tmp_err = 9.35e+00;
            width_tmp  = 3.3083e+01  ; width_tmp_err  = 2.97e+00;

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)

        ## 2Gaus+2Gaus for VV mj spectrum -> WZ and WW
        if in_model_name == "2_2Gaus":

            print "########### 2Gaus +2Gaus for mj fit  ############"
            mean1_tmp      = 9.0141e+01; mean1_tmp_err      = 1.63e-01;
            #mean1_tmp      = 7.3141e+01; mean1_tmp_err      = 1.63e-01; @@@ JEN            
	    deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
#            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;
            sigma1_tmp     = 7.3569e+00; sigma1_tmp_err     = 3.99e-02;

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.channel,"rrv_shift"+label+"_"+self.channel,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-14, mean1_tmp+14);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-14,sigma1_tmp+14 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,0.,-100,10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, 0.0, 1.0);
            #rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.channel,"rrv_frac1"+label+"_"+self.channel,frac_tmp, frac_tmp-frac_tmp_err*2, frac_tmp+frac_tmp_err*2);
            gausguas_1 =RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.channel,"@0-@1",RooArgList(rrv_mean1_gaus, rrv_shift));
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.channel,"@0-@1",RooArgList(rrv_mean2_gaus, rrv_shift));
            gaus3 = RooGaussian("gaus3"+label+"_"+self.channel,"gaus3"+label+"_"+self.channel, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
            gaus4 = RooGaussian("gaus4"+label+"_"+self.channel,"gaus4"+label+"_"+self.channel, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);
            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.channel+mass_spectrum,"gausguas_2"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,0.74,0.0,1.0)#,0.5,1.0)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)
            
        ## return the pdf
        getattr(self.workspace4fit_,"import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)


    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500, interpolate=0):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      ## call the make RooAbsPdf method
      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,area_init_value,0.,1e7);
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList,ismc_wjet)
      print "######## Model Pdf ########"        
      model_pdf.Print();
      rrv_number.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print "######## Model Extended Pdf ########"        

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4fit_,"import")(rrv_number)
      getattr(self.workspace4fit_,"import")(model) 
      self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
      ## return the total extended pdf
      return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);

    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0, massscale=""):

        print "\n\n############### Fit mlvj in mj sideband: ",label," ",mlvj_region,"  ",mlvj_model," ##################"
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_data_mlvj = self.workspace4fit_.data("rdataset_data_xww%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_mlvj_Model("_sb_lo");
        number_VV_sb_lo_mlvj    = self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel))
	print "------> number_VV_sb_lo_mlvj = ",number_VV_sb_lo_mlvj.Print()
        model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_sb_lo");
        number_TTbar_sb_lo_mlvj = self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel))
	print "------> number_TTbar_sb_lo_mlvj = ",number_TTbar_sb_lo_mlvj.Print()
        model_STop_backgrounds  = self.get_STop_mlvj_Model("_sb_lo");
        number_STop_sb_lo_mlvj  = self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel))
	print "------> number_STop_sb_lo_mlvj = ",number_STop_sb_lo_mlvj.Print()
	print "\n\n","*"*20,"\n\n"
	print "\t\tCheck normalization of MJ and MLVJ"
	print "\n\n","*"*20,"\n\n"
	number_VV_sb_lo_mj = self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel))
	number_TTbar_sb_lo_mj = self.workspace4fit_.var("rrv_number_TTbar_xww%s_%s_mj"%(massscale,self.channel))
	number_STop_sb_lo_mj = self.workspace4fit_.var("rrv_number_STop_xww%s_%s_mj"%(massscale,self.channel))
	print "------> number_VV_sb_lo_mj = ",number_VV_sb_lo_mj.Print()
	print "------> number_TTbar_sb_lo_mj = ",number_TTbar_sb_lo_mj.Print()
	print "------> number_STop_sb_lo_mj = ",number_STop_sb_lo_mj.Print()
	print "\n\n","*"*20,"\n\n"


        self.workspace4fit_.var("rrv_number_TTbar_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_STop_xww_sb_lo_%s_mlvj"%(self.channel)).Print();
        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_lo_from_fitting"%(label), mlvj_model,"_mlvj");
	print "\n\n ===\t W+jet model \t== \n"
        model_pdf_WJets.Print();
	model_pdf_WJets.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "==="
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb_lo = self.workspace4fit_.var("rrv_number%s_sb_lo_%s_mlvj"%(label,self.channel)).clone("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel));
        model_WJets =RooExtendPdf("model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),"model%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel),model_pdf_WJets,number_WJets_sb_lo);
	print "DEBUG : 1: Print W+jet model before fit: ==="
        model_pdf_WJets.Print();
	print "\n\n====\n\n"
	model_pdf_WJets.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "\n\n====\n\n"
	print "normalization of W+jet : "
        number_WJets_sb_lo.Print()
	print "\n\n====\n\n"
	print "\n\n====\t Print exteneded w+jet model before fit"
	model_WJets.Print()
	model_WJets.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "\n\n====\n\n"

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_mlvj"%(label,mlvj_region),"model_data%s%s_mlvj"%(label,mlvj_region),RooArgList(model_WJets,model_VV_backgrounds ));
	print "\n\n====\n\n"
	print "\n\n====\t Print mWW model before fit (data and all mc added together)"
	model_data.Print()
	model_data.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "\n\n====\n\n"
        
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.NumCPU(4));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.NumCPU(4));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.NumCPU(4));
	print "\n\n=== \t Print results after fit \t ==="
        rfresult.Print();
	print "\n\n=== \t Print covariance matrix \t ==="
        rfresult.covarianceMatrix().Print();
	print "\n\n====\n\n"
	print "\n\n====\t Print mWW model after fit (data and all mc added together)"
	model_data.Print()
	model_data.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "\n\n====\n\n"
	fitresultsmlvj.append(rfresult)
        getattr(self.workspace4fit_,"import")(model_data)

	print "\n\n===\t Print W+jet model : \n\n"
        model_WJets.Print();
	print "\n\n=== \n\n"
	print "=== \t Print parameters (WJets) \t==="
        model_WJets.getParameters(rdataset_data_mlvj).Print("v");
	print "== GetVariables ="
	params = model_WJets.getVariables();
	params.Print("v");
	print "=== print model model_pdf_WJets ==="
	model_pdf_WJets.Print();
	print "pars..."
	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
	#rrv_c_Exp_WJets0_xww_sb_lo_from_fitting_em = params.find("rrv_c_Exp_WJets0_xww_sb_lo_from_fitting_em")
	#rrv_c_Exp_WJets0_xww_sb_lo_from_fitting_em.setVal(rrv_c_Exp_WJets0_xww_sb_lo_from_fitting_em.getVal()+rrv_c_Exp_WJets0_xww_sb_lo_from_fitting_em.getError())
	print "\n\n=== \n\n"
        model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
	print "\n\n=== \n\n"
	print "\n\n== Print all MC model after Fit == \n\n"
	print "#### VV model "
	model_VV_backgrounds.Print();
	model_VV_backgrounds.getParameters(ROOT.RooArgSet(rrv_mass_lvj)).Print("v");
	print "\n\n----------\tFinish----\n\n"


        self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label,self.channel)).getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),"rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel),
                                                 self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getVal()+
                                                 self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getVal() );

        rrv_number_data_sb_lo_mlvj.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number%s_sb_lo_from_fitting_%s_mlvj"%(label,self.channel)).getError()+
                                                        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()*
                                                        self.workspace4fit_.var("rrv_number_VV_xww_sb_lo_%s_mlvj"%(self.channel)).getError()));

        getattr(self.workspace4fit_,"import")(rrv_number_data_sb_lo_mlvj)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label,self.channel,self.channel,self.channel,self.channel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(self.channel,self.channel,self.channel)),RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(self.color_palet["STop"]),RooFit.VLines());

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_sb_lo_from_fitting_%s_mlvj,model_VV_xww_sb_lo_%s_mlvj"%(label,self.channel,self.channel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_VV_xww_sb_lo_%s_mlvj"%(self.channel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_xww_sb_lo_%s_mlvj,model_STop_xww_sb_lo_%s_mlvj"%(self.channel,self.channel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_STop_xww_sb_lo_%s_mlvj"%(self.channel)), RooFit.Name("STop_line_invisible"), RooFit.LineColor(self.color_palet["STop"]), RooFit.LineWidth(2), RooFit.VLines());
 
            ### draw the error band 
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4fit_.var("rrv_number_data_xww_sb_lo_%s_mlvj"%(self.channel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            self.getData_PoissonInterval(rdataset_data_mlvj,mplot);

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);  

            ### Add the legend to the plot 
            self.leg=self.legend4Plot(mplot,0,1,0., 0.06, 0.16, 0.);
            mplot.addObject(self.leg)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            print mplot.chiSquare();
            print "#################### nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            ### write the result in the output
            self.file_out.write("\n fit_mlvj_in_Mj_sideband: nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof ) );

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
            parameters_list = model_data.getParameters(rdataset_data_mlvj);
            datahist = rdataset_data_mlvj.binnedClone( rdataset_data_mlvj.GetName()+"_binnedClone",rdataset_data_mlvj.GetName()+"_binnedClone" )
                
            self.draw_canvas_with_pull( rrv_mass_lvj,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), "m_lvj_sb_lo%s"%(label),"",1,1)

	"""
        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sb_lo_from_fitting_mlvj"%(label));
	purity = self.wtagger_label[0]+self.wtagger_label[1]
        Deco      = PdfDiagonalizer("Deco%s_sb_lo_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,purity),wsfit_tmp,rfresult);
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        print"#################### print parameters "
        model_pdf_WJets_deco.Print("v");
	print "=> print diagonalized parameters..."
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("v");
	print "=> Done...."
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_deco);

        #### Call the alpha evaluation in automatic
        #self.get_WJets_mlvj_correction_sb_lo_to_signal_region(label,mlvj_model);
        #self.get_WJets_mlvj_correction_sb_lo_to_signal_region_MCOnly(label,mlvj_model);

        ### Fix the pdf of signal, TTbar, STop and VV in the signal region 
        #if (options.interpolate == False):
        #    self.fix_Model("_%s_xww"%(self.signal_sample),"_signal_region","_mlvj")
        self.fix_Model("_VV_xww","_signal_region","_mlvj")

        ### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, STop, and WJets after the extrapolation via alpha
        #if (options.interpolate == True):
        #    self.get_mlvj_normalization_insignalregion("_%s_xww"%(self.signal_sample),"",1);
        #else:
        #    self.get_mlvj_normalization_insignalregion("_%s_xww"%(self.signal_sample));
        self.get_mlvj_normalization_insignalregion("_VV_xww");
        #self.get_mlvj_normalization_insignalregion(label,"model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel));    

	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),rrv_mass_lvj)
	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_auto.root")
	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),132)
	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_132bin.root")
	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),80)
	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_80bin.root")
	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),60)
	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_60bin.root")


	##
	#
	#		Get decorrelated parameters
	#
	##
	initialPars = []
	initialParNames = []
	initialParErrors = []
        parameters = model_pdf_WJets.getParameters(RooArgSet(rrv_mass_lvj));
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
	    initialPars.append(param.getVal())
	    initialParErrors.append(param.getError())
	    #initialParNames.append(param.getName())
	    print "Initial pars = ",param.getVal(), " +/- ", param.getError()
            param=par.Next()
        
	#while (param):
	cov = rfresult.covarianceMatrix();
	print "==> Covariance matrix :"
	cov.Print()

	parToVary = cov.GetNrows()

	print "Parameters to vary = ",parToVary

	eigen = ROOT.TMatrixDSymEigen(cov)

	EigenVector_matrix = eigen.GetEigenVectors()
	EigenValues = eigen.GetEigenValues()

	EigenVector_matrix.Print()
	EigenValues.Print()

	eigenvectors =[]
	for i in xrange(EigenVector_matrix.GetNrows()):
    		eigenvector = ROOT.TVectorD(EigenVector_matrix.GetNrows())
    		for j in xrange(EigenVector_matrix.GetNrows()):
    	    		eigenvector[j] = EigenVector_matrix[j][i]
    		#eigenvector.Print()
    		eigenvectors.append(eigenvector)

	systFunctions = []
	names = []
	
	eigenvector = ROOT.TVectorD()
	
	for k in xrange(EigenVector_matrix.GetNrows()):
	    if k != 1:
	    	eigenvector = eigenvectors[k]
	    	norm = eigenvector.Norm2Sqr()
	    	eigenvalue = EigenValues[k]
	    	sigma = math.sqrt(ROOT.TMath.Abs(eigenvalue))

	    	# compute unit vector in direction of i-th Eigenvector
	    	# eigenvector_unit = (1.0/float(norm))*eigenvector
	    	eigenvector_unit = eigenvector
	    	
	    	upPars = []
	    	downPars = []
	
	    	#print "initial pars = ",len(initialPars)
	    	print "DEBUG: 1: ",k,"=="*10
	    	print "eigenvector:"
	    	eigenvector.Print()
	    	print "Norm = ",norm
	    	print "Eigenvalue : ",eigenvalue
	    	#eigenvalue.Print()
	    	print "=="*10
	    	print "sigma = ",sigma
	    	for i in xrange(len(initialPars)):
	    	    newParUp = initialPars[i] + sigma*eigenvector_unit[2*i]
	    	    print "initialPars[",i,"] = ",initialPars[i],"\t sigma = ",sigma,"\teigenvector_unit = ",eigenvector_unit[2*i]
	    	    upPars.append(newParUp)
	    	    print "==> UP: ",initialPars[i],"\t",newParUp,"\t",sigma,"\t",eigenvector_unit[i]
	    	for i in xrange(len(initialPars)):
	    	    newParDown = initialPars[i] - sigma*eigenvector_unit[i]
	    	    print "initialPars[",i,"] = ",initialPars[i],"\t sigma = ",sigma,"\teigenvector_unit = ",eigenvector_unit[2*i]
	    	    downPars.append(newParDown)
	    	    print "==> Down: ",initialPars[i],"\t",newParDown,"\t",sigma,"\t",eigenvector_unit[i]
	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
	    	print "--"*21
	    	print "Up pars..."
	    	print "--"*21
	    	parameters_list = model_pdf_WJets.getParameters(rdataset_data_mlvj);
	    	par=parameters_list.createIterator(); par.Reset();
	    	param=par.Next()
	    	parCount=0
	    	print "==== Par Setting up ===="
	    	while (param):
	    	    param.setVal(upPars[parCount])
	    	    #print param.getName(),"\t",upPars[parCount]
	    	    print "upPars = ", upPars[parCount]
	    	    param.Print()
	    	    parCount+=1
	    	    param=par.Next()
	    	print "=="*10
	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");

	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),rrv_mass_lvj)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_auto_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),108)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_108bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),47)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_82bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),40)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_52bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),30)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_36bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),20)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_24bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),10)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_12bin_Up_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),4)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_4bin_Up_"+str(k)+".root")

	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
	    	print "--"*21
	    	print "Down pars..."
	    	print "--"*21
	    	parameters_list = model_pdf_WJets.getParameters(rdataset_data_mlvj);
	    	par=parameters_list.createIterator();
	    	par.Reset();
	    	param=par.Next()
	    	parCount=0
	    	while (param):
	    	    param.setVal(downPars[parCount])
	    	    parCount+=1
	    	    param=par.Next()
	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");

	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),rrv_mass_lvj)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_auto_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),108)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_108bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),47)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_82bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),40)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_52bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),30)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_36bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),20)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_24bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),10)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_12bin_Down_"+str(k)+".root")
	    	hist = model_pdf_WJets.createHistogram(rrv_mass_lvj.GetName(),4)
	    	hist.SaveAs("wjetmodel_Ex_"+label+"_"+mlvj_region+"_"+mlvj_model+"_4bin_Down_"+str(k)+".root")
	    	#

	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");
	    	print "--"*21
	    	print "Reset pars..."
	    	print "--"*21
	    	parameters_list = model_pdf_WJets.getParameters(rdataset_data_mlvj);
	    	par=parameters_list.createIterator();
	    	par.Reset();
	    	param=par.Next()
	    	parCount=0
	    	while (param):
	    	    param.setVal(initialPars[parCount])
	    	    param=par.Next()
	    	    parCount+=1
	    	model_pdf_WJets.getParameters(rdataset_data_mlvj).Print("v");




	    	print "\n\n** \t",i


    ##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
    def get_mlvj_normalization_insignalregion(self, label, model_name="", interpolate=0):
        
        print "############### get mlvj normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signal_region"+"_"+self.channel+"_mlvj");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj");

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj) );
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("signal_region"));
        highMassInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("high_mass"));
        
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        highMassInt_val = highMassInt.getVal()/fullInt_val 

        ## integal in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        if (interpolate==0):
            print "####### Events Number in MC Dataset:"
            self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").Print();
            self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").Print();
            print "########## Events Number get from fit:"
            rrv_tmp=self.workspace4fit_.var("rrv_number"+label+"_signal_region"+"_"+self.channel+"_mlvj");
            print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)
            print "Events Number in Signal Region from fitting 1: %s"%(rrv_tmp.getVal())
            print "Events Number in Signal Region from fitting 2: %s"%(signalInt_val)
        else:
            rrv_tmp=self.workspace4fit_.var("rrv_number"+label+"_signal_region"+"_"+self.channel+"_mlvj");
            print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal())

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        if (interpolate==0):
            self.file_out.write( "\nEvents Number in All Region from dataset : %s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal(),self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getError()) )
            self.file_out.write( "\nEvents Number in Signal Region from dataset: %s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal(),self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getError()) )
#            self.file_out.write( "\nRatio signal_region/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj").getVal() ) )
            self.file_out.write( "\nEvents Number in All Region from fitting : %s +/- %s\n"%(rrv_tmp.getVal(),rrv_tmp.getError()) )
	    # check below error if its correct (CHECK 1)
            self.file_out.write( "\nEvents Number in Signal Region from fitting: %s +/- %s\n"%(rrv_tmp.getVal()*signalInt_val,rrv_tmp.getError()) )
            self.file_out.write( "\nEvents Number in High Mass Region from fitting: %s +/- %s\n"%(rrv_tmp.getVal()*highMassInt_val,rrv_tmp.getError()) )
            self.file_out.write( "\nRatio signal_region/all_range from fitting :%s "%(signalInt_val ) )
        else:
            self.file_out.write( "\nEvents Number in Signal Region from fitting: %s +/- %s\n"%(rrv_tmp.getVal(),rrv_tmp.getError()) )

#LUCA
        if (interpolate==0):            
            if not self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj"):
                rrv_number_fitting_signal_region_mlvj = RooRealVar("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_fitting_signal_region"+label+"_"+
                                                                self.channel+"_mlvj", rrv_tmp.getVal()*signalInt_val );
                getattr(self.workspace4fit_,"import")(rrv_number_fitting_signal_region_mlvj);
            else :
                self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val);

        else:
                rrv_number_fitting_signal_region_mlvj = RooRealVar("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_fitting_signal_region"+label+"_"+
                                                                self.channel+"_mlvj", rrv_tmp.getVal() );
                getattr(self.workspace4fit_,"import")(rrv_number_fitting_signal_region_mlvj);

        self.workspace4fit_.var("rrv_number_fitting_signal_region"+label+"_"+self.channel+"_mlvj").Print();

    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_mlvj_correction_sb_lo_to_signal_region_MCOnly(self,label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",mlvj_model," ###############";          
        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        rdataset_WJets_sb_lo_mlvj = self.workspace4fit_.data("rdataset4fit%s_sb_lo_%s_mlvj"%(label,self.channel))
        rdataset_WJets_signal_region_mlvj = self.workspace4fit_.data("rdataset4fit%s_signal_region_%s_mlvj"%(label,self.channel))

	print "="*20
	print type(rdataset_WJets_sb_lo_mlvj)
	print type(rdataset_WJets_signal_region_mlvj)
	print "="*20
	alpha = RooFormulaVar("alpha","log(@0)-log(@1)",RooArgList(rdataset_WJets_signal_region_mlvj,rdataset_WJets_sb_lo_mlvj))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("alpha_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor))) ;
        mplot.GetYaxis().SetTitle("F_{W+jets}^{SR,MC},F_{W+jets}^{SB,MC} (arbitrary units)");

        alpha.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha") );

        self.draw_canvas(mplot,"%s/other/"%(self.plotsDir),"alpha_pdf_%s_%s_M_lvj_signal_region_to_sideband"%(label,mlvj_model),0,1,0,1);



    def get_WJets_mlvj_correction_sb_lo_to_signal_region(self,label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",mlvj_model," ###############";          
        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        rdataset_WJets_sb_lo_mlvj = self.workspace4fit_.data("rdataset4fit%s_sb_lo_%s_mlvj"%(label,self.channel))
        rdataset_WJets_signal_region_mlvj = self.workspace4fit_.data("rdataset4fit%s_signal_region_%s_mlvj"%(label,self.channel))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor))) ;
        mplot.GetYaxis().SetTitle("F_{W+jets}^{SR,MC},F_{W+jets}^{SB,MC} (arbitrary units)");

        if mlvj_model=="Exp":
            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Exp%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Exp%s_%s"%(label,self.channel),"rrv_delta_c_Exp%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal(),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                                      self.workspace4fit_.var("rrv_c_Exp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )

            correct_factor_pdf = RooExponential("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);
            
        if mlvj_model=="Pow":

            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Pow%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Pow%s_%s"%(label,self.channel),"rrv_delta_c_Pow%s_%s"%(label,self.channel),0., -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            correct_factor_pdf = RooPowPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);

        if mlvj_model == "ExpSlowFastFall":
            
	    rrv_c_ExpSlowFastFall_sb = self.workspace4fit_.var("rrv_c_ExpSlowFastFall%s_sb_lo_%s"%(label,self.channel));
	    rrv_n_ExpSlowFastFall_sb = self.workspace4fit_.var("rrv_n_ExpSlowFastFall%s_sb_lo_%s"%(label,self.channel));

	    rrv_delta_c_ExpSlowFastFall = RooRealVar("rrv_delta_c_ExpSlowFastFall%s_%s"%(label,self.channel), "rrv_delta_c_ExpSlowFastFall%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpSlowFastFall_sb.getVal(),
					self.workspace4fit_.var("rrv_c_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpSlowFastFall_sb.getVal()-4*rrv_c_ExpSlowFastFall_sb.getError(),
					self.workspace4fit_.var("rrv_c_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpSlowFastFall_sb.getVal()+4*rrv_c_ExpSlowFastFall_sb.getError())
	    rrv_delta_n_ExpSlowFastFall = RooRealVar("rrv_delta_n_ExpSlowFastFall%s_%s"%(label,self.channel), "rrv_delta_n_ExpSlowFastFall%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_n_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpSlowFastFall_sb.getVal(),
	    				self.workspace4fit_.var("rrv_n_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpSlowFastFall_sb.getVal()-4*rrv_n_ExpSlowFastFall_sb.getError(),
	    				self.workspace4fit_.var("rrv_n_ExpSlowFastFall%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpSlowFastFall_sb.getVal()+4*rrv_n_ExpSlowFastFall_sb.getError())
            model_pdf_1 = ROOT.RooExpDan("model_pdf_1"+label+"_"+self.channel,"model_pdf_1"+label+"_"+self.channel,rrv_x,rrv_delta_c_ExpSlowFastFall, rrv_delta_n_ExpSlowFastFall);

            #rrv_c_Gaus = RooRealVar("rrv_c_Gaus"+label+"_"+self.channel,"rrv_c_Gaus"+label+"_"+self.channel,-5.32e-3,-1e-1,-1e-6);
	    rrv_c_Gaus_sb = self.workspace4fit_.var("rrv_c_Gaus%s_sb_lo_%s"%(label,self.channel));
	    rrv_n_Gaus_sb = self.workspace4fit_.var("rrv_n_Gaus%s_sb_lo_%s"%(label,self.channel));
            #rrv_n_Gaus = RooRealVar("rrv_n_Gaus"+label+"_"+self.channel,"rrv_n_Gaus"+label+"_"+self.channel, -64.3, -1e4, 1e5);
	    rrv_delta_c_Gaus = RooRealVar("rrv_delta_c_Gaus%s_%s"%(label,self.channel),"rrv_delta_c_Gaus%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_Gaus_sb.getVal(),
	    				self.workspace4fit_.var("rrv_c_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_Gaus_sb.getVal()-4*rrv_c_Gaus_sb.getError(),
	    				self.workspace4fit_.var("rrv_c_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_Gaus_sb.getVal()+4*rrv_c_Gaus_sb.getError())
	    rrv_delta_n_Gaus = RooRealVar("rrv_delta_n_Gaus%s_%s"%(label,self.channel),"rrv_delta_n_Gaus%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_n_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_Gaus_sb.getVal(),
	    				self.workspace4fit_.var("rrv_n_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_Gaus_sb.getVal()-4*rrv_n_Gaus_sb.getError(),
	    				self.workspace4fit_.var("rrv_n_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_Gaus_sb.getVal()+4*rrv_n_Gaus_sb.getError())

            model_pdf_2 = RooGaussian("model_pdf_2"+label+"_"+self.channel,"model_pdf_2"+label+"_"+self.channel,rrv_x,rrv_delta_c_Gaus, rrv_delta_n_Gaus);

	    rrv_frac_Gaus_sb = self.workspace4fit_.var("rrv_frac_Gaus%s_sb_lo_%s"%(label,self.channel));
	    rrv_delta_frac_Gaus = RooRealVar("rrv_delta_frac_Gaus%s_%s"%(label,self.channel),"rrv_delta_frac_Gaus%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_frac_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_Gaus_sb.getVal(),
	    				self.workspace4fit_.var("rrv_frac_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_Gaus_sb.getVal()-4*rrv_frac_Gaus_sb.getError(),
	    				self.workspace4fit_.var("rrv_frac_Gaus%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_Gaus_sb.getVal()+4*rrv_frac_Gaus_sb.getError())

            correct_factor_pdf = ROOT.RooAddPdf("correct_factor_pdf","correct_factor_pdf",ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(rrv_frac_Gaus_sb));
            #correct_factor_pdf = ROOT.RooAddPdf("correct_factor_pdf","correct_factor_pdf",ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(frac),1);

        if mlvj_model == "Landau":
            
	    rrv_MPV_Landau_sb = self.workspace4fit_.var("rrv_MPV_Landau%s_sb_lo_%s"%(label,self.channel));
	    rrv_Sigma_Landau_sb = self.workspace4fit_.var("rrv_Sigma_Landau%s_sb_lo_%s"%(label,self.channel));

	    rrv_delta_MPV_Landau = RooRealVar("rrv_delta_MPV_Landau%s_%s"%(label,self.channel), "rrv_delta_MPV_Landau%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_MPV_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_MPV_Landau_sb.getVal(),
					self.workspace4fit_.var("rrv_MPV_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_MPV_Landau_sb.getVal()-8*rrv_MPV_Landau_sb.getError(),
					self.workspace4fit_.var("rrv_MPV_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_MPV_Landau_sb.getVal()+8*rrv_MPV_Landau_sb.getError())
	    rrv_delta_Sigma_Landau = RooRealVar("rrv_delta_Sigma_Landau%s_%s"%(label,self.channel), "rrv_delta_Sigma_Landau%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_Sigma_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_Sigma_Landau_sb.getVal(),
	    				self.workspace4fit_.var("rrv_Sigma_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_Sigma_Landau_sb.getVal()-8*rrv_Sigma_Landau_sb.getError(),
	    				self.workspace4fit_.var("rrv_Sigma_Landau%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_Sigma_Landau_sb.getVal()+8*rrv_Sigma_Landau_sb.getError())
            #model_pdf_1 = ROOT.RooLandau("model_pdf_1"+label+"_"+self.channel,"model_pdf_1"+label+"_"+self.channel,rrv_x,rrv_delta_MPV_Landau, rrv_delta_Sigma_Landau);

	    rrv_c_ExpN_sb = self.workspace4fit_.var("rrv_c_ExpN%s_sb_lo_%s"%(label,self.channel));
	    rrv_n_ExpN_sb = self.workspace4fit_.var("rrv_n_ExpN%s_sb_lo_%s"%(label,self.channel));
	    rrv_delta_c_ExpN = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label,self.channel),"rrv_delta_c_ExpN%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpN_sb.getVal(),
	    				self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpN_sb.getVal()-8*rrv_c_ExpN_sb.getError(),
	    				self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_ExpN_sb.getVal()+8*rrv_c_ExpN_sb.getError())
	    rrv_delta_n_ExpN = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label,self.channel),"rrv_delta_n_ExpN%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpN_sb.getVal(),
	    				self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpN_sb.getVal()-8*rrv_n_ExpN_sb.getError(),
	    				self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_ExpN_sb.getVal()+8*rrv_n_ExpN_sb.getError())

            #model_pdf_2 = ROOT.RooExpNPdf("model_pdf_2"+label+"_"+self.channel,"model_pdf_2"+label+"_"+self.channel,rrv_x,rrv_delta_c_ExpN, rrv_delta_n_ExpN);

            rrv_frac_ExpN_sb = self.workspace4fit_.var("rrv_frac_ExpN%s_sb_lo_%s"%(label,self.channel)); 
	    rrv_delta_frac_ExpN = RooRealVar("rrv_delta_frac_ExpN%s_%s"%(label,self.channel),"rrv_delta_frac_ExpN%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_frac_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_ExpN_sb.getVal(),
	    				self.workspace4fit_.var("rrv_frac_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_ExpN_sb.getVal()-8*rrv_frac_ExpN_sb.getError(),
	    				self.workspace4fit_.var("rrv_frac_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_frac_ExpN_sb.getVal()+8*rrv_frac_ExpN_sb.getError())

	    rrv_MPV_Landau_sr = RooFormulaVar("rrv_MPV_Landau_sr%s_%s"%(label,self.channel), "@0+@1", RooArgList(rrv_MPV_Landau_sb, rrv_delta_MPV_Landau));
	    rrv_Sigma_Landau_sr = RooFormulaVar("rrv_Sigma_Landau_sr%s_%s"%(label,self.channel), "@0+@1", RooArgList(rrv_Sigma_Landau_sb, rrv_delta_Sigma_Landau));
	    rrv_c_ExpN_sr = RooFormulaVar("rrv_c_ExpN_sr%s_%s"%(label,self.channel), "@0+@1", RooArgList(rrv_c_ExpN_sb, rrv_delta_c_ExpN));
	    rrv_n_ExpN_sr = RooFormulaVar("rrv_n_ExpN_sr%s_%s"%(label,self.channel), "@0+@1", RooArgList(rrv_n_ExpN_sb, rrv_delta_n_ExpN));
	    rrv_frac_ExpN_sr = RooFormulaVar("rrv_frac_ExpN_sr%s_%s"%(label,self.channel), "@0+@1", RooArgList(rrv_frac_ExpN_sb, rrv_delta_frac_ExpN));
            #correct_factor_pdf = ROOT.RooAddPdf("correct_factor_pdf","correct_factor_pdf",ROOT.RooArgList(model_pdf_1,model_pdf_2),ROOT.RooArgList(rrv_delta_frac_ExpN),1);
	    correct_factor_pdf = ROOT.RooAlpha4LandauPlusGaus("correct_factor_pdf","correct_factor_pdf",rrv_x, rrv_frac_ExpN_sr, rrv_MPV_Landau_sr, rrv_Sigma_Landau_sr, rrv_c_ExpN_sr, rrv_n_ExpN_sr, rrv_frac_ExpN_sb, rrv_MPV_Landau_sb, rrv_Sigma_Landau_sb, rrv_c_ExpN_sb, rrv_n_ExpN_sb );

        if mlvj_model=="ExpN":
            rrv_c_sb  = self.workspace4fit_.var("rrv_c_ExpN%s_sb_lo_%s"%(label,self.channel));
            rrv_n_sb  = self.workspace4fit_.var("rrv_n_ExpN%s_sb_lo_%s"%(label,self.channel));
            rrv_delta_c = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label,self.channel),"rrv_delta_c_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()-21*rrv_c_sb.getError(),
                                      self.workspace4fit_.var("rrv_c_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c_sb.getVal()+21*rrv_c_sb.getError() )
            rrv_delta_n = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label,self.channel),"rrv_delta_n_ExpN%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()-21*rrv_n_sb.getError(),
                                      self.workspace4fit_.var("rrv_n_ExpN%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_n_sb.getVal()+21*rrv_n_sb.getError() )

            correct_factor_pdf = RooExpNPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c, rrv_delta_n);
 
        if mlvj_model=="ExpTail":
            rrv_s_sb =self.workspace4fit_.var("rrv_s_ExpTail%s_sb_lo_%s"%(label,self.channel));
            rrv_a_sb =self.workspace4fit_.var("rrv_a_ExpTail%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_s = RooRealVar("rrv_delta_s_ExpTail%s_%s"%(label,self.channel),"rrv_delta_s_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()-4*rrv_s_sb.getError(),
                                      self.workspace4fit_.var("rrv_s_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_s_sb.getVal()+4*rrv_s_sb.getError() )
            rrv_delta_a = RooRealVar("rrv_delta_a_ExpTail%s_%s"%(label,self.channel),"rrv_delta_a_ExpTail%s_%s"%(label,self.channel),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()-4*rrv_a_sb.getError(),
                                      self.workspace4fit_.var("rrv_a_ExpTail%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_a_sb.getVal()+4*rrv_a_sb.getError() )
                     
            rrv_a_sr = RooFormulaVar("rrv_a_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_a_sb, rrv_delta_a ) );
            rrv_s_sr = RooFormulaVar("rrv_s_ExpTail_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_s_sb, rrv_delta_s ) );

            correct_factor_pdf = RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_s_sr, rrv_a_sr, rrv_s_sb, rrv_a_sb);

        if mlvj_model=="ErfExp_v1":

            rrv_c_sb       = self.workspace4fit_.var("rrv_c_ErfExp%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb  = self.workspace4fit_.var("rrv_offset_ErfExp%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb   = self.workspace4fit_.var("rrv_width_ErfExp%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfExp%s_%s"%(label,self.channel),"rrv_delta_c_ErfExp%s_%s"%(label,self.channel),0.,
                                           -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfExp%s_%s"%(label,self.channel),"rrv_delta_offset_ErfExp%s_%s"%(label,self.channel),0.,
                                           -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width = RooRealVar("rrv_delta_width_ErfExp%s_%s"%(label,self.channel),"rrv_delta_width_ErfExp%s_%s"%(label,self.channel),0.,
                                          -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb, rrv_x.getMin(), rrv_x.getMax());
            
        if mlvj_model=="ErfPow":

            rrv_c_sb      = self.workspace4fit_.var("rrv_c_ErfPow%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfPow%s_%s"%(label,self.channel),"rrv_delta_c_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPow%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width  = RooRealVar("rrv_delta_width_ErfPow%s_%s"%(label,self.channel),"rrv_delta_width_ErfPow%s_%s"%(label,self.channel),0.,
                                          -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());
            
            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPow2":

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_c2_sb     = self.workspace4fit_.var("rrv_c2_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow2%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow2%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c0      = RooRealVar("rrv_delta_c0_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_c0_ErfPow2%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal(),
					self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()-8*rrv_c0_sb.getError(),
					self.workspace4fit_.var("rrv_c0_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()+8*rrv_c0_sb.getError())
            rrv_delta_c1      = RooRealVar("rrv_delta_c1_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_c1_ErfPow2%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal(),
					self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()-8*rrv_c1_sb.getError(),
					self.workspace4fit_.var("rrv_c1_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()+8*rrv_c1_sb.getError())
            rrv_delta_c2      = RooRealVar("rrv_delta_c2_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_c2_ErfPow2%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_c2_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c2_sb.getVal(),
					self.workspace4fit_.var("rrv_c2_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c2_sb.getVal()-8*rrv_c2_sb.getError(),
					self.workspace4fit_.var("rrv_c2_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c2_sb.getVal()+8*rrv_c2_sb.getError())
            rrv_delta_offset  = RooRealVar("rrv_delta_offset_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPow2%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal(),
					self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()-8*rrv_offset_sb.getError(),
					self.workspace4fit_.var("rrv_offset_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()+8*rrv_offset_sb.getError())
            rrv_delta_width   = RooRealVar("rrv_delta_width_ErfPow2%s_%s"%(label,self.channel),"rrv_delta_width_ErfPow2%s_%s"%(label,self.channel),
	    				self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal(),
					self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()-8*rrv_width_sb.getError(),
					self.workspace4fit_.var("rrv_width_ErfPow2%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()+8*rrv_width_sb.getError())
            
            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_c2_sr     = RooFormulaVar("rrv_c2_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c2_sb, rrv_delta_c2 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = ROOT.RooAlpha4ErfPow3Pdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr,rrv_c2_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb,rrv_c2_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPowExp": ## take initial value from what was already fitted in the SR

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPowExp%s_sb_lo_%s"%(label,self.channel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPowExp%s_sb_lo_%s"%(label,self.channel));

            rrv_delta_c0  = RooRealVar("rrv_delta_c0_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_c0_ErfPowExp%s_%s"%(label,self.channel),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal(),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                                        self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError() )

            rrv_delta_c1 = RooRealVar("rrv_delta_c1_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_c1_ErfPowExp%s_%s"%(label,self.channel),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal(),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                                       self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError() )

            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_offset_ErfPowExp%s_%s"%(label,self.channel),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal(),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()-4*rrv_offset_sb.getError(),
                                       self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_offset_sb.getVal()+4*rrv_offset_sb.getError())

            rrv_delta_width = RooRealVar("rrv_delta_width_ErfPowExp%s_%s"%(label,self.channel),"rrv_delta_width_ErfPowExp%s_%s"%(label,self.channel),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal(),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()-4*rrv_width_sb.getError(),
                                         self.workspace4fit_.var("rrv_width_ErfPowExp%s_signal_region_%s"%(label,self.channel)).getVal()-rrv_width_sb.getVal()+4*rrv_width_sb.getError() )

            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.channel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowExpPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb);
 
        ### define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label,self.channel));

        model_pdf_sb_lo_WJets         = self.workspace4fit_.pdf("model_pdf%s_sb_lo_%s_mlvj"%(label,self.channel));
        model_pdf_signal_region_WJets = RooProdPdf("model_pdf%s_signal_region_%s_mlvj"%(label,self.channel),"model_pdf%s_signal_region_%s_mlvj"%(label,self.channel) ,model_pdf_sb_lo_WJets,correct_factor_pdf);

        simPdf = RooSimultaneous("simPdf","simPdf",data_category);
        simPdf.addPdf(model_pdf_sb_lo_WJets,"sideband");
        simPdf.addPdf(model_pdf_signal_region_WJets,"signal_region");
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.NumCPU(4));
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.NumCPU(4));
	print "####### print rfresult=simPdf "
        rfresult.Print();
	print "####### print rfresult=simPdf covariance matrix"
        rfresult.covarianceMatrix().Print();
	fitresultsfinal.append(rfresult)
	#fitresultsfinal.append("####### print rfresult=simPdf covariance matrix")
	fitresultsfinal.append(rfresult.covarianceMatrix())
	print "--------------------\n"

        ### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_mlvj"%(label));
        print "############### diagonalizer alpha ";
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label),wsfit_tmp,rfresult);
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf);
        correct_factor_pdf_deco.Print();
        correct_factor_pdf_deco.getParameters(rdataset_WJets_signal_region_mlvj).Print("v");
        getattr(self.workspace4fit_,"import")(correct_factor_pdf_deco);
                     
        ## in case of default Wjets with default shape
        if TString(label).Contains("_WJets0"):

            ### only mc plots in the SB region
            mplot_sb_lo = rrv_x.frame(RooFit.Title("WJets sb low"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
            
            rdataset_WJets_sb_lo_mlvj.plotOn(mplot_sb_lo, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_sb_lo_WJets.plotOn(mplot_sb_lo);
            mplot_pull_sideband = self.get_pull(rrv_x,mplot_sb_lo);
            parameters_list     = model_pdf_sb_lo_WJets.getParameters(rdataset_WJets_sb_lo_mlvj);
            mplot_sb_lo.GetYaxis().SetRangeUser(1e-2,mplot_sb_lo.GetMaximum()*1.2);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot_sb_lo.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            print mplot_sb_lo.chiSquare();
            print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot_sb_lo.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            datahist = rdataset_WJets_sb_lo_mlvj.binnedClone( rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone",rdataset_WJets_sb_lo_mlvj.GetName()+"_binnedClone" )

            self.draw_canvas_with_pull( rrv_x,datahist,mplot_sb_lo, mplot_pull_sideband,ndof,parameters_list,"%s/other/"%(self.plotsDir), "m_lvj%s_sb_lo_sim"%(label),"",1,1)
            ### only mc plots in the SR region
            mplot_signal_region = rrv_x.frame(RooFit.Title("WJets sr"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));

            rdataset_WJets_signal_region_mlvj.plotOn(mplot_signal_region, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_signal_region_WJets.plotOn(mplot_signal_region);
            mplot_pull_signal_region = self.get_pull(rrv_x, mplot_signal_region);
            parameters_list = model_pdf_signal_region_WJets.getParameters(rdataset_WJets_signal_region_mlvj);
            mplot_signal_region.GetYaxis().SetRangeUser(1e-2,mplot_signal_region.GetMaximum()*1.2);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot_signal_region.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            print mplot_signal_region.chiSquare();
            print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot_signal_region.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );

            self.draw_canvas_with_pull( rrv_x,datahist,mplot_signal_region, mplot_pull_signal_region,ndof,parameters_list,"%s/other/"%(self.plotsDir), "m_lvj%s_signal_region_sim"%(label),"",1,1);

        ### Total plot shape in sb_lo, sr and alpha
        model_pdf_sb_lo_WJets.plotOn(mplot,RooFit.Name("Sideband"));
        model_pdf_signal_region_WJets.plotOn(mplot, RooFit.LineColor(kRed), RooFit.Name("Signal Region"));
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha") );

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha: Alternate PS") );

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha: Alternate Function") );

        paras=RooArgList();

        if mlvj_model=="ErfExp_v1" or mlvj_model=="ErfPow" or mlvj_model=="2Exp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig4"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig5"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="ErfPow2" or mlvj_model=="ErfPowExp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig4"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig5"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig6"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig7"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="Exp" or mlvj_model=="Pow":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));

        if mlvj_model=="ExpN" or mlvj_model=="ExpTail" or mlvj_model=="Pow2":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig0"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig1"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig2"%(label,self.channel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_13TeV_eig3"%(label,self.channel, self.wtagger_label) ));
        
        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band ar 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,2 ,mplot,kGreen+2,"F",3002,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha_invisible #pm",20,400);
            
        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha_invisible") );

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_xww_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        elif TString(label).Contains("_WJets01") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj_13TeV"%(self.channel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        ### Add the legend
        self.leg=self.legend4Plot(mplot,0,0, -0.01, -0.14, 0.01, -0.06, 0., True);
        mplot.addObject(self.leg);
        
        ## set the Y axis in arbitrary unit 
        if self.signal_sample=="ggH600" or self.signal_sample=="ggH700": tmp_y_max=0.25
        else: tmp_y_max=0.28
        mplot.GetYaxis().SetRangeUser(0,tmp_y_max);

        #### Draw another axis with the real value of alpha
        model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x))
        model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))
        correct_factor_pdf_deco.getVal(RooArgSet(rrv_x))
        tmp_alpha_ratio = ( model_pdf_signal_region_WJets.getVal(RooArgSet(rrv_x))/model_pdf_sb_lo_WJets.getVal(RooArgSet(rrv_x)) );
        tmp_alpha_pdf   = correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)) * mplot.getFitRangeBinW(); ## value of the pdf in each point
        tmp_alpha_scale = tmp_alpha_ratio/tmp_alpha_pdf;
	print "CHECK 1: tmp_alpha_scale = ",tmp_alpha_scale,"\ttmp_alpha_pdf = ",tmp_alpha_pdf


        #add alpha scale axis
        axis_alpha=TGaxis( rrv_x.getMax(), 0, rrv_x.getMax(), tmp_y_max, 0, tmp_y_max*tmp_alpha_scale, 510, "+L" ); #-,-+,+,L
        axis_alpha.SetTitle("#alpha");
        axis_alpha.SetTitleOffset(0.65);
        axis_alpha.SetTitleSize(0.05);
        axis_alpha.SetLabelSize(0.045);
        axis_alpha.SetTitleFont(42);
        axis_alpha.SetLabelFont(42);
        #axis_alpha.RotateTitle(1);
        mplot.addObject(axis_alpha);

        self.draw_canvas(mplot,"%s/other/"%(self.plotsDir),"correction_pdf%s_%s_M_lvj_signal_region_to_sideband"%(label,mlvj_model),0,1,0,1);

        correct_factor_pdf_deco.getParameters(rdataset_WJets_sb_lo_mlvj).Print("v");
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco = self.workspace4fit_.pdf("model_pdf%s_sb_lo_from_fitting_%s_mlvj_Deco%s_sb_lo_from_fitting_%s_%s_mlvj_13TeV"%(label,self.channel,label, self.channel,self.wtagger_label[0]+self.wtagger_label[1]));
        model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco.Print("v");

        ### Wjets shape in the SR correctedfunction * sb 
        model_pdf_WJets_signal_region_after_correct_mlvj = RooProdPdf("model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel),"model_pdf%s_signal_region_%s_after_correct_mlvj"%(label,self.channel),model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj_13TeV"%(label,self.channel,self.wtagger_label)) );
        model_pdf_WJets_signal_region_after_correct_mlvj.Print("v")
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_signal_region_after_correct_mlvj)

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).Print();
        self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).Print();
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setVal(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).getVal());
        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setError(self.workspace4fit_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)).getError());

        self.workspace4fit_.var("rrv_number%s_signal_region_%s_mlvj"%(label,self.channel)).setConstant(kTRUE);

    ##### Counting of the events of each component in the signal region taking the lavel for the model
    def get_mj_normalization_insignalregion(self, label):
        print "################## get mj normalization ",label," ################## ";
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        model      = self.workspace4fit_.pdf("model"+label+"_"+self.channel+"_mj");

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        sb_loInt  = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("sb_lo"));
        signalInt = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
        sb_hiInt  = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("sb_hi"));
        
        fullInt_val     = fullInt.getVal()
        sb_loInt_val    = sb_loInt.getVal()/fullInt_val
        sb_hiInt_val    = sb_hiInt.getVal()/fullInt_val
        signalInt_val   = signalInt.getVal()/fullInt_val

        print "########### Events Number in MC Dataset: #############"
        self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").Print();

        print "########### Events Number get from fit: ##############"
        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj");
        rrv_tmp.Print();
        print "Events Number in sideband_low :%s"%(rrv_tmp.getVal()*sb_loInt_val)
        print "Events Number in Signal Region:%s"%(rrv_tmp.getVal()*signalInt_val)
        print "Events Number in sideband_high:%s"%(rrv_tmp.getVal()*sb_hiInt_val)
        print "Total Number in sidebands :%s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val) )
        print "Ratio signal_region/sidebands :%s"%(signalInt_val/(sb_loInt_val+sb_hiInt_val) )

        ##### Save numbers in the output text file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in sideband_low from dataset:%s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal() , self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getError()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset:%s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal() ,self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getError() ) )
        self.file_out.write( "\nEvents Number in sideband_high from dataset:%s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal() ,self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getError()) )
        self.file_out.write( "\nTotal Number in sidebands from dataset:%s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal() ,self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getError()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getError()) )
        if (self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()) != 0:
	   self.file_out.write( "\nRatio signal_region/sidebands from dataset:%s "%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()) ) )
	   #self.file_out.write( "\nRatio signal_region/sidebands from dataset:%s +/- %s"%(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()) ,self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal())*(self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getError()/self.workspace4fit_.var("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj").getVal()+(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getError()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getError())/(self.workspace4fit_.var("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj").getVal()))) )

        self.file_out.write( "\nEvents Number in sideband_low from fitting:%s +/- %s"%(rrv_tmp.getVal()*sb_loInt_val, rrv_tmp.getError()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting:%s +/- %s"%(rrv_tmp.getVal()*signalInt_val, rrv_tmp.getError()) )
        self.file_out.write( "\nEvents Number in sideband_high from fitting:%s +/- %s"%(rrv_tmp.getVal()*sb_hiInt_val, rrv_tmp.getError()) )
        self.file_out.write( "\nTotal Number in sidebands from fitting:%s +/- %s"%(rrv_tmp.getVal()*(sb_loInt_val+sb_hiInt_val) , rrv_tmp.getError() ) )
        self.file_out.write( "\nRatio signal_region/sidebands from fitting:%s "%(signalInt_val/(sb_loInt_val+sb_hiInt_val)) )

    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in signal_region

        print "############### Fit mj Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0_xww");
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets01_xww");

        ## take the normalization numbers
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_xww_in_mj_signal_region_from_fitting_%s"%(self.channel));
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_xww_in_mj_signal_region_from_fitting_%s"%(self.channel));
        rrv_WJets0.Print();
        rrv_WJets01.Print();
        
        ### total uncertainty combining the result with two different shapes
        total_uncertainty = TMath.Sqrt( TMath.Power(rrv_WJets0.getError(),2) + TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2) );
        rrv_WJets0.setError(total_uncertainty);
	print "Total uncertainty combining the result of two different shapes : "
        rrv_WJets0.Print();
	print "*"*20

    #### make the mj sideband fit on data ti get the Wjets normaliztion 
    def fit_WJetsNormalization_in_Mj_signal_region(self,label,massscale=""): 

        print "############### Fit mj Normalization: ",label," ",massscale," ##################"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj=self.workspace4fit_.data("rdataset_data_xww_%s_mj"%(self.channel))

        ### Fix TTbar, VV and STop
        model_VV    = self.get_VV_mj_Model("_VV_xww"+massscale);
	#sys.exit()
        ## only two parameters are fix, offset and width while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets = self.get_WJets_mj_Model(label);

        ## Total Pdf and fit only in sideband 
        model_data = RooAddPdf("model_data_xww%s_%s_mj"%(massscale,self.channel),"model_data_xww%s_%s_mj"%(massscale,self.channel),RooArgList(model_WJets,model_VV));
	print "\n\n","="*20,"==\t MODEL PRINT ","\n\n"
	model_data.Print("v");
	print "\n\n","="*20,"\n\n"
	model_data.getParameters(rdataset_data_mj).Print("v");
	print "\n\n","="*20,"\n\n"
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(4) );
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(4), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
	fitresultsfinal.append(rfresult)
        getattr(self.workspace4fit_,"import")(model_data);

        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),"rrv_number_data_xww%s_%s_mj"%(massscale,self.channel),
                                         self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getVal()+
                                         self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal());

        rrv_number_data_mj.setError(TMath.Sqrt(
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number_VV_xww%s_%s_mj"%(massscale,self.channel)).getError()+
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()*
                                               self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getError()));
	print "\n\n","==== \t prit number : \n\n"
	rrv_number_data_mj.Print("v")
	print "\n\n","="*20
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);
        
        ## if fit on Wjets default with the default shape
        if TString(label).Contains("_WJets0"):

            ## make the final plot
            mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)));
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ## plot solid style 
            model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
            
	    model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(self.color_palet["STop"]),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
    
            ## plot "dashed" style area
            model_data.plotOn(mplot,RooFit.Name("VV_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("TTbar_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj,model_TTbar_xww_%s_mj"%(label,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("STop_invisible"), RooFit.Components("model%s_%s_mj,model_STop_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(self.color_palet["STop"]),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("WJets_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
    
            ### solid line
            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            
            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_VV_%s_mj"%(label,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            ### dash line
            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_VV_xww_%s_mj"%(label,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

            rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### draw the error band using the sum of all the entries component MC + fit         
            draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### Get the pull and plot it 
            mplot_pull=self.get_pull(rrv_mass_j,mplot);

            ### signal window zone with vertical lines
            lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,mplot.GetMaximum()*0.9); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
            upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
            #lowerLine = TLine(65,0.,65,mplot.GetMaximum()*0.9); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
            middleLine1 = TLine(65,0.,65,mplot.GetMaximum()*0.9); middleLine1.SetLineWidth(2); middleLine1.SetLineColor(kBlack); middleLine1.SetLineStyle(9);
            middleLine2 = TLine(105,0.,105,mplot.GetMaximum()*0.9); middleLine2.SetLineWidth(2); middleLine2.SetLineColor(kBlack); middleLine2.SetLineStyle(9);
            middleLine3 = TLine(125,0.,125,mplot.GetMaximum()*0.9); middleLine3.SetLineWidth(2); middleLine3.SetLineColor(kBlack); middleLine3.SetLineStyle(9);
            #upperLine = TLine(95,0.,95,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
	    mplot.addObject(middleLine1);
	    mplot.addObject(middleLine2);
	   #mplot.addObject(middleLine3);
            mplot.addObject(lowerLine);
            mplot.addObject(upperLine);

	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.5728643,0.1008991,"NDC")
#	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.52,0.1008991,"NDC")
#	    pt = ROOT.TPaveText(0.3592965,0.02847153,0.50,0.1008991,"NDC")
	    pt.SetTextFont(42)
	    pt.SetTextSize(0.04995005)
	    #pt.SetTextAlign(12)
	    pt.SetFillColor(self.color_palet["WJets"])
	    pt.SetBorderSize(0)
	    text = pt.AddText("#leftarrow signal region #rightarrow")
	    text.SetTextFont(62)
	    mplot.addObject(pt);

	    pt1 = ROOT.TPaveText(0.3555276,0.1183816,0.4535176,0.1908092,"NDC")
	    pt1.SetTextFont(42)
	    pt1.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt1.SetFillColor(0)
	    pt1.SetFillStyle(0)
	    pt1.SetBorderSize(0)
	    text = pt1.AddText("ZV")
	    text.SetTextFont(62)
	    text = pt1.AddText("category")
	    text.SetTextFont(62)
#	    mplot.addObject(pt1);

	    pt2 = ROOT.TPaveText(0.4723618,0.1183816,0.5678392,0.1908092,"NDC")
	    pt2.SetTextFont(42)
	    pt2.SetTextSize(0.037)
	    #pt.SetTextAlign(12)
	    pt2.SetFillColor(0)
	    pt2.SetBorderSize(0)
	    pt2.SetFillStyle(0)
	    text = pt2.AddText("WZ")
	    text.SetTextFont(62)
	    text = pt2.AddText("category")
	    text.SetTextFont(62)
#	    mplot.addObject(pt2);
	    	    
            ### legend of the plot
            self.leg = self.legend4Plot(mplot,0,1,0.,0.,0.13,0.02,0,0,1);
            #self.leg = self.legend4Plot(mplot,0,1,-0.10,-0.01,0.10,0.01);
            mplot.addObject(self.leg);
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.8);

            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            print mplot.chiSquare();
            print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            datahist = rdataset_data_mj.binnedClone( rdataset_data_mj.GetName()+"_binnedClone",rdataset_data_mj.GetName()+"_binnedClone" )

            parameters_list = model_data.getParameters(rdataset_data_mj);
            self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), "m_j_sideband%s"%(label),"",1,0,1)

            ### call the function for getting the normalizatio in signal region for data, TTbar, STop, VV and W+jets = label -> store in a output txt file
            self.get_mj_normalization_insignalregion("_data_xww");
            self.get_mj_normalization_insignalregion("_VV_xww");
            self.get_mj_normalization_insignalregion(label);

        #### to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data
        fullInt   = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        signalInt = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
        rrv_number_WJets_in_mj_signal_region_from_fitting = RooRealVar("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),"rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.channel)).getVal()*signalInt_val);

        #### Error on the normalization --> from a dedicated function taking into account shape uncertainty
        rrv_number_WJets_in_mj_signal_region_from_fitting.setError( Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult,"signal_region") );
        print "########## error on the normalization due to shape + norm = %s"%(rrv_number_WJets_in_mj_signal_region_from_fitting.getError());
        self.file_out.write("\n########## error on the normalization due to shape + norm = %s \n"%(rrv_number_WJets_in_mj_signal_region_from_fitting.getError()));
        getattr(self.workspace4fit_,"import")(rrv_number_WJets_in_mj_signal_region_from_fitting);
        rrv_number_WJets_in_mj_signal_region_from_fitting.Print();

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def interpolate_signal_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Interpolate mlvj single MC sample ",in_file_name," ",label,"  ",mlvj_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc,500,1); #call with option "interpolate"
#        getattr(self.workspace4fit_,"import")(model) #needed??

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label,"  ",mlvj_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.channel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc);
	print "\n\n","="*20,"==\t MODEL MLVJ PRINT ","\n\n"
	model.Print("v");
	print "\n\n","="*20,"\n\n"
	model.getParameters(rdataset).Print("v");
	print "\n\n","="*20,"\n\n"

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) , RooFit.NumCPU(4));
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") , RooFit.NumCPU(4));
	print "\n\n","_"*20,"\n\n"
        rfresult.Print();
        rfresult.Print("v");
	fitresultsmlvj.append(rfresult)
	print "\n\n\t\t SAVE HISTOGRAM \n",in_file_name,label,in_range,mlvj_model,"\n\n"
	hist = model.createHistogram(rrv_mass_lvj.GetName(),rrv_mass_lvj)
	hist.SaveAs(in_file_name+"_"+label+"_"+in_range+"_"+mlvj_model+"_auto.root")
	print "\n\n","_"*20,"\n\n"

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{llj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,6,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        nPar = rfresult.floatParsFinal().getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-nPar;
        print mplot.chiSquare();
        print "#################### JENchi2 nPar=%s, chiSquare=%s/%s"%(nPar ,mplot.chiSquare(nPar)*ndof, ndof );
        datahist = rdataset.binnedClone( rdataset.GetName()+"_binnedClone",rdataset.GetName()+"_binnedClone" )
	rdataset.Print()
	
        ## get the pull 
        mplot_pull      = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
        
        self.draw_canvas_with_pull( rrv_mass_lvj, datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, logy);

        for i in range(datahist.numEntries()):
	 datahist.get(i)
	 datahist.weightError(RooAbsData.SumW2)
         print "Bin %i x=%f w = %f we = %f"%(i,datahist.get(i).getRealValue("rrv_mass_lvj"),datahist.weight(),datahist.weightError(RooAbsData.SumW2))
         
        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.channel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) , RooFit.NumCPU(4));
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"), RooFit.NumCPU(4));
            rfresult_pdf.Print();
	    fit_DecoResults.append(rfresult_pdf)


            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.channel+"_mlvj");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.channel+"_"+self.wtagger_label+"_mlvj_13TeV",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
            print "##################### diagonalize ";
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize            
            print "##################### workspace for decorrelation ";
            wsfit_tmp.Print("v");
            print "##################### original  parameters ";
            model_pdf.getParameters(rdataset).Print("v");
            print "##################### original  decorrelated parameters ";
            model_pdf_deco.getParameters(rdataset).Print("v");
            print "##################### original  pdf ";
            model_pdf.Print();
            print "##################### decorrelated pdf ";
            model_pdf_deco.Print();

            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_,"import")(model_pdf_deco);

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(int(rrv_mass_lvj.getBins()/self.narrow_factor)));
            
            if label=="_TTbar_xww" and in_range=="_signal_region":
                
                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name("TTbar_Powheg"),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset = RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_xww_signal_region_%s_mlvj"%(self.channel)).Print();
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.leg = self.legend4Plot(mplot_deco,0); ## add the legend                
            mplot_deco.addObject(self.leg);

            self.draw_canvas( mplot_deco, "%s/other/"%(self.plotsDir), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.channel+"_mlvj").Print();

    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj");
        rdataset_mj.Print();

        ## make the extended model
        model = self.make_Model(label,in_model_name);
	print "\n\n","="*20,"==\t MODEL PRINT ","\n\n"
	model.Print("v");
	print "\n\n","="*20,"\n\n"
	model.getParameters(rdataset_mj).Print("v");
	print "\n\n","="*20,"\n\n"
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE) , RooFit.NumCPU(4));
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") , RooFit.NumCPU(4));
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") , RooFit.NumCPU(4));
	print "\n\n","="*20,"==\t FIT RESULT","\n\n"
        rfresult.Print("v");
	print "\n\n","="*20,"\n\n"
	fitresultsmj.append(rfresult)
	
        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins()/self.narrow_factor)) );
        rdataset_mj.plotOn( mplot, RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,6,"L");
        ## re-draw the dataset
        rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## draw the function
        model.plotOn( mplot );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_mass_j, mplot); 
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        nPar = rfresult.floatParsFinal().getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-nPar;
        print "#################### nPar=%s, nBinX=%s , chiSquare=%s/%s"%(nPar,nBinX,mplot.chiSquare(nPar)*ndof,ndof);

        datahist = rdataset_mj.binnedClone( rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone" )
        parameters_list = model.getParameters(rdataset_mj);
        self.draw_canvas_with_pull( rrv_mass_j,datahist,mplot, mplot_pull,ndof,parameters_list,"%s/m_j_fitting/"%(self.plotsDir), label+in_file_name, in_model_name)

	self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();
	self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).Print();
        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

        #datahist = rdataset_mj.binnedClone( rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone" )
        #Nbin = rrv_mass_j.getBins()
        #rresult_param = rfresult.floatParsFinal()
        #nparameters =  rresult_param.getSize()
        #ChiSquare = model.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson))
        #chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters);
     	#print "nPar=%s, chiSquare=%f/%f"%(nparameters, ChiSquare.getVal()*(Nbin - nparameters), (Nbin - nparameters) )

        #datahist = rdataset_mj.binnedClone(rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone")
        #Nbin = int(rrv_mass_j.getBins());
        #nparameters = rfresult.floatParsFinal().getSize();
        #ChiSquare = model.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
        #chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters);
	#print "nPar=%s, chiSquare=%s/%s"%(nparameters, ChiSquare.getVal()*(Nbin - nparameters), (Nbin - nparameters) );
	
        ##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV 
        par=parameters_list.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
            param=par.Next()

    def IsGoodEvent(self,tree):
    
       keepEvent = False
       if ((options.jetalgo).find('Puppi'))!= -1:
           if ((tree.PuppiAK8_jet_tau2tau1) < 0.55) : keepEvent = True
       else:
           if self.wtagger_label.find('HP') != -1:  
               if tree.jet_tau2tau1 < 0.45 : keepEvent = True
           if self.wtagger_label.find('LP') != -1:  
               if tree.jet_tau2tau1 > 0.45 and tree.jet_tau2tau1 < 0.75: keepEvent = True


       """
       #if the category has a specific purity and a specific nbtag
       if (self.wtagger_label.find('HP') != -1 or self.wtagger_label.find('LP') != -1) and (self.wtagger_label.find('0') != -1 or self.wtagger_label.find('1') != -1 or self.wtagger_label.find('2') != -1):
          if tree.channel == self.categoryID: keepEvent = True
          
       #all HP nbtag categories == HP only
       if self.wtagger_label.find('HP') != -1 and (self.wtagger_label.find('0') == -1 and self.wtagger_label.find('1') == -1 and self.wtagger_label.find('2') == -1):  
          if self.channel == 'el':
             if tree.channel == 0 or tree.channel == 6 or tree.channel == 8 or tree.channel == 10: keepEvent = True
          elif self.channel == 'mu':   
             if tree.channel == 1 or tree.channel == 7 or tree.channel == 9 or tree.channel == 11: keepEvent = True

       #all LP nbtag categories == LP only
       if self.wtagger_label.find('LP') != -1 and (self.wtagger_label.find('0') == -1 and self.wtagger_label.find('1') == -1 and self.wtagger_label.find('2') == -1):  
          if self.channel == 'el':
             if tree.channel == 2 or tree.channel == 12 or tree.channel == 14 or tree.channel == 16: keepEvent = True
          elif self.channel == 'mu':   
             if tree.channel == 3 or tree.channel == 13 or tree.channel == 15 or tree.channel == 17: keepEvent = True
           
       #all purities                     
       if self.wtagger_label.find('ALLP') != -1:   
          if self.channel == 'el':
             #all purities but specific nbtag category
             if self.wtagger_label.find('0') != -1:
                if tree.channel == 6 or tree.channel == 12 or tree.channel == 18: keepEvent = True
             elif self.wtagger_label.find('1') != -1:
                if tree.channel == 8 or tree.channel == 14 or tree.channel == 20: keepEvent = True
             elif self.wtagger_label.find('2') != -1:
                if tree.channel == 10 or tree.channel == 16 or tree.channel == 22: keepEvent = True
             #all purities and all nbtag categories
             else:
                if tree.channel == 0 or tree.channel == 2 or tree.channel == 6 or tree.channel == 12 or tree.channel == 8 or tree.channel == 14 or tree.channel == 10 or tree.channel == 16: keepEvent = True
          elif self.channel == 'mu':
             #all purities but specific nbtag category
             if self.wtagger_label.find('0') != -1:
                if tree.channel == 7 or tree.channel == 13: keepEvent = True
             elif self.wtagger_label.find('1') != -1:
                if tree.channel == 9 or tree.channel == 15: keepEvent = True
             elif self.wtagger_label.find('2') != -1:
                if tree.channel == 11 or tree.channel == 17: keepEvent = True
             #all purities and all nbtag categories
             else:
                if tree.channel == 1 or tree.channel == 3 or tree.channel == 7 or tree.channel == 13 or tree.channel == 9 or tree.channel == 15 or tree.channel == 11 or tree.channel == 17: keepEvent = True
       """
       return keepEvent
                          
    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass ):# to get the shape of m_lvj,jet_mass="jet_mass_pr"

        print "################### get_mj_and_mlvj_dataset : ",in_file_name,"  ",label,"  ##################";

	fileIn_name = TString("root://cmseos.fnal.gov/")+TString(self.file_Directory+in_file_name);
	#fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn = TFile.Open(fileIn_name.Data());
        treeIn = fileIn.Get("otree");
        
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rrv_weight   = RooRealVar("rrv_weight","rrv_weight",-1000. ,10000000.)

        self.BinWidth_mlvj=self.BinWidth_mlvj/self.narrow_factor;
        print self.BinWidth_mlvj,self.narrow_factor;
        print rrv_mass_lvj;
        nbins_mlvj=int((rrv_mass_lvj.getMax()-rrv_mass_lvj.getMin())/self.BinWidth_mlvj);

        rrv_mass_lvj.setBins(nbins_mlvj);
	
	if TString(label).Contains("_RS1G_WW_") or TString(label).Contains("_BulkGraviton_") or TString(label).Contains("_Wprime_WZ_"):
	   bw = self.BinWidth_mlvj
	   nbins_mlvj=int((rrv_mass_lvj.getMax()-rrv_mass_lvj.getMin())/bw);
	   rrv_mass_lvj.setBins(nbins_mlvj)
	   
        ##### dataset of m_j -> scaleed and not scaled to lumi 
        rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        ##### dataset of m_lvj -> scaled and not scaled to lumi in different region
        rdataset_sb_lo_mlvj = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_sb_hi_mlvj = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_sb_lo_mlvj = RooDataSet("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_signal_region_mlvj = RooDataSet("rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_sb_hi_mlvj = RooDataSet("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### categorize the event in sideband and signal region --> combined dataset 

        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");
        combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4fit = RooDataSet("combData4fit"+label+"_"+self.channel,"combData4fit"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        print "###### N entries: ", treeIn.GetEntries()
        hnum_4region=TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_2region=TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total

        tmp_lumi=35867.06;
        tmp_scale_to_lumi=1.;
            
	nnevents = treeIn.GetEntries()
	#nnevents = 30000
	if nnevents > treeIn.GetEntries():
		nnevents = treeIn.GetEntries()
	print "Number of events to run = ",nnevents
        #for i in range(treeIn.GetEntries()):
        for i in range(nnevents):
            if i % 100000 == 0: print "iEntry: ",i
            treeIn.GetEntry(i);

            tmp_scale_to_lumi=treeIn.wSampleWeight*tmp_lumi;

            tmp_jet_mass=getattr(treeIn, jet_mass);

            self.isGoodEvent = 0;   

            # Analysis selection here
	    
            if ((options.jetalgo).find('Puppi') != -1): #Puppi
                if treeIn.mass_llj_PuppiAK8> rrv_mass_lvj.getMin() and treeIn.mass_llj_PuppiAK8<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax():
                	#if tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax():
                	self.isGoodEvent = 1;   
		if (treeIn.l_pt2<=30) : self.isGoodEvent = 0;
		if (treeIn.l_pt1<=50): self.isGoodEvent = 0;
		if (treeIn.dilep_m<76 or treeIn.dilep_m>107): self.isGoodEvent = 0;
		if ((treeIn.type == 0 and abs(treeIn.l_eta1)>=2.4) or (treeIn.type==1 and ((abs(treeIn.l_eta1)>=2.5) or (abs(treeIn.l_eta1)>=1.4442  and abs(treeIn.l_eta1)<=1.566))) ): self.isGoodEvent = 0;
		if ((treeIn.type == 0 and abs(treeIn.l_eta2)>=2.4) or (treeIn.type==1 and ((abs(treeIn.l_eta2)>=2.5) or (abs(treeIn.l_eta2)>=1.4442  and abs(treeIn.l_eta2)<=1.566))) ): self.isGoodEvent = 0;
		if (treeIn.l_charge1 > 0 and treeIn.l_charge2 > 0): self.isGoodEvent = 0;
		if (treeIn.l_charge1 < 0 and treeIn.l_charge2 < 0): self.isGoodEvent = 0;
                if (treeIn.ungroomed_PuppiAK8_jet_pt<=200.) : self.isGoodEvent = 0;
                if (abs(treeIn.ungroomed_PuppiAK8_jet_eta)>=2.4): self.isGoodEvent = 0;
                if (treeIn.PuppiAK8_jet_tau2tau1>=0.55) : self.isGoodEvent = 0;

                if (treeIn.nBTagJet_loose!=0) : self.isGoodEvent = 0;
		
                if (label=="_data_xww") and ((treeIn.PuppiAK8_jet_mass_so_corr > 65.) and (treeIn.PuppiAK8_jet_mass_so_corr < 105.)) : self.isGoodEvent = 0; #BLINDING

            #VBF SELECTION
            if ((options.type).find('vbf') != -1 and treeIn.vbf_maxpt_jj_m<800): self.isGoodEvent=0;
            if ((options.type).find('vbf') != -1 and abs(treeIn.vbf_maxpt_j1_eta-treeIn.vbf_maxpt_j2_eta)<4.0): self.isGoodEvent=0;
            if ((options.type).find('vbf') != -1 and (treeIn.vbf_maxpt_j1_pt<30 or treeIn.vbf_maxpt_j2_pt<30)): self.isGoodEvent=0;
	    
	    #if ((label =="_data" or label =="_data_xww") and treeIn.jet_mass_pr >105 and treeIn.jet_mass_pr < 135 ) : self.isGoodEvent = 0; 

            #VBF SELECTION
            if ((options.type).find('vbf') != -1 and treeIn.vbf_maxpt_jj_m<=800): self.isGoodEvent=0;
            if ((options.type).find('vbf') != -1 and abs(treeIn.vbf_maxpt_j1_eta-treeIn.vbf_maxpt_j2_eta)<=4.0): self.isGoodEvent=0;
            if ((options.type).find('vbf') != -1 and (treeIn.vbf_maxpt_j1_pt<=30 or treeIn.vbf_maxpt_j2_pt<=30)): self.isGoodEvent=0;
	    
            if self.isGoodEvent == 1:
	    	#print "== Good event"
                ### weigh MC events              
		#tmp_event_weight     = treeIn.genWeight*treeIn.wSampleWeight*tmp_lumi*treeIn.pu_Weight*treeIn.trig_eff_Weight*treeIn.id_eff_Weight*treeIn.btag0Wgt; 
                #tmp_event_weight4fit = treeIn.genWeight;
		#tmp_event_weight4fit = tmp_event_weight4fit*treeIn.pu_Weight*treeIn.trig_eff_Weight*treeIn.id_eff_Weight*treeIn.wSampleWeight*tmp_lumi*treeIn.btag0Wgt;
		tmp_event_weight     = treeIn.totalEventWeight_2Lep*treeIn.btag0Wgt*treeIn.pu_Weight*treeIn.L1_Prefweight;
		tmp_event_weight4fit = treeIn.totalEventWeight_2Lep*treeIn.btag0Wgt*treeIn.pu_Weight*treeIn.L1_Prefweight*tmp_scale_to_lumi;
	
                if label =="_data" or label =="_data_xww" :
                    tmp_event_weight=1.;
                    tmp_event_weight4fit=1.;                    
                else:
                    if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
                        tmp_event_weight=tmp_event_weight*self.rrv_wtagger_eff_reweight_forT.getVal();
                    else:
                        tmp_event_weight=tmp_event_weight*self.rrv_wtagger_eff_reweight_forV.getVal();
                
                rrv_mass_lvj.setVal(treeIn.mass_llj_PuppiAK8);

                #if ((tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max)):
                if ((tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max) or (tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max)):
                    rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

                    data_category.setLabel("sideband");
                    combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                    rdataset_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                    
                    data_category.setLabel("signal_region");
                    combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                    hnum_2region.Fill(1,tmp_event_weight);

                    if treeIn.mass_llj_PuppiAK8 >=self.mlvj_signal_min and treeIn.mass_llj_PuppiAK8 <self.mlvj_signal_max:
                        hnum_2region.Fill(0,tmp_event_weight);

                if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                    rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                    
                rrv_mass_j.setVal( tmp_jet_mass );
                rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
                    hnum_4region.Fill(1,tmp_event_weight);

                hnum_4region.Fill(2,tmp_event_weight);

        if not label=="_data" and not label =="_data_xww": ## correct also because events in 4fit dataset were not rescaled in the cycle
	    tmp_scale_to_lumi=tmp_scale_to_lumi;
            #if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
            #    tmp_scale_to_lumi=tmp_scale_to_lumi*self.rrv_wtagger_eff_reweight_forT.getVal();
            #else:
            #    tmp_scale_to_lumi=tmp_scale_to_lumi*self.rrv_wtagger_eff_reweight_forV.getVal();

	tmp_scale_to_lumi = 1.0;
        ### scaler to lumi for MC in 4fit datasets
        rrv_scale_to_lumi=RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi)
        rrv_scale_to_lumi.Print()
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)

        ### prepare m_lvj dataset to be compared with the fit results
        rrv_number_dataset_signal_region_mlvj=RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(1));
        rrv_number_dataset_AllRange_mlvj=RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(2));
        
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj)

        ### import the dataser       
        getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj);
        getattr(self.workspace4fit_,"import")(combData);
        getattr(self.workspace4fit_,"import")(combData4fit);

        ### write in the output 
        self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label,rdataset_signal_region_mlvj.sumEntries()))

        ### prepare m_j dataset
        rrv_number_dataset_sb_lo_mj=RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj=RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_sb_hi_mj=RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj)
                
        getattr(self.workspace4fit_,"import")(rdataset_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj)

        #### print everything
        
        rdataset_sb_lo_mlvj.Print();
        rdataset_signal_region_mlvj.Print();
        rdataset_sb_hi_mlvj.Print();
        rdataset_mj.Print();
        rdataset4fit_sb_lo_mlvj.Print();
        rdataset4fit_signal_region_mlvj.Print();
        rdataset4fit_sb_hi_mlvj.Print();
        rdataset4fit_mj.Print();
        rrv_number_dataset_signal_region_mlvj.Print()
        rrv_number_dataset_AllRange_mlvj.Print()
        rrv_number_dataset_sb_lo_mj.Print()
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_sb_hi_mj.Print()
        combData.Print();
        combData4fit.Print();
        
    #### function to run the selection on data to build the datasets 
    def get_data(self):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.file_data,"_data_xww", self.jetalgo)
        #self.get_mj_and_mlvj_dataset(self.file_data,"_data_xww", "Mjsoftdrop")
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_signal_region_data_xww_%s_mlvj"%(self.channel)).clone("observation_for_counting_xww"))

    #### Define the steps to fit STop MC in the mj and mlvj spectra
    def fit_STop(self):
        print "############################## fit_STop  #################################"
        self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop_xww", self.jetalgo)
        #self.fit_mj_single_MC(self.file_STop_mc,"_STop_xww","2Gaus_ErfExp");
        self.fit_mj_single_MC(self.file_STop_mc,"_STop_xww","ExpGaus");
        #self.fit_mj_single_MC(self.file_STop_mc,"_STop_xww","User1");
        #self.fit_mj_single_MC(self.file_STop_mc,"_STop_xww","ErfExp");
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop_xww","_sb_lo","ExpN", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop_xww","_signal_region","ExpN", 1, 0, 1);
        print "________________________________________________________________________"

    ##### Define the steps to fit VV MC in the mj and mlvj spectra
    def fit_VV(self):
        print "############################# fit_VV ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV_xww", self.jetalgo)
        ### fitting shape as a function of the mlvj region -> signal mass
        #self.fit_mj_single_MC(self.file_VV_mc,"_VV_xww","2Gaus_ErfExp");
        self.fit_mj_single_MC(self.file_VV_mc,"_VV_xww","ExpGaus");
        #self.fit_mj_single_MC(self.file_VV_mc,"_VV_xww","ErfExp");
        self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV_xww","_sb_lo","ExpN", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV_xww","_signal_region",self.MODEL_4_mlvj, 1, 0, 1); 
        print "________________________________________________________________________"

    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):
        print "################################ fit_TTbar #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar_xww", self.jetalgo)# to get the shape of m_lvj
        self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_xww","2Gaus_ErfExp");
	self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar_xww","_sb_lo","ExpN");
        self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar_xww","_signal_region","ExpN",1, 0, 1);
        print "________________________________________________________________________"

    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0_xww", self.jetalgo)# to get the shape of m_lvj
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets01_xww", self.jetalgo)# to get the shape of m_lvj
	
	### Fit in mj depends on the mlvj lower limit -> fitting the turn on at low mass or not
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_xww","User1");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_xww","ExpGaus");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","ErfExp");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","ExpGaus");
        #self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","2Gaus_ErfExp");
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets01_xww","2_2Gaus");
	#sys.exit()
            
        #### Fit the mlvj in sb_lo, signal region using two different model as done in the mj
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0_xww","_sb_lo",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0_xww","_signal_region",self.MODEL_4_mlvj, 0, 0, 1, 1);
        #self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets01_xww","_sb_lo",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets01_xww","_sb_lo",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);
        self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets01_xww","_signal_region",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);       
        print "________________________________________________________________________"

    #### Define the steps to fit signal distribution in the mj and mlvj spectra
    def fit_Signal(self,model_narrow="DoubleCB_v1",model_width="BWDoubleCB"):
        print "############# fit_Signal #################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_signal,"_%s_xww"%(self.signal_sample), self.jetalgo)# to get the shape of m_lvj
        self.fit_mlvj_model_single_MC(self.file_signal,"_%s_xww"%(self.signal_sample),"_signal_region",model_narrow, 0, 0, 0, 0);
        print "________________________________________________________________________"

    #### Interpolate rate and shape for the signal from external files
    def interpolate_Signal(self,model_narrow="DoubleCB_v1",model_width="BWDoubleCB"):
        print "############# interpolate_Signal #################"
        ### Build the dataset
        self.interpolate_signal_MC(self,"_%s_xww"%(self.signal_sample),"_signal_region",model_narrow, 0, 0, 0, 0);
        print "________________________________________________________________________"

    ##### Fit of all the MC in both mj and mlvj : Signal, TTbar, STop, VV and Wjets
    def fit_AllSamples_Mj_and_Mlvj(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.fit_WJets()
        self.fit_VV()
        self.fit_TTbar()
        self.fit_STop()
	#sys.exit()
        print "________________________________________________________________________"

    ##### Analysis with sideband alpha correction 
    def analysis_sideband_correction_method1(self):
        print "##################### Start sideband correction full analysis ##############";
        ### Fit all MC components in both mj and mlvj
        self.fit_AllSamples_Mj_and_Mlvj();
        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm();
        ### fit data in the mlvj low sideband with two different models
        self.fit_mlvj_in_Mj_sideband("_WJets01_xww","_sb_lo",self.MODEL_4_mlvj_alter,1)
        #self.fit_mlvj_in_Mj_sideband("_WJets01_xww","_sb_lo","ExpTail",1)
        #self.fit_mlvj_in_Mj_sideband("_WJets0_xww","_sb_lo","ErfExp",1)
        self.fit_mlvj_in_Mj_sideband("_WJets0_xww","_sb_lo",self.MODEL_4_mlvj,1)
        ### Prepare the workspace and datacards     
        #self.prepare_limit("sideband_correction_method1",1,0,0)
        ### finale plot and check of the workspace
        #self.read_workspace(1)
        
    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    def prepare_limit(self,mode, isTTbarFloating=0, isVVFloating=0, isSTopFloating=0):
        print "####################### prepare_limit for %s method ####################"%(mode);

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_mass_lvj"));

        ### whole number of events from the considered signal sample, WJets, VV, TTbar, STop -> couting
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_WJets0_xww_%s_mlvj"%(self.channel)).clone("rate_WJets_xww_for_counting"))
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_VV_xww_%s_mlvj"%(self.channel)).clone("rate_VV_xww_for_counting"))

        ### number of signal, Wjets, VV, TTbar and STop --> unbin
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_WJets0_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_WJets_xww_for_unbin"));
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_VV_xww_signal_region_%s_mlvj"%(self.channel)).clone("rate_VV_xww_for_unbin"));

        ### Set the error properly -> taking into account lumi, Vtagger and theoretical uncertainty on XS -> for VV, TTbar and STop
        self.workspace4limit_.var("rate_VV_xww_for_unbin").setError(self.workspace4limit_.var("rate_VV_xww_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal() +self.XS_VV_NLO_uncertainty*self.XS_VV_NLO_uncertainty ) );

        ### Get the dataset for data into the signal region
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.data("rdataset_data_xww_signal_region_%s_mlvj"%(self.channel)).Clone("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)))
        ### Take the corrected pdf from the alpha method for the WJets
        if mode=="sideband_correction_method1":
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WJets0_xww_signal_region_%s_after_correct_mlvj"%(self.channel)).clone("WJets_xww_%s_%s"%(self.channel, self.wtagger_label)));

        if isVVFloating :    
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj_Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV"%(self.channel, self.channel, self.wtagger_label)).clone("VV_%s_%s"%(self.channel,self.wtagger_label)))
        else:
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_xww_signal_region_%s_mlvj"%(self.channel)).clone("VV_xww_%s_%s"%(self.channel,self.wtagger_label)))


        ### Fix all the Pdf parameters 
        rrv_x = self.workspace4limit_.var("rrv_mass_lvj");

        self.fix_Pdf(self.workspace4limit_.pdf("VV_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x));
        self.fix_Pdf(self.workspace4limit_.pdf("WJets_xww_%s_%s"%(self.channel,self.wtagger_label)), RooArgSet(rrv_x));
        
        print " ############## Workspace for limit ";
        parameters_workspace = self.workspace4limit_.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()
        
        params_list = [];
        ### main modality for the alpha function method
        if mode=="sideband_correction_method1":

            if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));

                
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)))

                if isTTbarFloating !=0 :
                 self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                 params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));

            if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])));

                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)))


                ### TTbar use exp
                if isTTbarFloating !=0:
                    print "##################### TTbar will float in the limit procedure + final plot ######################";
                    self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
                    params_list.append(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));

                ### VV use ExpTail:
                if isVVFloating !=0:
                  print "##################### VV will float in the limit procedure + final plot ######################";
                  self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_VV);
                  self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_VV);
                  params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                  params_list.append(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)));
 
                ### STop use Exp:
                if isSTopFloating !=0:
                  print "##################### STop will float in the limit procedure + final plot ######################";
                  self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)).setError(self.shape_para_error_STop);
                  params_list.append(self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                                       
                
        #### add signal shape parameters' uncertainty -> increase the uncertainty on the mean and the sigma since we are using a CB or a Double CB or a BWxDB or BWxCB
        if self.workspace4limit_.var("rrv_mean_CB_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)):

           self.workspace4limit_.var( "rrv_mean_shift_scale_lep_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)).setError(self.mean_signal_uncertainty_lep_scale);
           self.workspace4limit_.var( "rrv_mean_shift_scale_jes_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)).setError(self.mean_signal_uncertainty_jet_scale);
           self.workspace4limit_.var( "rrv_sigma_shift_lep_scale_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_lep_scale);
           self.workspace4limit_.var( "rrv_sigma_shift_jes_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_jet_scale);
           self.workspace4limit_.var( "rrv_sigma_shift_res_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_jet_res);

           if self.channel == "mu":
            self.workspace4limit_.var("CMS_sig_p1_scale_m_13TeV").setError(1);
            self.workspace4limit_.var("CMS_sig_p2_scale_m_13TeV").setError(1);
            params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_m_13TeV"));
            params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_m_13TeV"));
           elif self.channel == "el":
            self.workspace4limit_.var("CMS_sig_p1_scale_e_13TeV").setError(1);
            self.workspace4limit_.var("CMS_sig_p2_scale_e_13TeV").setError(1);
            params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_e_13TeV"));
            params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_e_13TeV"));
           elif self.channel == "em":
            self.workspace4limit_.var("CMS_sig_p1_scale_em_13TeV").setError(1);
            self.workspace4limit_.var("CMS_sig_p2_scale_em_13TeV").setError(1);
            params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_em_13TeV"));
            params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_em_13TeV"));
               
           self.workspace4limit_.var("CMS_sig_p1_jes_13TeV").setError(1);
           self.workspace4limit_.var("CMS_sig_p2_jes_13TeV").setError(1);
           self.workspace4limit_.var("CMS_sig_p2_jer_13TeV").setError(1);

           params_list.append(self.workspace4limit_.var("CMS_sig_p1_jes_13TeV"));
           params_list.append(self.workspace4limit_.var("CMS_sig_p2_jes_13TeV"));
           params_list.append(self.workspace4limit_.var("CMS_sig_p2_jer_13TeV"));

        
        ### calculate the shape uncertainty for cut-and-counting
        self.rrv_counting_uncertainty_from_shape_uncertainty = RooRealVar("rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel),"rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.channel),0);
        self.rrv_counting_uncertainty_from_shape_uncertainty.setError( Calc_error("WJets_xww_%s_%s"%(self.channel,self.wtagger_label), "rrv_mass_lvj" ,self.FloatingParams,self.workspace4limit_,"signal_region") );
        self.rrv_counting_uncertainty_from_shape_uncertainty.Print();

        print " param list ",params_list ;
            
        ### Print the datacard for unbin and couting analysis
        self.print_limit_datacard("unbin",params_list);
        self.print_limit_datacard("counting");

        if mode=="sideband_correction_method1":

          if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)) );

                if isTTbarFloating!=0:
                 self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));


          if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sb_lo_from_fitting_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label[0]+self.wtagger_label[1])) );


                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig2"%(self.channel, self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_xww_sim_%s_%s_mlvj_13TeV_eig3"%(self.channel, self.wtagger_label)) );

                if isTTbarFloating!=0:
                    self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));

                if isVVFloating!=0:     
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel, self.wtagger_label)));
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_VV_xww_signal_region_%s_%s_mlvj_13TeV_eig1"%(self.channel, self.wtagger_label)));

                if isSTopFloating!=0:
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_STop_xww_signal_region_%s_%s_mlvj_13TeV_eig0"%(self.channel,self.wtagger_label)));
                    

          if self.workspace4limit_.var("rrv_mean_CB_%s_xww_signal_region_%s_%s"%(self.signal_sample, self.channel, self.wtagger_label)):

             if self.channel == "mu":
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p1_scale_m_13TeV"));
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p2_scale_m_13TeV"));
             elif self.channel == "el":
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p1_scale_e_13TeV"));
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p2_scale_e_13TeV"));
             elif self.channel == "em":
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p1_scale_em_13TeV"));
              self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p2_scale_em_13TeV"));
                 

             self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p1_jes_13TeV"));
             self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p2_jes_13TeV"));
             self.FloatingParams.add(self.workspace4limit_.var("CMS_sig_p2_jer_13TeV"));

        ### Add the floating list to the combiner --> the pdf which are not fixed are floating by default
        getattr(self.workspace4limit_,"import")(self.FloatingParams);

        ### Save the workspace
        self.save_workspace_to_file();

    #### Method used to print the general format of the datacard for both counting and unbinned analysis
    def print_limit_datacard(self, mode, params_list=[]):
        print "############## print_limit_datacard for %s ################"%(mode)
        if not (mode == "unbin" or mode == "counting"):
            print "print_limit_datacard use wrong mode: %s"%(mode);raw_input("ENTER");

        ### open the datacard    
        datacard_out = open(getattr(self,"file_datacard_%s"%(mode)),"w");

        ### start to print inside 
        datacard_out.write( "imax 1" )
        datacard_out.write( "\njmax 4" )
        datacard_out.write( "\nkmax *" )
        datacard_out.write( "\n--------------- ")

        if mode == "unbin":
            fnOnly = ntpath.basename(self.file_rlt_root) ## workspace for limit --> output file for the workspace
            if TString(self.signal_sample).Contains("RS1G_WW"):
             datacard_out.write("\nshapes RSWW_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            elif TString(self.signal_sample).Contains("BulkGraviton"):
             datacard_out.write("\nshapes BulkWW_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            elif TString(self.signal_sample).Contains("Wprime_WZ"):
             datacard_out.write("\nshapes WprimeWZ_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            else:
             datacard_out.write("\nshapes %s_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.signal_sample,self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
                
            datacard_out.write("\nshapes WJets_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes TTbar_xww  CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes STop_xww   CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes VV_xww     CMS_xww_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write("\nshapes data_obs   CMS_xww_%s1J%s  %s %s:$PROCESS_xww_%s_%s"%(self.channel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.channel, self.wtagger_label));
            datacard_out.write( "\n--------------- ")
            
        datacard_out.write( "\nbin CMS_xww_%s1J%s "%(self.channel,self.wtagger_label));    
        if mode == "unbin":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)).sumEntries()) )
        if mode == "counting":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting_xww").getVal()) )
            
        datacard_out.write( "\n------------------------------" );

        datacard_out.write( "\nbin CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s CMS_xww_%s1J%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label,self.channel,self.wtagger_label));

        if TString(self.signal_sample).Contains("RS1G_WW"):
         datacard_out.write( "\nprocess RSWW_xww WJets_xww TTbar_xww STop_xww VV_xww"); ## just one signal sample
        elif TString(self.signal_sample).Contains("BulkGraviton"):
         datacard_out.write( "\nprocess BulkWW_xww WJets_xww TTbar_xww STop_xww VV_xww"); ## just one signal sample
        elif TString(self.signal_sample).Contains("Wprime_WZ"):
         datacard_out.write( "\nprocess WprimeWZ_xww WJets_xww TTbar_xww STop_xww VV_xww"); ## just one signal sample
        else:
         datacard_out.write( "\nprocess %s_xww WJets_xww TTbar_xww STop_xww VV_xww "%(self.signal_sample)); ## just one signal sample

        datacard_out.write( "\nprocess -1 1 2 3 4" );

        ### rates for the different process
        datacard_out.write( "\n-------------------------------- " )

        ### luminosity nouisance
        datacard_out.write( "\nlumi_13TeV lnN %0.3f - %0.3f %0.3f %0.3f"%(1.+self.lumi_uncertainty, 1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty) )

        ### STop XS  nouisance in boosted regime
        datacard_out.write( "\nCMS_xww_XS_STop_13TeV lnN - - - %0.3f -"%(1+self.XS_STop_NLO_uncertainty) )

        ### VV XS  nouisance in boosted regime
        datacard_out.write( "\nCMS_xww_XS_VV_13TeV lnN - - - - %0.3f"%(1+self.XS_VV_NLO_uncertainty) )

        ### WJets Normalization from data fit -> data driven
        if self.number_WJets_insideband >0:
            datacard_out.write( "\nCMS_xww_WJ_norm_13TeV gmN %0.3f %0.3f - - -"%(self.number_WJets_insideband, getattr(self, "datadriven_alpha_WJets_xww_%s"%(mode)) ) )
        else:
            datacard_out.write( "\nCMS_xww_WJ_norm_%s_%s_13TeV lnN - %0.3f - - -"%(self.channel, self.wtagger_label, 1+ self.workspace4limit_.var("rate_WJets_xww_for_unbin").getError()/self.workspace4limit_.var("rate_WJets_xww_for_unbin").getVal() ) );

        ### Top normalization due to SF in the ttbar CR
        datacard_out.write( "\nCMS_xww_Top_norm_%s_%s_13TeV lnN - %0.3f/%0.3f %0.3f/%0.3f %0.3f/%0.3f -"%(self.channel, self.wtagger_label, 1-0.6*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+0.6*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1-self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() ) );

        ### V-Tagger SF nouisance
        datacard_out.write( "\nCMS_eff_vtag_tau21_sf_13TeV lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError()));
            

        ### btag scale factor on the MC background
        datacard_out.write( "\n#CMS_eff_vtag_model_13TeV lnN %0.3f - - - %0.3f"%(1+self.eff_vtag_model,1+self.eff_vtag_model) );

        ### jet Mass effect only if available -> shapes changing due to the jet mass uncertainty (JEC for CA8/AK7) -> affects also WJets
        if ( self.workspace4fit_.var("rrv_number_WJets0_xww_massup_in_mj_signal_region_from_fitting_%s"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_WJets0_xww_massdn_in_mj_signal_region_from_fitting_%s"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_xww_massdown_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_xww_massdn_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massup_%s_mj"%(self.channel)) and
             self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_xww_massdn_%s_mj"%(self.channel))) :
                                                                                                                                                         
         datacard_out.write( "\nJetMass_%s lnN - %0.3f %0.3f %0.3f %0.3f"%(self.channel, 1+self.WJets_normlization_uncertainty_from_jet_mass, 1+self.TTbar_normlization_uncertainty_from_jet_mass, 1+self.STop_normlization_uncertainty_from_jet_mass, 1+self.VV_normlization_uncertainty_from_jet_mass ) )

        if self.channel == "mu":
            self.channel_short = "m"
        elif self.channel =="el":
            self.channel_short = "e"
        elif self.channel =="em":
            self.channel_short = "em"
            
        ### trigger efficiency
        datacard_out.write( "\nCMS_xww_trigger_%s_13TeV lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty ) );

        ### Lepton SF
        datacard_out.write( "\nCMS_eff_%s_13TeV lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty ) );

        ############ Evaluated just for signal, in principle also on all the backgrounds with the same topology

        ### Lepton Energy scale
        datacard_out.write( "\nCMS_scale_%s_13TeV lnN %0.3f - - - -"%(self.channel_short,1+self.signal_lepton_energy_scale_uncertainty));

        ### Lepton Energy Resolution
        datacard_out.write( "\nCMS_res_%s_13TeV lnN %0.3f - - - -"%(self.channel_short,1+self.signal_lepton_energy_res_uncertainty));

        ### CA8 jet energy scale
        #datacard_out.write( "\nCMS_scale_j_13TeV lnN %0.3f - - - -"%(1+self.signal_jet_energy_scale_uncertainty));
        datacard_out.write( "\nCMS_scale_j_13TeV lnN %0.3f/%0.3f - - - -"%(1+self.signal_jet_energy_scale_uncertainty_up,1+self.signal_jet_energy_scale_uncertainty_down));

        ### jet mass scale
        datacard_out.write( "\nCMS_mass_scale_j_13TeV lnN %0.3f/%0.3f - - - -"%(1+self.signal_jet_mass_scale_uncertainty_up,1+self.signal_jet_mass_scale_uncertainty_down));

        ### jet mass resolution
        datacard_out.write( "\nCMS_mass_res_j_13TeV lnN %0.3f/%0.3f - - - -"%(1+self.signal_jet_mass_res_uncertainty_up,1+self.signal_jet_mass_res_uncertainty_down));

        ### CA8 jet energy resolution
        datacard_out.write( "\nCMS_res_j_13TeV lnN %0.3f - - - -"%(1+self.signal_jet_energy_res_uncertainty));


        ### btag on the signal
        datacard_out.write( "\nCMS_xww_btag_eff_13TeV lnN %0.3f - - - -"%(1+self.signal_btag_uncertainty));


        ### print shapes parameter to be taken int account
        if mode == "unbin":
            for ipar in params_list:
              print "Name %s",ipar.GetName();
              if TString(ipar.GetName()).Contains("Deco_TTbar_xww_signal_region"):
               datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )
              else:
               datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )
        if mode == "counting":
           datacard_out.write( "\nShape_%s_%s_13TeV lnN - - %0.3f - - -"%(self.channel, self.wtagger_label, 1+self.rrv_counting_uncertainty_from_shape_uncertainty.getError()))

    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root);
        self.file_out.close()

    #### Read the final workspace and produce the latest plots 
    def read_workspace(self, logy=0):
        print "--------------------- read_workspace -------------------------";

        ### Taket the workspace for limits  
        file = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;
        workspace.Print()

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------";
        parameters_workspace = workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next()
        print "---------------------------------------------";

        workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label)).Print()

        print "----------- Pdf in the Workspace -------------";
        pdfs_workspace = workspace.allPdfs();
        par = pdfs_workspace.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            param.Print();
            param = par.Next()
        print "----------------------------------------------";

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_xww_%s_%s"%(self.channel,self.wtagger_label));
            
        model_pdf_WJets  = workspace.pdf("WJets_xww_%s_%s"%(self.channel,self.wtagger_label));
        model_pdf_VV     = workspace.pdf("VV_xww_%s_%s"%(self.channel,self.wtagger_label));
        #model_pdf_TTbar  = workspace.pdf("TTbar_xww_%s_%s"%(self.channel,self.wtagger_label));
        #model_pdf_STop   = workspace.pdf("STop_xww_%s_%s"%(self.channel,self.wtagger_label));

        model_pdf_WJets.Print();
        model_pdf_VV.Print();
        #model_pdf_TTbar.Print();
        #model_pdf_STop.Print();

        rrv_number_WJets  = workspace.var("rate_WJets_xww_for_unbin");
        rrv_number_VV     = workspace.var("rate_VV_xww_for_unbin");
        #rrv_number_TTbar  = workspace.var("rate_TTbar_xww_for_unbin");
        #rrv_number_STop   = workspace.var("rate_STop_xww_for_unbin");

        rrv_number_WJets.Print();
        rrv_number_VV.Print();
        #rrv_number_TTbar.Print();
        #rrv_number_STop.Print();

        #### Prepare the final plot starting from total background 
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww","rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()* rrv_number_WJets.getError()+
                rrv_number_VV.getError()* rrv_number_VV.getError()));

        #### Total pdf 
        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww","model_Total_background_MC_xww",RooArgList(model_pdf_WJets,model_pdf_VV),RooArgList(rrv_number_WJets,rrv_number_VV));

        if data_obs.sumEntries() != 0:
         #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
         scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries()
        else:
	 scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()   

        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        mplotP = rrv_x.frame(RooFit.Title("check_workspaceP"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_xww_%s_%s"%(self.channel,self.wtagger_label)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());


        #solid line
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s"%(self.channel,self.wtagger_label,self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_xww_%s_%s"%(self.channel,self.wtagger_label)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());
 
        #### plot the observed data using poissonian error bar
        self.getData_PoissonInterval(data_obs,mplot);
        
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible());

        mplot_pull=self.get_pull(rrv_x,mplot);

        datahist = data_obs.binnedClone( data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone" )
        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band2(model_Total_background_MC, rrv_x.GetName(), rrv_number_Total_background_MC,self.FloatingParams,workspace ,mplot,mplotP,datahist,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.leg = self.legend4Plot(mplot,0,1,-0.01,-0.05,0.11,0.);
        #self.leg.SetTextSize(0.036);
        mplot.addObject(self.leg);
	pt1 = ROOT.TPaveText(0.6180905,0.4355644,0.8291457,0.507992,"NDC")
	pt1.SetTextFont(42)
	pt1.SetTextSize(0.05)
	#pt.SetTextAlign(12)
	pt1.SetFillColor(0)
	pt1.SetFillStyle(0)
	pt1.SetBorderSize(0)
	text = pt1.AddText("")
	if options.category.find('Z') != -1: text = pt1.AddText("WZ category")
	elif options.category.find('W') != -1: text = pt1.AddText("WW category")
	text.SetTextFont(62)
#	mplot.addObject(pt1)
	            
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
            
        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal());
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-self.nPar_float_in_fitTo;
        print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );

        parameters_list = RooArgList();
        self.draw_canvas_with_pull2( rrv_x,datahist,mplot, mplotP, mplot_pull,ndof,parameters_list,"%s/m_lvj_fitting/"%(self.plotsDir),"check_workspace_for_limit","",0,1);
        

### funtion to run the complete alpha analysis
def pre_limit_sb_correction(method, channel, signal_sample="BulkG_c0p2_M1000", jetalgo="Mjsoftdrop",in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700,
                            in_mj_min=30, in_mj_max=140, in_mlvj_min=400, in_mlvj_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1",interpolate=False): 

    print "#################### pre_limit_sb_correction: channel %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel,signal_sample,in_mlvj_signal_region_min,in_mlvj_signal_region_max,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);
    print method    
    boostedW_fitter=doFit_wj_and_wlvj(channel, signal_sample, jetalgo,in_mlvj_signal_region_min,in_mlvj_signal_region_max,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter,interpolate);
    getattr(boostedW_fitter,"analysis_sideband_correction_%s"%(method) )();
                                                

#### Main Code
if __name__ == '__main__':

    clock = TStopwatch()
    clock.Start()

    channel=options.channel;
    mass=options.mass;
    sample = options.sample+str(int(mass))
    
    lomass = 170;
    himass = 3000; 
    
    os.system('echo "Deleting plot directories...";rm -r plots_em_HP cards_em_HP')
    #pre_limit_sb_correction("method1",channel,sample,options.jetalgo, 600,3000,40,150, 600,3000,"Exp","ExpTail",options.interpolate) 
    pre_limit_sb_correction("method1",channel,sample,options.jetalgo, 600,3000,40,150, 600,3000,"ExpTail","Exp",options.interpolate) 

    print "\n\n","_"*30,"\n\n\t Fit results of mj","\n\n"
    for i in fitresultsmj:
	print "\n\n=====\t",i,"\t======"
        i.Print()
    print "\n\n","_"*30,"\n\n\t Fit results of mlvj","\n\n"
    for i in fitresultsmlvj:
	print "\n\n=====\t",i,"\t======"
        i.Print()
    print "\n\n","_"*30,"\n\n\t Decorrelated results of mlvj","\n\n"
    for i in fit_DecoResults:
	print "\n\n=====\t",i,"\t======"
    	i.Print()
    print "\n\n","_"*30,"\n\n\t Final results","\n\n"
    for i in fitresultsfinal:
	print "\n\n=====\t",i,"\t======"
        i.Print()

    clock.Stop()
    print 'Tree loop profiling stats:'
    print 'Real Time used:', clock.RealTime()/60,"minutes"
    print 'CPU Time used:', clock.CpuTime()/60,"minutes"
