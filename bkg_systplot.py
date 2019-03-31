import ROOT as r
import sys
import Utils

fileIn = r.TFile("WV_bkg_estimation_4Bins_50GeVLepCut.root","read");
#fileIn = r.TFile("Check_bkg_bug.root","read");
#fileIn = r.TFile("Check_bkg_bug_2500.root","read");
#r.gStyle.SetErrorX(0.0001);

c1 = r.TCanvas("c1");
c1, pad1, pad2 = Utils.createCanvasPads();
legend = r.TLegend(0.5,0.7,0.9,0.9);

NominalHist =fileIn.Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_Down");
for i in range(1,5):
   NominalHist.SetBinError(i,0.001)
NominalHist.SetMarkerColor(1);
NominalHist.SetMarkerSize(1);
NominalHist.SetLineColor(1);
NominalHist.SetFillStyle(0);
NominalHist.SetTitle("")
NominalHist.GetXaxis().SetTitle("M_{WW} (GeV)")

lUp = fileIn.Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_ParVar_Up");
lUp.SetMarkerColor(2);
lUp.SetMarkerSize(1);
lUp.SetLineColor(2);
lUp.SetFillStyle(0);
lDown = fileIn.Get("WjetFitSyst_SignalRegion_Corr_Hist_From_Data_4bins_AlternateShape_ParVar_Down");
lDown.SetMarkerColor(3);
lDown.SetMarkerSize(1);
lDown.SetLineColor(4);
lDown.SetFillStyle(0);

NominalHist.SetTitle("")
c1.cd()
pad1.cd()
NominalHist.Draw();
lUp.Draw("same hist")
lDown.Draw("same hist")

legend.AddEntry(NominalHist,"nominial","lpe");
legend.AddEntry(lUp,"+1#sigma shift","lp");
legend.AddEntry(lDown,"-1#sigma shift","lp");
legend.Draw()

pad2.cd()
hratio1 = Utils.createRatio(NominalHist,  lUp, "M_{WW} (GeV)", 2)
hratio2 = Utils.createRatio(NominalHist,  lDown, "M_{WW} (GeV)", 4)

for i in range(1,5):
  hratio1.SetBinError(i,0.001)
  hratio2.SetBinError(i,0.001)

hratio1.Draw("")
hratio2.Draw("same ")
l = r.TLine(600,1.0,2500,1.0);
l.SetLineColor(1);
l.SetLineStyle(3);
l.SetLineWidth(1);
l.Draw();


c1.SaveAs("WjetFitSyst_AlternateFun_Par0_Wjets.png")
c1.SaveAs("WjetFitSyst_AlternateFun_Par0_Wjets.pdf")
c1.SetLogy(1)
c1.SaveAs("WjetFitSyst_AlternateFun_Par0_Wjets_log.png")
c1.SaveAs("WjetFitSyst_AlternateFun_Par0_Wjets_log.pdf")
c1.SetLogy(0)
c1.Clear()
