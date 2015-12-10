import os,commands
import sys
import random
from array import array
from optparse import OptionParser
from optparse import OptionGroup
import ROOT
from ROOT import gROOT, gStyle, gSystem, TLatex, TH1D, TString, TPaveText, TGaxis
import subprocess
from subprocess import Popen
from optparse import OptionParser
import CMS_lumi, tdrstyle
import time

#$1: outfile
#$2: jobid
#$3: cmsswdir
#$4: datacard
#$5: singlepoint
#$6: seed
#$7: mass
#$8: suffix
#$9: seoutdir
#$10: toys
#$11: iterations
#$12: subid

def get_canvas(text):

   tdrstyle.setTDRStyle()
   CMS_lumi.lumi_8TeV = "19.7 fb^{-1}," + text
   CMS_lumi.writeExtraText = 1
   CMS_lumi.extraText = "Preliminary"

   iPos = 0
   if( iPos==0 ): CMS_lumi.relPosX = 0.12

   H_ref = 600 
   W_ref = 800 
   W = W_ref
   H  = H_ref

   T = 0.08*H_ref
   B = 0.12*H_ref 
   L = 0.12*W_ref
   R = 0.04*W_ref

   canvas = ROOT.TCanvas("c2","c2",50,50,W,H)
   canvas.SetFillColor(0)
   canvas.SetBorderMode(0)
   canvas.SetFrameFillStyle(0)
   canvas.SetFrameBorderMode(0)
   canvas.SetLeftMargin( L/W )
   canvas.SetRightMargin( R/W )
   canvas.SetTopMargin( T/H )
   canvas.SetBottomMargin( B/H )
   canvas.SetGrid()
   canvas.SetLogy()
   
   return canvas

def getAsymLimits(file):
    
    
    f = ROOT.TFile(file)
    t = f.Get("limit")
    entries = t.GetEntries()
    
    lims = [0,0,0,0,0,0]
    
    for i in range(entries):
        
        t.GetEntry(i)
        t_quantileExpected = t.quantileExpected
        t_limit = t.limit
        
        #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected
        
        if t_quantileExpected == -1.: lims[0] = t_limit
        elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit
        elif t_quantileExpected >= 0.15  and t_quantileExpected <= 0.17:  lims[2] = t_limit
        elif t_quantileExpected == 0.5: lims[3] = t_limit
        elif t_quantileExpected >= 0.83  and t_quantileExpected <= 0.85:  lims[4] = t_limit
        elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit
        else: print "Unknown quantile!"
    
    return lims

def doLimitPlot( ):

    nPoints = len(mass)
    
    xbins        = array('d', [])
    xbins_env    = array('d', [])
    ybins_exp    = array('d', [])
    ybins_obs    = array('d', [])
    ybins_1s     = array('d', [])
    ybins_2s     = array('d', [])
    ybins_xs_HVT = array('d', [])
    ybins_xs_LH  = array('d', [])
    
    br_lnubb = 0.1057*0.577
    
    if channel_ =="el" :
        text = "W#rightarrow e#nu"
        br_lnubb = 0.1075*0.577
    if channel_ =="mu" :
        text = "W#rightarrow #mu#nu"
        br_lnubb = 0.1057*0.577
    if channel_ =="em" :   
        text = "W#rightarrow l#nu"
        br_lnubb = 0.108*0.577
    
    
    for i in range(len(mass)):
        curFile = "m%i-fullCLS/higgsCombine_%s_.HybridNew.mH%i.root"%(mass[i],channel_,mass[i])
        xsHVT = xsDict2[mass[i]]
        xsLH  = xsDictLH[mass[i]]
        curAsymLimits = getAsymLimits(curFile)
        xbins.append( mass[i] )
        xbins_env.append( mass[i] )
        ybins_exp.append( curAsymLimits[3]*xsHVT )
        ybins_obs.append( curAsymLimits[0]*xsHVT )
        ybins_2s.append( curAsymLimits[1]*xsHVT )
        ybins_1s.append( curAsymLimits[2]*xsHVT )
	ybins_xs_HVT.append(xsHVT)
        ybins_xs_LH.append(xsLH)
	print "mass %i exp %.6f obs %.6f 2s %.6f 1s %.6f" %(mass[i],curAsymLimits[3],curAsymLimits[0],curAsymLimits[1],curAsymLimits[2])
    
    for i in range( len(mass)-1, -1, -1 ):
        curFile = "m%i-fullCLS/higgsCombine_%s_.HybridNew.mH%i.root"%(mass[i],channel_,mass[i])
        xsHVT = xsDict2[mass[i]]
        xsLH  = xsDictLH[mass[i]]
        curAsymLimits = getAsymLimits(curFile)
        xbins_env.append( mass[i] )
        ybins_2s.append( curAsymLimits[5]*xsHVT )
        ybins_1s.append( curAsymLimits[4]*xsHVT )
	print "mass %i exp %.6f obs %.6f 2s %.6f 1s %.6f" %(mass[i],curAsymLimits[3],curAsymLimits[0],curAsymLimits[5],curAsymLimits[4])
    
    canv = get_canvas(text)
    canv.cd()
        
    curGraph_exp    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp)
    curGraph_obs    = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs)
    curGraph_xs_HVT = ROOT.TGraph(nPoints,xbins,ybins_xs_HVT)
    curGraph_xs_LH  = ROOT.TGraph(nPoints,xbins,ybins_xs_LH)
    curGraph_1s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s)
    curGraph_2s     = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s)
    
    curGraph_obs.SetMarkerStyle(20)
    curGraph_obs.SetLineWidth(3)
    curGraph_obs.SetLineStyle(1)
    curGraph_obs.SetMarkerSize(1.6)
    curGraph_exp.SetMarkerSize(1.3)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_exp.SetLineStyle(2)
    curGraph_exp.SetLineWidth(3)
    curGraph_exp.SetMarkerSize(2)
    curGraph_exp.SetMarkerStyle(24)
    curGraph_exp.SetMarkerColor(ROOT.kBlack)

    curGraph_xs_HVT.SetLineStyle(ROOT.kSolid)
    curGraph_xs_HVT.SetFillStyle(3344)
    curGraph_xs_HVT.SetLineWidth(2)
    curGraph_xs_HVT.SetMarkerSize(2)
    curGraph_xs_HVT.SetLineColor(ROOT.kRed)

    curGraph_xs_LH.SetLineStyle(ROOT.kSolid)
    curGraph_xs_LH.SetFillStyle(3344)
    curGraph_xs_LH.SetLineWidth(2)
    curGraph_xs_LH.SetMarkerSize(2)
    curGraph_xs_LH.SetLineColor(ROOT.kBlue)
    
    curGraph_1s.SetFillColor(ROOT.kGreen)
    curGraph_1s.SetFillStyle(1001)
    curGraph_1s.SetLineStyle(ROOT.kDashed)
    curGraph_1s.SetLineWidth(3)

    curGraph_2s.SetFillColor(ROOT.kYellow)
    curGraph_2s.SetFillStyle(1001)
    curGraph_2s.SetLineStyle(ROOT.kDashed)
    curGraph_2s.SetLineWidth(3)

    hrl_SM = canv.DrawFrame(750,10e-8, 3050, 10)
    hrl_SM.GetYaxis().SetTitle("#sigma_{95%} #times BR(W' #rightarrow WH)(pb)")
    hrl_SM.GetYaxis().CenterTitle()
    hrl_SM.GetYaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetTitleSize(0.06)
    hrl_SM.GetXaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetLabelSize(0.045)
    hrl_SM.GetYaxis().SetTitleOffset(0.9)
    hrl_SM.GetXaxis().SetTitleOffset(0.9)
    hrl_SM.GetXaxis().SetTitle("M_{W'} (GeV)")
    hrl_SM.GetXaxis().CenterTitle()
    hrl_SM.SetMinimum(0.001)
    hrl_SM.SetMaximum(10)
    hrl_SM.GetYaxis().SetNdivisions(505)
        
    curGraph_2s.Draw("F")
    curGraph_1s.Draw("Fsame")
    curGraph_exp.Draw("Lsame")
    curGraph_obs.Draw("LPsame")
    curGraph_xs_HVT.Draw("Csame")
    curGraph_xs_LH.Draw("Csame")
       
    leg2 = ROOT.TLegend(0.52,0.68,0.93,0.88)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.03)

    leg2.AddEntry(curGraph_1s,"Asympt. CL_{S}  Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s,"Asympt. CL_{S}  Expected #pm 2#sigma","LF")
    leg2.AddEntry(curGraph_xs_HVT,"HVT B(gv=3)","L")
    leg2.AddEntry(curGraph_xs_LH,"Little Higgs","L")
     
    leg2.Draw()
             
    canv.Update()   
    canv.cd()
    CMS_lumi.CMS_lumi(canv, 2, 0)	   
    canv.cd()
    canv.Update()
    canv.RedrawAxis()
    canv.RedrawAxis("g")
    frame = canv.GetFrame()
    frame.Draw()   
    canv.cd()
    canv.Update()    
    
    time.sleep(1000)
    
    os.system('mkdir LimitResult')
    os.system('mkdir LimitResult/Limit_ExpTail')
    
    canv.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.png"%(channel_))
    canv.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.pdf"%(channel_))
    canv.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.root"%(channel_))
    canv.SaveAs("./LimitResult/Limit_ExpTail/Lim%s.C"%(channel_))
       
parser = OptionParser()
parser.add_option("--outdir", dest="outdir", action="store",
                  help="SE output subdir (default = None)", metavar="OUTDIR", type="string", 
                  default="")
parser.add_option("--cmsswdir", dest="cmsswdir", action="store",
                  help="cmssw dir", metavar="CMSSWDIR", type="string", 
                  default="")
parser.add_option("--datadir", dest="datadir", action="store",
                  help="datacard dir", metavar="DATADIR", type="string", 
                  default="")
parser.add_option("--channel", dest="channel", action="store",
                  help="channel (el,mu,combo)", metavar="CHANNEL", type="string", 
                  default="")
parser.add_option("--toys", dest="toys", action="store",
                  help="number of toys", metavar="TOYS", type="int", 
                  default=1)
parser.add_option("--i", dest="iterations", action="store",
                  help="number of iterations", metavar="TOYS", type="int", 
                  default=3)
parser.add_option("--merge", dest="mergefiles", action="store_true",
                  help="merge job outputs (default = False)", metavar="MERGEFILES", 
                  default=False)
parser.add_option("--runtoys", dest="runtoys", action="store_true",
                  help="run toys (default = False)", metavar="MERGEFILES", 
                  default=False)	
parser.add_option("--runfullcls", dest="runfullcls", action="store_true",
                  help="run full cls (default = False)", metavar="MERGEFILES", 
                  default=False)
parser.add_option("--plotlimits", dest="plotlimits", action="store_true",
                  help="plot limits (default = False)", metavar="PLOTLIMITS", 
                  default=False)
		  		  		  	  		  		  		  
(options, args) = parser.parse_args()

cmsswdir_ = options.cmsswdir
outdir_ = options.outdir
datadir_ = options.datadir
channel_ = options.channel
toys_ = options.toys
it_ = options.iterations

points = []		  
for p in range(1,10):
   points+=[float(p/10.)]
   points+=[float(p/10.+0.05)]
   points+=[float(p/1.)]
   points+=[float(p/1.+0.5)]
   points+=[float(p*10.)]
   points+=[float(p*10.+5.)]

mass = [1900]
#mass = [800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500]

xsDict2 = {800:0.337368865,900:0.248172934,
          1000:0.170538599,1100:0.115925052,1200:0.080506243,1300:0.055909013,
          1400:0.038827022,1500:0.025069330,1600:0.018725724,1700:0.013004416,
          1800:0.009031152,1900:0.006271846,2000:0.004251021,2100:0.003024823,
	  2200:0.002100643,2300:0.001458829, 2400:0.001013110, 2500:0.000731277,
          2600:0.000488609,2700:0.000339323,2800:0.000235649,2900:0.000163651,3000:0.000120961}

xsDictLH = {800:0.5087,900:0.3026,
            1000:0.1866,1100:0.1181,1200:0.07653,1300:0.05061,
            1400:0.03386,1500:0.0229,1600:0.01562,1700:0.01075,
            1800:0.007426,1900:0.005172,2000:0.003612,2100:0.002527,
            2200:0.001763,2300:0.001237,2400:0.000867,2500:0.0006071,
            2600:0.0004247,2700:0.0002964,2800:0.0002069,2900:0.0001438,3000:0.00009972} 
	  
if options.runtoys:

   for m in mass:

      for p in range(len(points)):

         point = points[p]
	 for i in range(10):
            seed = int(random.random()*pow(10,7))
            outfile = "higgsCombine_%s_.HybridNew.mH%i.%i.root" %(channel_,m,seed)
            datacard = "%s/cards_EXO_%s_ALLP_g1/whlvj_MWp_%i_bb_%s_ALLP_unbin.txt" %(datadir_,channel_,m,channel_)
            cmd = "qsub -q all.q submitJobsOnT3batch.sh %s %s %s %s %f %i %i _%s_ %s %i %i %i" %(outfile,p,cmsswdir_,datacard,point,seed,m,channel_,outdir_,toys_,it_,i)
            print cmd
            os.system(cmd)

if options.mergefiles:

   user = os.popen('whoami').read()
   user = user.split('\n')[0]
   
   for m in mass: 
        
      status,ls_la = commands.getstatusoutput( 'ls -l m%i-fullCLS'%(m) )  													      
      if status:																				      
         os.system('mkdir m%i-fullCLS'%(m) )
      #else:
         #os.system('rm -rf m%i-fullCLS'%(m) )
	 #os.system('mkdir m%i-fullCLS'%(m) )	 
      
      for p in xrange(45,len(points)):
      
         for i in range(10):
            outdir = "jobtmp"
            cmd = "uberftp t3se01.psi.ch 'ls /pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i'" %(user,outdir_,m,p,i)
            status,ls_la = commands.getstatusoutput( cmd )
            filename = ''
            list_ = ls_la.split(os.linesep)

            for a in list_:
              b = a.split(" ")
              if b[-1:][0].find("root") != -1:
            	 filename = b[len(b)-1].split('\r')[0]
       
            cmd = "lcg-cp -bD srmv2 srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i/%s m%i-fullCLS/%s" %(user,outdir_,m,p,i,filename,m,filename) 
            #print cmd
	    #os.system(cmd)
	    cmd = "srmrm srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i/%s" %(user,outdir_,m,p,i,filename)
	    print cmd
	    os.system(cmd)
	    cmd = "srmrmdir srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/%s/%s/m%i-%i-%i/" %(user,outdir_,m,p,i)
            print cmd
	    os.system(cmd)
	    
         hadd = "hadd -f m%i-fullCLS/grid_mX%i_p%i_%s.root m%i-fullCLS/higgsCombine*" %(m,m,p,channel_,m)
         #print hadd
         #os.system(hadd)
	 rm = "rm m%i-fullCLS/higgsCombine*" %(m)
         #print rm
	 #os.system(rm)
	    
      hadd = "hadd -f m%i-fullCLS/grid_mX%i_%s.root m%i-fullCLS/grid_mX*" %(m,m,channel_,m)
      #print hadd
      #os.system(hadd)

if options.runfullcls:

   for m in mass:
   
      datacard = "%s/cards_EXO_%s_ALLP_g1/whlvj_MWp_%i_bb_%s_ALLP_unbin.txt" %(datadir_,channel_,m,channel_)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.5" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.16" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.84" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)                  
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.025" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)
      cmd = "combine %s -M HybridNew --testStat LHC --grid m%i-fullCLS/grid_mX%i_%s.root -m %i -n _%s_ --expectedFromGrid 0.975" %(datacard,m,m,channel_,m,channel_)
      print cmd
      os.system(cmd)
      
if options.plotlimits:

   doLimitPlot()      
      
#python runFullCLS.py --channel el --cmsswdir /shome/jngadiub/CMSSW_6_1_1 --datadir /shome/jngadiub/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit --outdir jobtmp --runtoys --toys 1 --i 3
#python runFullCLS.py --channel el --outdir jobtmp --merge
#python runFullCLS.py --channel el --datadir /shome/jngadiub/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit --runfullcls
#python runFullCLS.py --channel el --plotlimits
