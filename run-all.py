import os,commands
import sys
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")
parser.add_option('-s', action="store",type="string",dest="sample",default="BulkG_WW")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="Mjpruned")
(options, args) = parser.parse_args()

#masses = [800,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500]
masses = [3000,3500,4000,4500]
for m in masses:
   cmd = "python g1_exo_doFit_class.py -b -c %s --mass %i --category %s --sample %s_lvjj --jetalgo %s > log/%s_M%i_%s_%s.log" %(options.channel,m,options.category,options.sample,options.jetalgo,options.sample,m,options.channel,options.category)
   print cmd
   os.system(cmd)

#python run-all.py --channel mu -s Wprime_WZ --jetalgo Mjsoftdrop --category HP
#python run-all.py -c mu -s BulkG_WW --category HPW
