#!/bin/python
'''Giles Strong (giles.strong@outlook.com)'''

from __future__ import division
import os
import os.path
import optparse
from data import *

outDir = '/afs/cern.ch/work/g/gstrong/private/hhYRSkims/'
softDir = '$HOME/Simple_Delphes_Event_Selection/src'

def makeJOFile(inputFile, uid, opts):
    outputFile = outDir + opts.sample + "_" + str(uid)
    cmd = "./delphes_event_selection "
    cmd += "-i " + inputFile
    cmd += " -o " + outputFile
    cmd += " -d " + str(opts.debug)

    joName = opts.sample + "_" + str(uid) + ".sh"
    joFile = open(joName, "w")
    joFile.write("#!/bin/sh\n")
    joFile.write("echo Beginning\ job\n")
    joFile.write("export HOME=/afs/cern.ch/user/g/gstrong/\n")
    joFile.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_87/x86_64-slc6-gcc49-opt/setup.sh\n")
    joFile.write("source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.06.00-a82d7/x86_64-slc6-gcc49-opt/bin/thisroot.sh\n")
    joFile.write("export LD_LIBRARY_PATH=$HOME/programs/delphes/:$LD_LIBRARY_PATH\n")
    joFile.write("export X509_USER_PROXY=$HOME/x509up_u5020023\n")
    joFile.write("cd " + softDir + "\n")
    joFile.write("echo Paths\ set\n")
    joFile.write(cmd + "\n")
    joFile.close()

    os.system("chmod 744 " + joName)

def make_sub_file(opts):
    subName = "condor.sub"
    subFile = open(subName, "w")
    subFile.write("executable = $(filename)\n")
    subFile.write("output = " + outDir + "/" + opts.sample + "_$(ProcId).out\n")
    subFile.write("error = " + outDir + "/" + opts.sample + "_$(ProcId).err\n")
    subFile.write("log = " + outDir + "/" + opts.sample + "_$(ProcId).log\n") 
    subFile.write('requirements = (OpSysAndVer =?= "SLCern6")\n')
    subFile.write('+JobFlavour = "' + opts.queue + '"\n')
    subFile.write("queue filename matching (" + opts.sample + "_*.sh)")
    subFile.close()

    sub = "condor_submit " + subName
    print "Submitting: " + sub
    os.system("chmod 744 " + subName)
    os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-s", "--sample", dest = "sample", action = "store", help = "Sample to analyse")
    parser.add_option("-d", "--debug", dest = "debug", action = "store", default = 0, help = "Run in debug mode. {0,1}, default: 0")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = "workday", help = "Queue to run jobs. Default: normal")
    parser.add_option("-f", "--first", dest = "first", action = "store", default = 0, help = "First job to run. Default: 0")
    opts, args = parser.parse_args()

    os.system("voms-proxy-init -voms cms -valid 72:00")
    os.system("cp /tmp/x509up_u61049 /afs/cern.ch/user/g/gstrong/x509up_u5020023")
    
    if opts.sample == "signal":
        files = signalFiles
        loc = signalLoc

    elif opts.sample == "ttbar":
        loc = ttbarLoc
        filesttbar = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
        
    elif opts.sample == "ttbar_DiLepton":
        files = ttbar_DiLeptonFiles
        loc = ttbar_DiLeptonLoc

    elif opts.sample == "htautau":
        files = htautauFiles
        loc = htautauLoc

    elif opts.sample == "hbb":
        loc = hbbLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    elif opts.sample == "dy70":
        loc = dy70Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy100":
        loc = dy100Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy200":
        loc = dy200Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy400":
        loc = dy400Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy600":
        loc = dy600Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy800":
        loc = dy800Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy1200":
        loc = dy1200Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "dy2500":
        loc = dy2500Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    elif opts.sample == "qcd0":
        files = qcd0Files
        loc = qcd0Loc
    elif opts.sample == "qcd50":
        files = qcd50Files
        loc = qcd50Loc
    elif opts.sample == "qcd80":
        files = qcd80Files
        loc = qcd80Loc
    elif opts.sample == "qcd120":
        files = qcd120Files
        loc = qcd120Loc
    elif opts.sample == "qcd170":
        files = qcd170Files
        loc = qcd170Loc
    elif opts.sample == "qcd300":
        files = qcd300Files
        loc = qcd300Loc
    elif opts.sample == "qcd470":
        loc = qcd470Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == "qcd600":
        files = qcd600Files
        loc = qcd600Loc
    elif opts.sample == "qcd800":
        files = qcd800Files
        loc = qcd800Loc
    elif opts.sample == "qcd1000":
        files = qcd1000Files
        loc = qcd1000Loc
    elif opts.sample == "qcdFlatA":
        files = qcdFlatAFiles
        loc = qcdFlatALoc
    elif opts.sample == "qcdFlatB":
        files = qcdFlatBFiles
        loc = qcdFlatBLoc
    elif opts.sample == "qcdFlatEMA":
        files = qcdFlatEMAFiles
        loc = qcdFlatEMALoc
    elif opts.sample == "qcdFlatEMB":
        files = qcdFlatEMBFiles
        loc = qcdFlatEMBLoc
    elif opts.sample == "qcdFlatBBEM":
        files = qcdFlatBBEMFiles
        loc = qcdFlatBBEMLoc

    elif opts.sample == 'singleTop':
        loc = singleTopLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'singleAntiTop':
        loc = singleAntiTopLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    elif opts.sample == 'tZqll':
        loc = tZqllLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    elif opts.sample == 'WW':
        loc = WWLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'ZZTo2L2Q':
        loc = ZZTo2L2QLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
        
    elif opts.sample == "wJets":
        files = wJetsFiles
        loc = wJetsLoc

    elif opts.sample == "VVTo2L2Nu":
        files = VVTo2L2NuFiles
        loc = VVTo2L2NuLoc

    elif opts.sample == "ttHTobb":
        files = ttHTobbFiles
        loc = ttHTobbLoc
    elif opts.sample == "ttHToNonbb":
        files = ttHToNonbbFiles
        loc = ttHToNonbbLoc

    elif opts.sample == 'ZH_HToBB_ZToQQ_M125':
        loc = ZH_HToBB_ZToQQ_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'WminusH_HToBB_WToQQ_M125':
        loc = WminusH_HToBB_WToQQ_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'WplusH_HToBB_WToQQ_M125':
        loc = WplusH_HToBB_WToQQ_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'VHToNonbb':
        loc = VHToNonbbLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'WplusH_HToBB_WToLNu_M125':
        loc = WplusH_HToBB_WToLNu_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'WminusH_HToBB_WToLNu_M125':
        loc = WminusH_HToBB_WToLNu_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
    elif opts.sample == 'ZH_HToBB_ZToLL_M125':
        loc = ZH_HToBB_ZToLL_M125Loc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    elif opts.sample == 'signal_13TeV':
        loc = signal_13TeVLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    print(len(files), " files found, begining job submission from file ", int(opts.first))

    for i, f in enumerate(files):
        if i >= int(opts.first):
            makeJOFile(loc+f, i, opts)
    make_sub_file(opts)
