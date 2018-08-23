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

    joName = "analysis_" + str(uid) + ".job"
    joFile = open(joName, "w")
    joFile.write("echo Beginning\ job\n")
    joFile.write("export HOME=/afs/cern.ch/user/g/gstrong/\n")
    joFile.write("source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh\n")
    joFile.write("source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
    joFile.write("export LD_LIBRARY_PATH=$HOME/programs/delphes/:$LD_LIBRARY_PATH\n")
    joFile.write("export X509_USER_PROXY=$HOME/x509up_u5020023\n")
    joFile.write("cd " + softDir + "\n")
    joFile.write("echo Paths\ set\n")
    joFile.write(cmd + "\n")
    joFile.close()

    sub = "bsub -q " + opts.queue + " " + joName
    print "Submitting: " + sub
    os.system("chmod 744 " + joName)
    os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-s", "--sample", dest = "sample", action = "store", help = "Sample to analyse")
    parser.add_option("-d", "--debug", dest = "debug", action = "store", default = 0, help = "Run in debug mode. {0,1}, default: 0")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = "8nh", help = "Queue to run jobs. Default: normal")
    parser.add_option("-f", "--first", dest = "first", action = "store", default = 0, help = "First job to run. Default: 0")
    opts, args = parser.parse_args()

    os.system("voms-proxy-init -voms cms -valid 72:00")
    os.system("cp /tmp/x509up_u61049 /afs/cern.ch/user/g/gstrong/x509up_u5020023")
    
    if opts.sample == "signal":
        files = signalFiles
        loc = signalLoc

    elif opts.sample == "ttbar":
        loc = ttbarLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]
        
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

    print(len(files), " files found, begining job submission from file ", int(opts.first))

    for i, f in enumerate(files):
        if i >= int(opts.first):
            makeJOFile(loc+f, i, opts)
        #break
