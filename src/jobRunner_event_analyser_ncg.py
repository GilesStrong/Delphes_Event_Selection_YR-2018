'''Giles Strong (giles.strong@outlook.com)'''

from __future__ import division
import os
import os.path
import optparse
from data import *

softDir = "/lstore/cms/giles/Simple_Delphes_Event_Selection/src/"
outDir = '/lstore/cms/giles/HLStudies/'

def makeJOFile(inputFile, uid, opts):
    outputFile = outDir + opts.sample + "_" + str(uid)
    cmd = "./delphes_event_selection "
    cmd += "-i " + inputFile
    cmd += " -o " + outputFile
    cmd += " -d " + str(opts.debug[-1])
    joName = "analysis_" + str(uid) + ".job"
    joFile = open(joName, "w")
    joFile.write("echo Beginning\ job\n")
    joFile.write("module load gcc-4.8.3\n")
    joFile.write("module load python-2.7.11\n")
    joFile.write("export PATH=/lstore/cms/giles/programs/bin:$PATH\n")
    joFile.write("export LD_LIBRARY_PATH=/lstore/cms/giles/programs/lib64/:/lstore/cms/giles/programs/lib/:/lstore/cms/giles/programs/lib/root/:/lstore/cms/giles/programs/delphes/:$LD_LIBRARY_PATH\n")
    joFile.write("source /lstore/cms/giles/programs/bin/thisroot.sh\n")
    joFile.write("export X509_USER_PROXY=/lstore/cms/giles/x509up_u5020023\n")
    joFile.write("cd " + softDir + "\n")
    joFile.write("echo Paths\ set\n")
    joFile.write(cmd + "\n")
    joFile.close()
    sub = "qsub " + joName
    print "Submitting: " + sub
    os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-s", "--sample", dest = "sample", action = "store", help = "Sample to analyse")
    parser.add_option("-d", "--debug", dest = "debug", action = "store", default = [0], help = "Run in debug mode. {0,1}, default: 0")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = ["normal"], help = "Queue to run jobs. Default: normal")
    opts, args = parser.parse_args()

    os.system("voms-proxy-init -voms cms -valid 72:00")
    os.system("cp /tmp/x509up_u5020023 /afs/cern.ch/user/g/gstrong/x509up_u5020023")
    
    if opts.sample == "signal":
        files = signalFiles
        loc = signalLoc
    elif opts.sample == "ttbar":
        files = ttbarFiles
        loc = ttbarLoc

    elif opts.sample == "ttbar_13TeV":
        loc = ttbar_13TeVLoc
        files = [x[x.rfind("/")+1:] for x in glob.glob(loc + "/*.root")]

    for i, f in enumerate(files):
        makeJOFile(loc+f, i, opts)
        #break
