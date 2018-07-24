'''Giles Strong (giles.strong@outlook.com)'''

from __future__ import division
import os
import os.path
import optparse
from data import *

softDir = "/lstore/cms/giles/Simple_Delphes_Event_Selection/src/"
userStorage = '/lstore/cms/giles/HLStudies/'

def makeJOFile(inputFile, uid, opts):
    outputFile = opts.sample + "_" + str(uid)
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
    joFile.write("source lstore/cms/giles/programs/bin/thisroot.sh\n")
    joFile.write("export X509_USER_PROXY=/lstore/cms/giles/x509up_uXXXX\n")
    joFile.write("cd " + softDir + "\n")
    joFile.write("echo Paths\ set\n")
    joFile.write(cmd + "\n")
    joFile.close()
    sub = "qsub " + joName
    print "Submitting: " + sub
    #os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-s", "--sample", dest = "sample", action = "store", help = "Sample to analyse")
    parser.add_option("-d", "--debug", dest = "debug", action = "store", default = [0], help = "Run in debug mode. {0,1}, default: 0")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = ["normal"], help = "Queue to run jobs. Default: normal")
    opts, args = parser.parse_args()
    
    if opts.sample == "signal":
        files = signalFiles
        loc = signalLoc
    elif opts.sample == "ttbar":
        files = ttbarFiles
        loc = ttbarLoc

    for i, f in enumerate(files):
        makeJOFile(loc+f, i, opts)
        break
