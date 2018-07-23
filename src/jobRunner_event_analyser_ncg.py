'''Giles Strong (giles.strong@outlook.com)'''

from __future__ import division
import os
import os.path
import optparse
from data import *

pwd = os.getcwd()
userDir = "/home/t3cms/giles/"
userStorage = '/lstore/cms/giles/HLStudies/'

def makeJOFile(inputFile, uid, opts):
    outputFile = userStorage + opts.sample + "/" + opts.sample + "_" + str(uid)
    cmd = "./delphes_event_selection "
    cmd += "-i " + inputFile
    cmd += " -o " + outputFile
    cmd += " -d " + str(opts.debug[-1])
    joName = name + "analysis_" + str(uid) + ".job"
    joFile = open(joName, "w")
    joFile.write("echo Beginning\ job\n")
    joFile.write("source " + userDir + ".bashrc\n")
    joFile.write("export X509_USER_PROXY=" + userDir + "x509up_uXXXX")
    joFile.write("echo Paths\ set\n")
    joFile.write("cd " + pwd + "\n")
    joFile.write(cmd + "\n")
    joFile.close()
    sub = "qsub " + joName + " -q " + opts.queue[-1]
    print "Submitting: " + sub
    os.system(sub)

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
