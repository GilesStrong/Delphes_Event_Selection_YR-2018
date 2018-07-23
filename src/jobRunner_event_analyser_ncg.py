'''Giles Strong (giles.strong@outlook.com)'''

from __future__ import division
import os
import os.path
import optparse

pwd = os.getcwd()
userDir = "/home/t3cms/giles/"

def makeJOFile(start, opts):
    stop = start+int(opts.population[-1])
    numbers = str(start) + "_to_" + str(stop)
    inputFile = opts.input[-1] + numbers + ".root"
    name = opts.input[-1].split("/")[-1]
    if not "default" in opts.output[-1]:
        outputFile = opts.output[-1]
    else:
        outputFile = name
    outputFile += numbers
    cmd = "./delphes_event_selection "
    cmd += "-i " + inputFile
    cmd += " -o " + outputFile
    cmd += " -t " + str(opts.truth[-1])
    cmd += " -r " + str(opts.response[-1])
    cmd += " -s " + str(opts.selection[-1])
    cmd += " -d " + str(opts.debug[-1])
    cmd += " -m " + str(opts.mva[-1])
    joName = name + "analysis_" + numbers + ".job"
    joFile = open(joName, "w")
    joFile.write("echo Beginning\ job\n")
    joFile.write("source " + userDir + ".bashrc\n")
    joFile.write("echo Paths\ set\n")
    joFile.write("cd " + pwd + "\n")
    joFile.write(cmd + "\n")
    joFile.close()
    sub = "qsub " + joName + " -q " + opts.queue[-1]
    print "Submitting: " + sub
    os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-i", "--input", dest = "input", action = "append", help = "Basic input name.")
    parser.add_option("-o", "--output", dest = "output", action = "append", default = ["default"], help = "Output name. Default: base on input name")
    parser.add_option("-n", "--nFiles", dest = "nFiles", action = "append", default = [100], help = "Number of files to process. [0, inf], default: 100")
    parser.add_option("-p", "--population", dest = "population", action = "append", default = [100000], help = "Number of events per file. [0, inf], default: 100000")
    parser.add_option("-a", "--first", dest = "first", action = "append", default = [0], help = "First event event to process. [0, inf], default: 0")
    parser.add_option("-t", "--truth", dest = "truth", action = "append", default = [0], help = "Use MC truth cut. {0,1}, default: 0")
    parser.add_option("-s", "--selection", dest = "selection", action = "append", default = [1], help = "Run event selection. {0,1}, default: 1")
    parser.add_option("-d", "--debug", dest = "debug", action = "append", default = [0], help = "Run in debug mode. {0,1}, default: 0")
    parser.add_option("-m", "--mva", dest = "mva", action = "append", default = [0], help = "Output information for MVA selection. {0,1}, default: 0")
    parser.add_option("-q", "--queue", dest = "queue", action = "append", default = ["normal"], help = "Queue to run jobs. Default: normal")
    opts, args = parser.parse_args()
    nEvents = int(opts.nFiles[-1])*opts.population[-1]
    runNumbers = range(int(opts.first[-1]), int(opts.first[-1])+nEvents, int(opts.population[-1]))
    print "Setting " + str(len(runNumbers)) + " jobs to run on queue " + opts.queue[-1] + " over " + str(nEvents/len(runNumbers)) + " events each,"
    print "between event numbers " + str(opts.first[-1]) + " and " + str(int(opts.first[-1])+nEvents) + " of file " + opts.input[-1]
    if raw_input("Continue? [y/n]: ").lower().strip() == "y":
        for i in runNumbers:
            makeJOFile(i, opts)
