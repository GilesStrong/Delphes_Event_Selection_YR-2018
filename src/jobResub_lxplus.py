import os, glob
failures = [x[:x.rfind("/")] for x in glob.glob("*/STDOUT") if "Data saved" not in open(x).read()]
print len(failures), failures