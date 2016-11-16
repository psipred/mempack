import sys
import re
from subprocess import call

topology_string = ''
with open(sys.argv[1], "r") as f:
    for line in f:
        if line.startswith("Topology:"):
            topology_string = line
f.close()
topology_string = re.sub("Topology:\s+", "", topology_string)
topology_string = topology_string.rstrip()
topology_string = re.sub("-", ",", topology_string)
exe = sys.argv[2]
mtx = sys.argv[3]
mdir = sys.argv[4]

args = [exe, "-mtx", "1", "-t", topology_string, '-w', mdir, mtx]
print(" ".join(args))
# call(args)
