import sys
import re
from subprocess import call

# Usage:
# python web_wrapper.py /tmp/21ee4184-ac09-11e6-affd-989096c13ee6/21ee4184-ac09-11e6-affd-989096c13ee6.memsat_svm /opt/Code/mempack/run_mempack.pl /tmp/21ee4184-ac09-11e6-affd-989096c13ee6/21ee4184-ac09-11e6-affd-989096c13ee6.mtx /tmp/21ee4184-ac09-11e6-affd-989096c13ee6/ /opt/Code/mempack/

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
out = sys.argv[4]
mdir = sys.argv[5]

args = [exe, "-mtx", "1", "-t", topology_string, '-w', mdir, "-j", out, mtx]
print(" ".join(args))
call(args)
