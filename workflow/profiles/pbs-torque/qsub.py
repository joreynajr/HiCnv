#!/mnt/BioHome/jreyna/software/anaconda3/envs/hic_tls/bin/python
"""
A custom qsub submission script for snakemake. To use use a snakemake
command like "qsub --cluster <path>/qsub.py.

Do not remove the python shebang, snakemake uses this internally (if I
remember this correctly).
"""

import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)


for k1 in job_properties:
    print(k1, job_properties[k1])
exit


def ifelse_default(dt, key, default):
    """
        If key in dictionary they return that value, else
        return a default value.
    """
    if key in dt:
        return(dt[key])
    else:
        return(default)


# access the resource properties defined for each rule
if 'resources' not in job_properties:
    mem = 4000
    nodes = 1
    ppn = 1
else:
    mem = ifelse_default(job_properties['resources'], 'mem_mb', 4000)
    nodes = ifelse_default(job_properties['resources'], 'nodes', 1)
    ppn = ifelse_default(job_properties['resources'], 'ppn', 1)

# checking for the presence of log files
try:
    out = job_properties['log'][0]

    # also make the log directory tree + log file
    bn = os.path.dirname(out)
    os.makedirs(bn, exist_ok=True)
    #open(out, 'a').close()

except KeyError:
    raise("Forgot to set out log file for current rule.")

# creating and submitting the final qsub command
cmd = 'qsub -l walltime=200:00:00,mem={}mb,nodes={}:ppn={} -o {} -e {} {}'
cmd = cmd.format(mem, nodes, ppn, out, out, jobscript)
os.system(cmd)
