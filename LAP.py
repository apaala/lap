import docker
import getopt, sys, os 
from optparse import OptionParser 

parser = OptionParser() 
parser.add_option("-t", "--target", dest="target",help="Path to target", metavar="FILE") 
parser.add_option("-r", "--reference", dest="reference",help="Path to reference", metavar="FILE") 
parser.add_option("-n", "--name", dest="name",help="Name for output file", metavar="NAME")
(options, args) = parser.parse_args() 

client = docker.from_env()
vol_dir="/Users/apaala/Docker/LAP_vol/"
def prep_mod1(target, fname):
    bgcmd="bgzip -c "+ target +">/tmp/py_test/"+fname+"_out.vcf.gz"
    print(bgcmd)
    bgzipf=target+".gz"
    tbxf=bgzipf+".tbi"
    client.containers.run("dockerbiotools/bcftools:latest", bgcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    tbx_cmd="tabix -p vcf "+bgzipf
    client.containers.run("dockerbiotools/bcftools:latest", tbx_cmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o /tmp/"+fname+"_out.vcf.gz "+ target
    client.containers.run("dockerbiotools/bcftools:latest", bcfcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})

prep_mod1(options.target, options.name)

#client.containers.run("apaala/beagle:0.1", "sh /home/beagle.example")
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example_ed', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})


