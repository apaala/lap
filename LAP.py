import docker
import getopt, sys, os 
from optparse import OptionParser 
import tabix

parser = OptionParser() 
parser.add_option("-t", "--target", dest="target",help="Path to target", metavar="FILE") 
parser.add_option("-r", "--reference", dest="reference",help="Path to reference", metavar="FILE") 
parser.add_option("-n", "--name", dest="name",help="Name for output file", metavar="NAME")
parser.add_option("-o", "--out", dest="outdir",help="output dir", metavar="Out")
parser.add_option("-c", "--chr", dest="chrom",help="chromosome", metavar="NUM") 
parser.add_option("-s", "--sample", dest="samples",help="samplelist", metavar="FILE")
parser.add_option("-m", "--map", dest="map",help="mapfile", metavar="FILE") 
(options, args) = parser.parse_args() 

client = docker.from_env()
vol_dir="/Users/apaala/Docker/LAP_vol/"
name=os.path.basename(options.name)


def prep_mod1(target, fname):
    bgcmd="bgzip -c "+ target+" > /tmp/py_test/"+fname+".vcf.gz"
    print(bgcmd)
    bgzipf=target+".gz"
    print(bgzipf)
    tbxf=bgzipf+".tbi"
    print(tbxf)
    ###Something wrong in env here
    #client.containers.run("dockerbiotools/bcftools:latest", "ls /tmp" , volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})
    client.containers.run("dockerbiotools/bcftools:latest", bgcmd, volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})
    tbx_cmd="tabix -p vcf "+bgzipf
    client.containers.run("dockerbiotools/bcftools:latest", tbx_cmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o /tmp/"+fname+"_out.vcf.gz "+ target
    client.containers.run("dockerbiotools/bcftools:latest", bcfcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})

#prep_mod1(options.target, options.name)
#prep_mod_local(options.target, options.name)

def prep_mod_local(target,ref, fname, outdir):
    bgcmd="bgzip -c "+target+"> "+outdir+"/"+fname+"_"+os.path.basename(target)+".gz" 
    bgzipf=outdir+"/"+fname+"_"+os.path.basename(target)+".gz"
    print(bgcmd)
#    os.system(cmd)

    tabxcmd="tabix -p vcf "+bgzipf
    print(tabxcmd)
#    os.system(tabxcmd)
    tbxf=bgzipf+".tbi"
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o "+outdir+"/"+fname+ "_out.vcf.gz "+ bgzipf
    print(bcfcmd)
#    os.system(bcfcmd)
    bcffile=outdir+"/"+fname+ "_out.vcf.gz"
    conformcmd="java -jar /usr/local/packages/conform-gt/conform-gt.jar ref="+ref+" gt="+bcffile+" chrom=1"+" out=conform"+fname
    print(conformcmd)
    conformf=outdir+ "/conform" +fname+".vcf.gz" 
    conformtbx="tabix "+conformf
    print(conformtbx)
    print(conformf)
#    $bcftools merge -m none -o $temp -Oz --threads 2 $conformed_file $peru_1KG
    bcfm="bcftools merge -m none -o mergeout_"+ fname+".vcf.gz -Oz --threads 2 "+ conformf+ " "+ref
    #$bcftools view -h $vcf | tail -1 | tr '\t' '\n' | grep "NWD" > temp.chr$chr.txt
    print(bcfm)
    mfile="mergeout_"+ fname+".vcf.gz"
    bcfc="bcftools view -h "+mfile+" | tail -1 | tr \'\\t\' \'\\n\' | grep 'NWD' > temp_"+mfile
    print(bcfc)
    #cat temp.chr$chr.txt /home/dloesch/WORKSPACE/rfmix.ref_samples.txt > samples.chr$chr.txt 
    ccmd="cat temp_"+mfile+ " "+samples+" > samples_chr"+fname+".txt"
    print(ccmd)
    #$bcftools view -S samples.chr$chr.txt -e 'GT[*] = "mis"' -o $merged -Oz $temp 
#    bccmd="
    #ref=/local/chib/toconnor_grp/LARGE-PD/1KG/HG38/chr$chr.1kg.HG38.vcf.gz 
    #java -jar $conform ref=$ref gt=$merged chrom=chr$chr out=$prefix.merged.chr$chr.CONFORM match=POS 
    #tabix $prefix.merged.chr$chr.CONFORM.vcf.gz 
    #print(bgcmd)
    #print(bgzipf)
    #print(tabxcmd)
    #print(bcfcmd)


def beagle_rfmix(target, ref, fname, chrn):
    #rfmix=/home/dloesch/WORKSPACE/rfmix.1KG_samples.txt 
    #$bcftools view -S^$rfmix  -e 'GT[*] = "mis"' -o $new_ref -Oz $ref 
    #tabix $new_ref #index if needed 
 
    #run beagle 
#   beagle=/usr/local/packages/beagle-5.0/beagle.jar 
 
    #out=$prefix.merged.chr$chr.PHASED 
    #map=/local/chib/toconnor_grp/LARGE-PD/genetic_maps/plink.chr$chr.v2.GRCh38.map 
 
    #java -Xmx50g -jar $beagle gt=$gt ref=$new_ref out=$out map=$map chrom=chr$chr impute=false 

prep_mod_local(options.target,options.reference, options.name, options.outdir) 


#client.containers.run("apaala/beagle:0.1", "sh /home/beagle.example")
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example_ed', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})


