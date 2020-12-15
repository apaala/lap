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
parser.add_option("-s", "--sample", dest="samples",help="sample target list", metavar="FILE")
parser.add_option("-v", "--rsample", dest="rsamples",help="sample reference list", metavar="FILE")
parser.add_option("-p", "--refmap", dest="refmap",help="Reference panel sample population classification map", metavar="FILE")
parser.add_option("-m", "--map", dest="mapf",help="mapfile", metavar="FILE") 
(options, args) = parser.parse_args() 

client = docker.from_env()
vol_dir="/local/devel/local_ancestry_pipeline/volume_lap/"
name=os.path.basename(options.name)


def prep_mod1(target,ref, fname,chrom,sample,rsample,refmap,mapf):
    bgcmd="bgzip -c tmp/"+ target+" > /tmp/"+fname+".vcf.gz"
    print("1: starting bgzip")
    bgzipf="/tmp/"+fname+".vcf.gz"
    print(bgzipf)
    tbxf=bgzipf+".tbi"
    #print(tbxf)
    client.containers.run("dockerbiotools/bcftools:latest", command=[bgcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    #tbx_cmd="tabix -p vcf "+bgzipf
    bgcmd="bgzip -d "+bgzipf
    client.containers.run("dockerbiotools/bcftools:latest", command=[bgcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    bgzipf="/tmp/"+fname+".vcf"
    tbx_cmd="tabix -p vcf "+bgzipf
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbx_cmd], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})

    print("2: running tabix")
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o /tmp/"+fname+"filt_out.vcf.gz "+ bgzipf

    print("3: filter with bcftools")
    client.containers.run("dockerbiotools/bcftools:latest", command=[bcfcmd], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    bcff="/tmp/"+fname+"filt_out.vcf.gz"
    ###Edit for chrom=1

    print("4: starting conform")
    conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref=/tmp/"+ref+" gt="+bcff+" chrom="+str(chrom)+" match=POS out=/tmp/conform"+fname
    #conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref="+ref+" gt="+bgzipf+" chrom="+str(chrom)+" out=/tmp/conform"+fname
    print(conformcmd) 
    client.containers.run("apaala/beagle:0.1", command=conformcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    conf="/tmp/conform"+fname+".vcf.gz"
    print("5: Tabix on conformed file")

    tbxcon="tabix -p vcf "+conf
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    
    print("6: merging with bcftools")
    mergec="bcftools merge -m none -o /tmp/"+fname+"_"+str(chrom)+"_merged.vcf.gz -Oz "+conf+" /tmp/"+ref
    client.containers.run("dockerbiotools/bcftools:latest", command=[mergec], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    mergef="/tmp/"+fname+"_"+str(chrom)+"_merged.vcf.gz"
    
    print("7: running tabix on merged file")
    tbxcon="tabix -p vcf "+mergef
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    
    print("8: conform on merged")
    conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref=/tmp/"+ref+" gt="+mergef+" chrom="+str(chrom)+" match=POS out=/tmp/merge_conform"+fname
    client.containers.run("apaala/beagle:0.1", command=conformcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    mcf="/tmp/merge_conform"+fname+".vcf.gz"
    print(mcf)
    tbxcon="tabix -p vcf "+mcf
    
    print("9: Tabix on merged conformed file")
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    ##Change cmd to reflect beagle map file
    #cmap="/home/map38/plink."+str(chrom)+".GRCh38.map"
    cmap="plink.chr22.GRCh38.map"
    
    print("10: running beagle")
    beagcmd="java -jar /home/beagle.24Mar20.5f5.jar gt="+mcf+" ref=/tmp/"+ref+" map=/tmp/"+cmap+" out=/tmp/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    #beagcmd="beagle gt="+conf+" ref="+ref+" map="+cmap+" out=/tmp/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    client.containers.run("apaala/beagle:0.1",command=beagcmd,volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    
    print("11: Tabix on phased")
    phased="/tmp/beagle_"+fname+".phased.vcf.gz"
    tbxcon="tabix -p vcf "+phased
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    ######make vcf files instead so we dont need bcftools in rfmix docker

    print("12:generating vcf files")
    prepcmd="bcftools view -S /tmp/"+sample+" -Ov -o /tmp/"+fname+"target_prepped.vcf "+phased
    client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    target_bcf="/tmp/"+fname+"target_prepped.vcf"
    prepcmd="bcftools view -S /tmp/"+rsample+" -Ov -o /tmp/"+fname+"_ref_prepped.vcf /tmp/"+ref
    client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    ref_bcf="/tmp/"+fname+"_ref_prepped.vcf"
    #prepcmd="bgzip -c "+target_bcf+" > "+target_bcf+".gz"
    #client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    #prepcmd="bcftools index "+target_bcf+".gz"
    #client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    #prepcmd="bgzip -c "+ref_bcf+" > "+ref_bcf+".gz"
    #client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    #prepcmd="bcftools index "+ref_bcf+".gz"
    #client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    print("13:Running RFMIX2")
    chromosome=22 ###make it into a loop!
    #/usr/local/packages/rfmix-2.0.0/bin/rfmix -f B1_prep.bcf -r B1_prep_ref.bcf -m rfmix_sample.txt -g plink.chr22.GRCh38.map -o oct_rfmix_B1.phased.thread.vcf.gz --chromosome=22 -e 1 --n-threads=8 --crf-weight=3
    prepcmd="rfmix -f "+target_bcf+" -r "+ref_bcf+" -g /tmp/"+refmap+" -m /tmp/"+mapf+" -o /tmp/"+fname+"_rfmix"+" --chromosome="+str(chrom)+" -e 1 --n-threads=8 --crf-weight=3"
    ######Do beagle through rfmix in loop with chromosomes.
#    client.containers.run("biocontainers/rfmix", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    client.containers.run("quay.io/biocontainers/rfmix:2.03.r0.9505bfa--he1b5a44_1", command=prepcmd , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})

tmp_reference=os.path.basename(options.reference)
tmp_target=os.path.basename(options.target)
prep_mod1(tmp_target, tmp_reference, options.name, options.chrom,options.samples,options.rsamples,options.refmap,options.mapf)




def prep_mod_local(target,ref, fname, outdir):
    bgcmd="bgzip -c "+target+"> "+outdir+"/"+fname+"_"+os.path.basename(target)+".gz" 
    bgzipf=outdir+"/"+fname+"_"+os.path.basename(target)+".gz"
    print(bgcmd)
    tabxcmd="tabix -p vcf "+bgzipf
    print(tabxcmd)
    tbxf=bgzipf+".tbi"
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o "+outdir+"/"+fname+ "_out.vcf.gz "+ bgzipf
    print(bcfcmd)
    bcffile=outdir+"/"+fname+ "_out.vcf.gz"
    conformcmd="java -jar /usr/local/packages/conform-gt/conform-gt.jar ref="+ref+" gt="+bcffile+" chrom=1"+" out=conform"+fname
    print(conformcmd)
    conformf=outdir+ "/conform" +fname+".vcf.gz" 
    conformtbx="tabix "+conformf
    print(conformtbx)
    print(conformf)
    bcfm="bcftools merge -m none -o mergeout_"+ fname+".vcf.gz -Oz --threads 2 "+ conformf+ " "+ref
    print(bcfm)
    mfile="mergeout_"+ fname+".vcf.gz"
    bcfc="bcftools view -h "+mfile+" | tail -1 | tr \'\\t\' \'\\n\' | grep 'NWD' > temp_"+mfile
    print(bcfc)
    ccmd="cat temp_"+mfile+ " "+samples+" > samples_chr"+fname+".txt"
    print(ccmd)
