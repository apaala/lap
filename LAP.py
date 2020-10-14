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
vol_dir="/local/devel/local_ancestry_pipeline/volume_lap/"
name=os.path.basename(options.name)


def prep_mod1(target,ref, fname,chrom):
    bgcmd="bgzip -c "+ target+" > /tmp/"+fname+".vcf.gz"
    print("1: starting bgzip")
    bgzipf="/tmp/"+fname+".vcf.gz"
    print(bgzipf)
    tbxf=bgzipf+".tbi"
    #print(tbxf)
    client.containers.run("dockerbiotools/bcftools:latest", command=[bgcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    tbx_cmd="tabix -p vcf "+bgzipf
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbx_cmd], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})

    print("2: running tabix")
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o /tmp/"+fname+"filt_out.vcf.gz "+ bgzipf

    print("3: filter with bcftools")
    client.containers.run("dockerbiotools/bcftools:latest", command=[bcfcmd], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    bcff="/tmp/"+fname+"filt_out.vcf.gz"
    ###Edit for chrom=1

    print("4: starting conform")
    ###FOR DEMO NOT USING FILTERED FILE
    conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref="+ref+" gt="+bcff+" chrom="+str(chrom)+" match=POS out=/tmp/conform"+fname
    #conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref="+ref+" gt="+bgzipf+" chrom="+str(chrom)+" out=/tmp/conform"+fname
    print(conformcmd) 
    client.containers.run("apaala/beagle:0.1", command=conformcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    conf="/tmp/conform"+fname+".vcf.gz"
    print("5: Tabix on conformed file")

    #print("5: merging with bcftools")
    ###ADD bcftools merge cmd 
    #    $bcftools merge -m none -o $temp -Oz --threads 2 $conformed_file $peru_1KG  
    tbxcon="tabix -p vcf "+conf
    #print("Tabix on conformed file")
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    ####Added --force-samples check with team to see if thats ok
    
    print("6: merging with bcftools")
    mergec="bcftools merge -m none -o /tmp/"+fname+"_"+str(chrom)+"_merged.vcf.gz -Oz "+conf+" "+ref
    client.containers.run("dockerbiotools/bcftools:latest", command=[mergec], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    mergef="/tmp/"+fname+"_"+str(chrom)+"_merged.vcf.gz"
    
    print("7: running tabix on merged file")
    tbxcon="tabix -p vcf "+mergef
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    
    print("8: conform on merged")
    conformcmd="java -jar /home/conform-gt.24May16.cee.jar ref="+ref+" gt="+mergef+" chrom="+str(chrom)+" match=POS out=/tmp/merge_conform"+fname
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
    ####EDITED FOR DEMO ONLY, NOT MERGED FILE
    beagcmd="java -jar /home/beagle.24Mar20.5f5.jar gt="+mcf+" ref="+ref+" map=/tmp/"+cmap+" out=/tmp/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    #beagcmd="beagle gt="+conf+" ref="+ref+" map="+cmap+" out=/tmp/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    client.containers.run("apaala/beagle:0.1",command=beagcmd,volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    
    print("11: Tabix on phased")
    phased="/tmp/beagle_"+fname+".phased.vcf.gz"
    tbxcon="tabix -p vcf "+phased
    client.containers.run("dockerbiotools/bcftools:latest", command=[tbxcon], volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})
    #$bcftools view -S ^$ref -Ob -o $prefix.target.chr$chr.bcf $vcf
    #$bcftools view -S $ref -Ob -o ref.chr$chr.bcf $vcf
    #$bcftools index ref.chr$chr.bcf    
    #$bcftools index $prefix.target.chr$chr.bcf
    #bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o /tmp/"+fname+"filt_out.vcf.gz "+ bgzipf
    print("12:generating bcf files")
    prepcmd="bcftools view -S "+options.samples+" -Ob -o /tmp/"+fname+"target_prepped.bcf "+phased
    client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    target_bcf="/tmp/"+fname+"target_prepped.bcf"
    prepcmd="bcftools view -S "+options.samples+" -Ob -o /tmp/"+fname+"_ref_prepped.bcf "+ref
    client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    ref_bcf="/tmp/"+fname+"_ref_prepped.bcf"
    prepcmd="bcftools index "+target_bcf
    client.containers.run("dockerbiotools/bcftools:latest", command=[prepcmd] , volumes={vol_dir: {'mode': 'rw', 'bind': '/tmp'}})
    prepcmd="bcftools index "+ref_bcf



tmp_reference="/tmp/"+os.path.basename(options.reference)
tmp_target="/tmp/"+os.path.basename(options.target)
prep_mod1(tmp_target, tmp_reference, options.name, options.chrom)
#prep_mod_local(options.target, options.name)
#def determine_map(chrom):
#    hold=client.containers.run("apaala/beagle:0.1", command='ls *.

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
    #conformcmd="java -jar /usr/local/packages/conform-gt/conform-gt.jar ref="+ref+" gt="+bcffile+" chrom=1"+" out=conform"+fname
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
#   bccmd="
    #ref=/local/chib/toconnor_grp/LARGE-PD/1KG/HG38/chr$chr.1kg.HG38.vcf.gz 
    #java -jar $conform ref=$ref gt=$merged chrom=chr$chr out=$prefix.merged.chr$chr.CONFORM match=POS 
    #tabix $prefix.merged.chr$chr.CONFORM.vcf.gz 
    #print(bgcmd)
    #print(bgzipf)
    #print(tabxcmd)
    #print(bcfcmd)


#def beagle_rfmix(target, ref, fname, chrn):
    #rfmix=/home/dloesch/WORKSPACE/rfmix.1KG_samples.txt 
    #$bcftools view -S^$rfmix  -e 'GT[*] = "mis"' -o $new_ref -Oz $ref 
    #tabix $new_ref #index if needed 
 
    #run beagle 
#   beagle=/usr/local/packages/beagle-5.0/beagle.jar 
 
    #out=$prefix.merged.chr$chr.PHASED 
    #map=/local/chib/toconnor_grp/LARGE-PD/genetic_maps/plink.chr$chr.v2.GRCh38.map 
 
    #java -Xmx50g -jar $beagle gt=$gt ref=$new_ref out=$out map=$map chrom=chr$chr impute=false 

#prep_mod_local(options.target,options.reference, options.name, options.outdir) 


#client.containers.run("apaala/beagle:0.1", "sh /home/beagle.example")
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})
#client.containers.run('apaala/beagle:0.1','sh /tmp/beagle.example_ed', volumes={'/Users/apaala/Docker/LAP_vol/':{'bind':'/tmp', 'mode':'rw'}})


