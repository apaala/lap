#import docker
import getopt, sys, os
from optparse import OptionParser
#import tabix
import os  
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
 
name=os.path.basename(options.name)

def prep_mod1(target,ref,fname,outd,chrom,sample,rsample,refmap,mapf):
    bgcmd="bgzip -c "+target+" > "+outd+"/"+fname+".vcf.gz"
    print("1: starting bgzip")
    print(bgcmd)
    os.system(bgcmd)
    bgzipf=outd+"/"+fname+".vcf.gz"
    print(bgzipf)
    tbxf=bgzipf+".tbi"
    bgcmd="bgzip -d "+bgzipf
    os.system(bgcmd)
    bgzipf=outd+"/"+fname+".vcf"
    tbx_cmd="tabix -p vcf "+bgzipf
    os.system(tbx_cmd)
    print("2: running tabix")
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o "+outd+fname+"filt_out.vcf.gz "+ bgzipf
    os.system(bcfcmd)
    print("3: filter with bcftools")
    bcff=outd+"/"+fname+"filt_out.vcf.gz"
    print(bcfcmd)
    print("4: starting conform")
    conformcmd="java -jar conform-gt.24May16.cee.jar ref="+ref+" gt="+bcff+" chrom="+str(chrom)+" match=POS out="+outd+"conform"+fname   
    print(conformcmd)
    os.system(conformcmd)                                                                                                                            #client.containers.run("apaala/beagle:0.1", command=conformcmd, volumes={vol_dir:{'bind':'/tmp', 'mode':'rw'}})                               
    conf=outd+"/conform"+fname+".vcf.gz"                                                                                                        
    print("5: Tabix on conformed file")                                                                                                          
    tbxcon="tabix -p vcf "+conf                                                                                                                  
    os.system(tbxcon)
    print("6: merging with bcftools")                                                                                                            
    mergec="bcftools merge -m none -o "+outd+"/"+fname+"_"+str(chrom)+"_merged.vcf.gz -Oz "+conf+" "+ref                                        
    mergef=outd+"/"+fname+"_"+str(chrom)+"_merged.vcf.gz"                                                                                         
    os.system(mergec)                                                                                                                              
    print("7: running tabix on merged file")                                                                                                     
    tbxcon="tabix -p vcf "+mergef                                                                                                                
    os.system(tbxcon)
    print("8: conform on merged")                                                                                                                
    conformcmd="java -jar conform-gt.24May16.cee.jar ref="+ref+" gt="+mergef+" chrom="+str(chrom)+" match=POS out="+outd+"/merge_conform"+fname         
    os.system(conformcmd)
    mcf=outd+"/merge_conform"+fname+".vcf.gz"                                                                                                     
    print(mcf)                                                                                                                                   
    tbxcon="tabix -p vcf "+mcf                                                                                                                   
    os.system(tbxcon)                                                                                                                                  
    print("9: Tabix on merged conformed file")                                                                                                   
    cmap=mapf                                                                                                                
    print("10: running beagle")                                                                                                                  
    beagcmd="beagle gt="+mcf+" ref="+ref+" map="+cmap+" out="+outd+"/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    os.system(beagcmd)
    print("11: Tabix on phased")                                                                                                                 
    phased=outd+"/beagle_"+fname+".phased.vcf.gz"                                                                                                 
    tbxcon="tabix -p vcf "+phased                                                                                                                
    os.system(tbxcon)
    print("12:generating vcf files")                                                                                                             
    prepcmd="bcftools view -S "+sample+" -Ov -o "+outd+"/"+fname+"target_prepped.vcf "+phased                                                   
    os.system(prepcmd)
    target_bcf=outd+"/"+fname+"target_prepped.vcf"                                                                                                
    prepcmd="bcftools view -S "+rsample+" -Ov -o "+outd+"/"+fname+"_ref_prepped.vcf "+ref
    ref_bcf=outd+"/"+fname+"_ref_prepped.vcf"
    os.system(prepcmd)
    print("13:Running RFMIX2")                                                                                                                   
    #chromosome=22 ###make it into a loop!                                                                                                        
    prepcmd="rfmix -f "+target_bcf+" -r "+ref_bcf+" -m "+refmap+" -g "+mapf+" -o "+outd+"/"+fname+"_rfmix"+" --chromosome="+str(chrom)+" -e 1 --n-threads=8 --crf-weight=3"
    print(prepcmd)
    os.system(prepcmd)
    filterList= [] 
    fileslist=os.listdir(outd) 
    for f in fileslist: 
        if 'rfmix' in f and str(chrom) in f: 
            print(f) 
            ff=outd+"/"+f 
            filterList.append(ff)
    lname=outd+"/"+"combo"+str(chrom)+".list" 
    with open(lname, mode='wt', encoding='utf-8') as myfile: 
        myfile.write('\n'.join(filterList))
        myfile.write('\n') 
        myfile.close() 
                                                                                                                                          
tmp_reference=os.path.basename(options.reference)
tmp_reference=os.path.basename(options.reference)                                                                                                
tmp_target=os.path.basename(options.target)                                                                                                      
prep_mod1(options.target, options.reference, options.name,options.outdir, options.chrom,options.samples,options.rsamples,options.refmap,options.mapf)

