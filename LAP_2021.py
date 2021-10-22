#import docker
import getopt, sys, os
from optparse import OptionParser
#import tabix
import os  
parser = OptionParser()
parser.add_option("-t", "--target", dest="target",help="Path to target VCF file", metavar="FILE")
parser.add_option("-r", "--reference", dest="reference",help="Path to reference VCF file", metavar="FILE")
parser.add_option("-n", "--name", dest="name",help="Name for output file", metavar="NAME")
parser.add_option("-o", "--out", dest="outdir",help="output dir", metavar="Out")
parser.add_option("-c", "--chr", dest="chrom",help="chromosome", metavar="NUM")
parser.add_option("-s", "--sample", dest="samples",help="sample target list", metavar="FILE")
parser.add_option("-v", "--rsample", dest="rsamples",help="sample reference list", metavar="FILE")
parser.add_option("-p", "--refmap", dest="refmap",help="Reference panel sample population classification map", metavar="FILE")
parser.add_option("-m", "--map", dest="mapf",help="mapfile", metavar="FILE")
parser.add_option("-w", "--crf", dest="crf",help="Assign CRF Weights", metavar="crf")
(options, args) = parser.parse_args()
 
name=os.path.basename(options.name)

def prep_mod1(target,ref,fname,outd,chrom,sample,rsample,refmap,mapf):
    #Step 1: BGZIP target file
    bgcmd="bgzip -c "+target+" > "+outd+"/"+fname+".vcf.gz"
    print("1: starting bgzip")
    #print(bgcmd)
    #os.system(bgcmd)
    #bgzipf=outd+"/"+fname+".vcf.gz"
    #print(bgzipf)
    #Check if file was created and has values.
    #if(os.stat(bgzipf).st_size == 0):
    #    sys.exit("Could not generate BGZIP file. Failed at Step 1/13")
    #decompress bgzip file
    #bgzipf=target
    #bgcmd="bgzip -d "+bgzipf
    #os.system(bgcmd)
    #bgzipf=outd+"/"+fname+".vcf"
    bgzipf=target
    #Check if decompressed file generated
    #if(os.stat(bgzipf).st_size == 0):
    #    sys.exit("Could not generate decompressed file. Failed at Step 1/13")

    #Step 2: Tabix
    tbxf=bgzipf+".tbi"
    tbx_cmd="tabix -p vcf "+bgzipf
    print("2: running tabix")
    os.system(tbx_cmd)
    #Check if tabix file generated
    if(os.stat(tbxf).st_size == 0):
        sys.exit("Could not generate Tabix file. Failed at Step 2/13")

    #Step 3: bcftools filter
    bcfcmd="bcftools filter --include 'AN=2*N_SAMPLES' -Oz -o "+outd+fname+"filt_out.vcf.gz "+ bgzipf
    print("3: filter with bcftools")
    print(bcfcmd)
    os.system(bcfcmd)
    bcff=outd+"/"+fname+"filt_out.vcf.gz"
    if(os.stat(bcff).st_size == 0):
        sys.exit("Could not generate bcftools filtered file. Failed at Step 3/13")
    print(bcfcmd)

    #Step 4:Conform on filtered file
    print("4: starting conform")
    conformcmd="java -jar conform-gt.24May16.cee.jar ref="+ref+" gt="+bcff+" chrom="+str(chrom)+" match=POS out="+outd+"conform"+fname   
    print(conformcmd)
    os.system(conformcmd)                                                      
    conf=outd+"/conform"+fname+".vcf.gz"
    if(os.stat(conf).st_size == 0):
        sys.exit("Could not generate conform output. Failed at Step 4/13")                                                                                                      
    #Step 5: Tabix on conform output
    tbxcon="tabix -p vcf "+conf
    print("5: Tabix on conformed file")
    os.system(tbxcon)
    tbfile=conf+".tbi"                                                                                          
    if(os.stat(tbfile).st_size == 0):
        sys.exit("Could not generate tabix on conform output. Failed at Step 5/13")

    #Step 6: Merge using bcftools
    mergec="bcftools merge -m none -o "+outd+"/"+fname+"_"+str(chrom)+"_merged.vcf.gz -Oz "+conf+" "+ref
    print("6: merging with bcftools")                                                                                                   
    mergef=outd+"/"+fname+"_"+str(chrom)+"_merged.vcf.gz"
    os.system(mergec)                                                                               
    if(os.stat(mergef).st_size == 0):
        sys.exit("Could not generate merged bcftools output. Failed at Step 6/13")            
                                                                                                       
    #Step 7: Tabix on merged file
    tbxcon="tabix -p vcf "+mergef
    print("7: running tabix on merged file")
    tbfile=mergef+".tbi"
    os.system(tbxcon)
    if(os.stat(tbfile).st_size == 0):
        sys.exit("Could not generate tabix on merged bcftools output. Failed at Step 7/13")
    
    #Step 8: Conform on merged file
    print("8: conform on merged")
    conformcmd="java -jar conform-gt.24May16.cee.jar ref="+ref+" gt="+mergef+" chrom="+str(chrom)+" match=POS out="+outd+"/merge_conform"+fname         
    os.system(conformcmd)
    mcf=outd+"/merge_conform"+fname+".vcf.gz"
    print(mcf)
    if(os.stat(mcf).st_size == 0):                                          
        sys.exit("Could not generate conform on merged bcftools output. Failed at Step 8/13")                                           
    
    #Step 9: Tabix on conformed file
    tbxcon="tabix -p vcf "+mcf                                                
    os.system(tbxcon)
    tbfile=mcf+".tbi"
    print("9: Tabix on merged conformed file")
    if(os.stat(tbfile).st_size == 0):                                                                                      
        sys.exit("Could not generate tabix on conform-gt output. Failed at Step 9/13")
    
    #Step 10: Beagle phasing
    #remove variable, uneccessary
    cmap=mapf                                                     
    print("10: running beagle")
    beagcmd="beagle gt="+mcf+" ref="+ref+" map="+cmap+" out="+outd+"/beagle_"+fname+".phased chrom="+str(chrom)+" impute=False"
    os.system(beagcmd)
    phased=outd+"/beagle_"+fname+".phased.vcf.gz"
    if(os.stat(phased).st_size == 0):
        sys.exit("Could not generate beagle output. Failed at Step 10/13")
    
    #Step 11: Tabix on phased file
    print("11: Tabix on phased")
    tbxcon="tabix -p vcf "+phased
    os.system(tbxcon)
    tbfile=phased+".tbi"
    if(os.stat(tbfile).st_size == 0):
        sys.exit("Could not generate tabix on beagle output. Failed at Step 11/13")

    #Step 12: Generate vcf files using bcftools
    print("12:generating vcf files")                                                                                                  
    prepcmd="bcftools view -S "+sample+" -Ov -o "+outd+"/"+fname+"target_prepped.vcf "+phased
    os.system(prepcmd)
    target_bcf=outd+"/"+fname+"target_prepped.vcf"
    if(os.stat(target_bcf).st_size == 0):
        sys.exit("Could not generate target vcf files. Failed at Step 12/13")
    prepcmd="bcftools view -S "+rsample+" -Ov -o "+outd+"/"+fname+"_ref_prepped.vcf "+ref
    ref_bcf=outd+"/"+fname+"_ref_prepped.vcf"
    os.system(prepcmd)
    if(os.stat(ref_bcf).st_size == 0):
        sys.exit("Could not generate reference vcf files. Failed at Step 12/13")
    
    #Step 13: Run RFMIX2. Testing incomplete due to resources?
    print("13:Running RFMIX2")                                                                                                    
    ###RFMIX Testing incomplete
    #Check to see if user gave crf weight value
    if(options.crf!=False):
        prepcmd="rfmix -f "+target_bcf+" -r "+ref_bcf+" -m "+refmap+" -g "+mapf+" -o "+outd+"/"+fname+"_rfmix"+" --chromosome="+str(chrom)+" -e 1 --n-threads=8 --crf-weight="+str(options.crf)
    else:
        prepcmd="rfmix -f "+target_bcf+" -r "+ref_bcf+" -m "+refmap+" -g "+mapf+" -o "+outd+"/"+fname+"_rfmix"+" --chromosome="+str(chrom)+" -e 1 --n-threads=8"
    print(prepcmd)
    os.system(prepcmd)

    #Make list files for gathering script
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

