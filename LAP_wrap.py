import numpy as np
import allel
import getopt, sys, os 
from optparse import OptionParser
import pathlib
import pandas
from zipfile import ZipFile
from os.path import exists

parser = OptionParser() 
parser.add_option("-t", "--target", dest="target",help="Path to target", metavar="FILE") 
parser.add_option("-r", "--reference", dest="reference",help="Path to reference", metavar="FILE") 
parser.add_option("-n", "--name", dest="name",help="Name for output file", metavar="NAME") 
parser.add_option("-o", "--out", dest="outdir",help="output dir", metavar="Out") 
#parser.add_option("-c", "--chr", dest="chrom",help="chromosome", metavar="NUM") 
parser.add_option("-s", "--sample", dest="samples",help="sample target list", metavar="FILE") 
parser.add_option("-v", "--rsample", dest="rsamples",help="sample reference list", metavar="FILE") 
parser.add_option("-p", "--refmap", dest="refmap",help="Reference panel sample population classification map", metavar="FILE") 
parser.add_option("-m", "--map", dest="mapf",help="2 column file with chr and path to corresponding map", metavar="FILE") 
parser.add_option("-q", "--qsubp", dest="qsubp",help="qsub project id", metavar="ID")
parser.add_option("-w", "--crf", dest="crf",help="Assign CRF Weights", metavar="crf")
(options, args) = parser.parse_args() 
rpath = pathlib.Path(__file__).resolve().parent
scpath = str(rpath)
if(options.crf):
    setcrf=options.crf
else:
    setcrf=False

#Get all chrs in vcf
def getchr(vcff):
    callset = allel.read_vcf(vcff)
    chrn=np.unique(callset['variants/CHROM']).tolist()
    print('got chr')
    return(chrn)
#call processing on each chr
def perchrcall(chrnum,target,ref,fname,outd,sample,rsample,refmap,mapf,qsubp):
    if(setcrf!="False"):
        scmd="python "+scpath+"/LAP_2021.py"+" -t "+target+" -r "+ref+" -n "+fname+" -o "+outd+" -c "+chrnum+" -s "+sample+" -v "+rsample+" -p "+refmap+" -m "+mapf+" -w "+str(setcrf)
    else:
        scmd="python "+scpath+"/LAP_2021.py"+" -t "+target+" -r "+ref+" -n "+fname+" -o "+outd+" -c "+chrnum+" -s "+sample+" -v "+rsample+" -p "+refmap+" -m "+mapf
    print(scmd)
    jid="JN_"+fname
    qcmd="qsub -cwd -b y -q all.q -P "+qsubp+" -l mem_free=16G -N "+jid+" -e "+outd+"/"+fname+".err"+" -o "+outd+"/"+fname+".out -V "+scmd
    os.system(qcmd)
    print(qcmd)
    return(qcmd)

def getmap(pathm, chrnum):
    print('Get Map')
    inf=pandas.read_table(pathm) 
    print("Getting map")
    print(inf)
    #print(inf.loc[inf['chr'])
    cmap=np.squeeze(inf.loc[inf['chr'].isin([int(chrnum)]),'path'])
    print(cmap)
    return(cmap)
    
#not using subroutine currently
def getref(pathm, chrnum): 
    print('Get Ref') 
    inf=pandas.read_table(pathm)  
    print(inf) 
    print(chrnum) 
    cmap=np.squeeze(inf.loc[inf['chr'].isin([chrnum]),'refpath']) 
    #cmap=(inf.loc[inf['chr']==chrnum]) 
    print(cmap)  
    return(cmap) 


def split_targetvcf(chrnum,target,outdir,name):
    print("Splitting VCF File")
    print(chrnum)
    #Tabix bgzip file if it does not exist:
    tbf=target+".tbi"
    if(exists(tbf)==False):
        tbcmd="tabix -p vcf "+target
        os.system(tbcmd)
    chrfile=outdir+"/"+name+"_target"+".vcf"
    splitcmd="bcftools view "+target+" --regions "+chrnum+" > "+chrfile
    os.system(splitcmd)
    #sf=chrfile+".sorted.vcf"
    #sortcmd="bcftools sort –temp-dir /local/scratch/achatterjee/ "+chrfile+" > "+sf
    #os.system(sortcmd)
    bgcmd="bgzip -c "+chrfile+" > "+chrfile+".gz"
    os.system(bgcmd)
    bgf=chrfile+".gz"
    print(bgf)
    return(bgf)

def split_refvcf(chrnum,target,outdir,name):
    print("Splitting VCF File")
    print(chrnum)
    #Tabix bgzip file if it does not exist:
    tbf=target+".tbi"
    if(exists(tbf)==False):
        tbcmd="tabix -p vcf "+target
        os.system(tbcmd)
    chrfile=outdir+"/"+name+"_ref"+".vcf"
    splitcmd="bcftools view "+target+" --regions "+chrnum+" > "+chrfile
    #splitcmd="tabix -h "+target+" "+chrnum+" > "+chrfile
    os.system(splitcmd)
    #sf=chrfile+".sorted.vcf"
    #sortcmd="bcftools sort "+chrfile+" > "+sf
    #os.system(sortcmd)
    bgcmd="bgzip -c "+chrfile+" > "+chrfile+".gz"
    os.system(bgcmd)
    bgf=chrfile+".gz"
    tbxcmd="tabix -p vcf "+bgf
    os.system(tbxcmd)
    print(bgf)
    return(bgf)

def combo(fname,outdir,chrn,qsubp):
    jid="JN_"+fname
    filterList= []
    fileslist=os.listdir(outdir)
    for f in fileslist:
        if 'rfmix' in f and chrn in f:
            print(f)
            ff=outdir+"/"+f
            filterList.append(ff)
    lname=outdir+"/"+"combo"+chrn+".list"
    #with open(lname, mode='wt', encoding='utf-8') as myfile:
    #    myfile.write('\n'.join(filterList))
    #    myfile.close()
    scmd="tar -cvf "+outdir+"/"+"combo_"+chrn+".tar -T "+lname
    qcmd="qsub -cwd -b y -q all.q -P "+qsubp+" -l mem_free=2G -N Combo"+chrn+" -hold_jid "+jid+" -e "+outdir+"/"+"combo"+chrn+".err"+" -o "+outdir+"/"+"combo"+chrn+".out -V "+scmd
    print(qcmd)
    os.system(qcmd)

chrns=getchr(options.target)
print(chrns)

for chrn in chrns:
    chrnum=str(chrn)
    n=options.name+"_"+chrnum
    mapp=getmap(options.mapf, chrn)
    print(mapp)
    ref=options.reference
    target=split_targetvcf(chrnum,options.target,options.outdir,n)
    reference=split_refvcf(chrnum,options.reference,options.outdir,n)
    perchrcall(chrnum,target, reference, n,options.outdir,options.samples,options.rsamples,options.refmap,mapp, options.qsubp)
    #combo(n,options.outdir,chrnum,options.qsubp)
#chrnum=str(22)
#perchrcall(chrnum,options.target, options.reference, options.name,options.outdir,options.samples,options.rsamples,options.refmap,options.mapf)
