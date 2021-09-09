import numpy as np
import allel
import getopt, sys, os 
from optparse import OptionParser
import pathlib
import pandas
from zipfile import ZipFile

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
(options, args) = parser.parse_args() 
rpath = pathlib.Path(__file__).resolve().parent
scpath = str(rpath)

#Get all chrs in vcf
def getchr(vcff):
    callset = allel.read_vcf(vcff)
    chrn=np.unique(callset['variants/CHROM']).tolist()
    print('got chr')
    return(chrn)
#call processing on each chr
def perchrcall(chrnum,target,ref,fname,outd,sample,rsample,refmap,mapf,qsubp):
    scmd="python "+scpath+"/LAP_2021.py"+" -t "+target+" -r "+ref+" -n "+fname+" -o "+outd+" -c "+chrnum+" -s "+sample+" -v "+rsample+" -p "+refmap+" -m "+mapf
    jid="JN_"+fname
    qcmd="qsub -cwd -b y -q all.q -P "+qsubp+" -l mem_free=16G -N "+jid+" -e "+outd+"/"+fname+".err"+" -o "+outd+"/"+fname+".out -V "+scmd
    os.system(qcmd)
    print(qcmd)
    return(qcmd)

def getmap(pathm, chrnum):
    print('Get Map')
    inf=pandas.read_table(pathm) 
    print(inf)
    print(chrnum)
    cmap=np.squeeze(inf.loc[inf['chr'].isin([chrnum]),'path'])
    #cmap=(inf.loc[inf['chr']==chrnum]['path'])
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
    scmd="tar -cvf "+outdir+"/"+combo_"+chrn+".tar -T "+lname
    qcmd="qsub -cwd -b y -q all.q -P "+qsubp+" -l mem_free=2G -N Combo"+chrn+" -hold_jid "+jid+" -e "+outdir+"/"+"combo"+chrn+".err"+" -o "+outdir+"/"+"combo"+chrn+".out -V "+scmd
    print(qcmd)
    os.system(qcmd)

chrns=getchr(options.target)

for chrn in chrns:
    chrnum=str(chrn)
    n=options.name+"_"+chrnum
    mapp=getmap(options.mapf, chrn)
    ref=options.reference
    perchrcall(chrnum,options.target, ref, n,options.outdir,options.samples,options.rsamples,options.refmap,mapp, options.qsubp)
    combo(n,options.outdir,chrnum,options.qsubp)
#chrnum=str(22)
#perchrcall(chrnum,options.target, options.reference, options.name,options.outdir,options.samples,options.rsamples,options.refmap,options.mapf)
