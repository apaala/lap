import numpy as np
import allel
import getopt, sys, os
from optparse import OptionParser
import pathlib
import pandas
from zipfile import ZipFile
import ntpath



parser = OptionParser()
parser.add_option("-o", "--out", dest="outdir",help="output dir", metavar="Out")
parser.add_option("-p", "--path", dest="pdir",help="path with outputfiles to combine", metavar="Path")
parser.add_option("-l", "--prefix", dest="prefix",help="sample output prefix list", metavar="FILE")
parser.add_option("-n", "--name", dest="name",help="Name for output file", metavar="NAME")
parser.add_option("-s", "--separate", dest="sep",help="Per chromosome tar", metavar="NAME", action='store_true')
#parser.add_option("-d", "--process", dest="dest",help="Process the outputs on grid or locally?", metavar="local/grid")
#parser.add_option("-i", "--projectid", dest="pid",help="Project ID for grid submission", metavar="FILE")
(options, args) = parser.parse_args()

rpath = pathlib.Path(__file__).resolve().parent
scpath = str(rpath)


def main():
    outdir=options.outdir
    #prefixlist=pandas.read_table(options.prefix)
    f=open(options.prefix,"r")
    lines=f.readlines()
    result=[]
    for x in lines:
        result.append(x.rstrip('\n'))
    f.close()
    dirpath=options.pdir
    print(result)
    #If separate files for each chromosome
    if(options.sep):
        for p in result:
            print("result")
            print(p)
            filelist=listdir(p,dirpath,p)
    #If combined files requested
    else:
        combinedchr(result, dirpath,outdir,options.name)        
        
#Combine list of files and archive
def combinedchr(results,dirpath,outdir, name):
    all_files=[]
    for r in results:
        file_list = [f for f in os.listdir(dirpath) if f.startswith(r)]
        all_files.extend(file_list)
    print(all_files)
    fullpath=[dirpath + i for i in all_files]
    print(fullpath)
    sis=[f for f in fullpath if f.endswith('sis.tsv')]
    msp=[f for f in fullpath if f.endswith('msp.tsv')]
    fb=[f for f in fullpath if f.endswith('fb.tsv')]
    qf=[f for f in fullpath if f.endswith('.Q')]
    print(sis)
    print(msp)
    print(fb)
    fname=outdir+"/"+name
    concatfiles(fname+'_sis.tsv',sis)
    concatfiles(fname+'_msp.tsv',msp)
    concatfiles(fname+'_fb.tsv',fb)
    final_files=[fname+'_sis.tsv',fname+'_msp.tsv',fname+'_fb.tsv']
    final_files.extend(qf)
    print(final_files)
    zipFileName=name+'.zip'
    zipFilesInDir("", zipFileName, final_files)

#concat all files
def concatfiles(fname,flist):
    with open(fname, 'w') as outfile:
        for f in flist:
            with open(f) as infile:
                outfile.write(infile.read())
    outfile.close()
    infile.close()

#Get list of files
def listdir(prefixlist,dirpath,name):
    print(prefixlist,dirpath,name)
    file_list = [f for f in os.listdir(dirpath) if f.startswith(prefixlist)]
    print(file_list)
    zipFilesInDir(dirpath, str(name)+'.zip', file_list)
    return(file_list)

# Zip the files from given directory that matches the filter
def zipFilesInDir(dirName, zipFileName, filenames):
    # create a ZipFile object
    with ZipFile(zipFileName, 'w') as zipObj:
        # Iterate over all the files in directory
        for f in filenames:
            filePath = os.path.join(dirName, f)
            # Add file to zip
            zipObj.write(filePath, os.path.basename(filePath))


if __name__ == '__main__':
    main()
