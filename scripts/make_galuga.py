#!/usr/bin/python
import re,os,sys


NEVENTS_PER_RUN=10000000
FILE_PREFIX="galuga_"


template = """#include "$ROOTIOROOT/share/jobOptions_ReadRec.txt"
#include "$VERTEXFITROOT/share/jobOptions_VertexDbSvc.txt"
#include "$MAGNETICFIELDROOT/share/MagneticField.txt"
#include "$ABSCORROOT/share/jobOptions_AbsCor.txt"
#include "$TAUTAUROOT/share/jobOptions_TauTau.txt"

// Input REC or DST file name 
EventCnvSvc.digiRootInputFile = {%s};

// Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
MessageSvc.OutputLevel = 5;

// Number of events to be rocessed (default is 10)
ApplicationMgr.EvtMax = %d;

ApplicationMgr.HistogramPersistency = "ROOT";
NTupleSvc.Output = { "FILE1 DATAFILE='%s' OPT='NEW' TYP='ROOT'"};
TauTau.CENTER_MASS_ENERGY=%f;
"""

#print template % ("test.dst",  100, "test.root", 3.53907)

#default file list for test
file_list = [
"e3.53907_EEee_1.dst",
"e3.55087_EEee_1.dst",
"e3.55287_EEee_1.dst",
"e3.55393_EEee_1.dst",
"e3.56036_EEee_1.dst",
"e3.59957_EEee_1.dst",
"e3.60151_EEee_1.dst",
"e3.53907_EEee_2.dst",
"e3.55087_EEee_2.dst",
"e3.55287_EEee_2.dst",
"e3.55393_EEee_2.dst",
"e3.56036_EEee_2.dst",
"e3.59957_EEee_2.dst",
"e3.60151_EEee_2.dst"
]


#create all file list (used for path.walk
def proceed_create_file_list(filelist, directory, files):
    for file in files:
        filelist += [os.path.join(directory,file)]
    filelist.sort()

#create all file list in directory (recursively)
def create_file_list(directory):
    files = []
    os.path.walk(directory, proceed_create_file_list, files)
    return files;

#filter file list with specified regular expression reg
def filter_file_list(files, reg):
    r = re.compile(reg)
    filtered_file_list = []
    for file in files:
        if re.match(r,file):
            filtered_file_list += [file]
    return filtered_file_list

print "Making initial file list"
file_list = filter_file_list(create_file_list("data"), ".*.dst")

#print file_list

input_file_dict={}
for f in file_list:
    n = re.findall(r"[-+]?\d*\.\d+|\d+", f)
    if n[0] in input_file_dict:
        input_file_dict[n[0]].append(f)
    else:
        input_file_dict[n[0]]=[f]


submit_file = open("submit.csh","w");

for W, flist in input_file_dict.items():
    files=""
    for f in flist:
        files=files+'"'+f+'",\n'
    files=files[:-2]
    cfg_file = FILE_PREFIX+W+".cfg"
    output_file = FILE_PREFIX+W+".root"
    f = open(cfg_file,'w')
    print "Creating ", cfg_file, "..."
    config = template % (files,  NEVENTS_PER_RUN, "galuga_"+W+".root", float(W))
    f.write(config)
    submit_file.write("boss.condor "+cfg_file)

print "To run signle file: boss.exe <filename.cfg>"
print "To run all file: source submit.csh"
#    print template % (files,  NEVENTS_PER_RUN, "galuga_"+W+".root", float(W))



