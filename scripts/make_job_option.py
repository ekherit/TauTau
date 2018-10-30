#!/usr/bin/python
import re,os,sys
import argparse
parser = argparse.ArgumentParser(description='Create tau tau selection configuration files.')
parser.add_argument('dirs',   default="data", nargs="+",help='Input directories where recursivily files will be searched. Last one will be output dir')
parser.add_argument('--output_dir', default='', help='Output dir')
parser.add_argument('--filter', default='.+.dst', help='regexp filter of input file name')
parser.add_argument('--combine', default=r'\d+\.\d+', help='regex template to combine several files into one job')
parser.add_argument('--prefix', default="", help='prefix for output files')
parser.add_argument('--N', type=int,default=1000000000, help='Number of event per job')

cfg = parser.parse_args()

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

#default file list for test algorithms...
test_file_list = [
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
file_list = []

if len(cfg.dirs) == 0:
    print "You must specify input dir with files"
    sys.exit(1)

if cfg.output_dir == '':
    if len(cfg.dirs) == 1:
        cfg.output_dir = os.path.abspath(os.curdir)
        print "Use ./ as output dir"
    else:
        cfg.output_dir = os.path.abspath(cfg.dirs[-1])
        cfg.dirs.pop();

print "Input dirs: ", cfg.dirs
print "Output dir: ", cfg.output_dir

if not os.path.exists(cfg.output_dir):
    os.mkdir(cfg.output_dir)
else:
    c = raw_input("Output dir exists. Owerwrite content? (y/N):")
    if c != 'y': 
        print "exit"
        sys.exit(0)



for dir in cfg.dirs:
    file_list += filter_file_list(create_file_list(dir), cfg.filter)

input_file_dict={}
for f in file_list:
    n = re.findall(cfg.combine, f)
    if len(n) != 0:
        if n[0] in input_file_dict:
            input_file_dict[n[0]].append(f)
        else:
            input_file_dict[n[0]]=[f]


submit_file = open(cfg.output_dir+'/submit.csh',"w");
submit_file2 = open(cfg.output_dir+'/submit.sh',"w");

for W, flist in input_file_dict.items():
    files=""
    for f in flist:
        files=files+'"'+os.path.abspath(f)+'",\n'
    files=files[:-2]
    cfg_file = cfg.output_dir+'/'+cfg.prefix+W+".cfg"
    output_file = cfg.output_dir+'/'+cfg.prefix + W + ".root"
    f = open(cfg_file,'w')
    print "Creating ", cfg_file, "..."
    config = template % (files, cfg.N, output_file, float(W))
    f.write(config)
    submit_file.write("boss.condor "+cfg_file+"\n")
    submit_file2.write("boss.condor "+cfg_file+"\n")

print "To run signle file: boss.exe <filename.cfg>"
print "To run all file: source submit.csh"
