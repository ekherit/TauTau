#!/usr/bin/python
import re,os,sys
import argparse
import ScanPoint
parser = argparse.ArgumentParser(description='Create tau tau selection configuration files.')
parser.add_argument('dirs',   default="data", nargs="+",help='Input directories where recursivily files will be searched. Last one will be output dir')
parser.add_argument('--output_dir', default='', help='Output dir')
parser.add_argument('--filter', default='.+.dst$', help='regexp filter of input file name')
parser.add_argument('--combine', default=r'\d{7}', help='regex template to combine several files into one job')
parser.add_argument('--prefix', default="", help='prefix for output files')
parser.add_argument('--N', type=int,default=-1, help='Number of event per job')
parser.add_argument('--config', default='')
parser.add_argument('--W', help='c.m. energy in filename')


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



for dir in cfg.dirs:
    file_list += filter_file_list(create_file_list(dir), cfg.filter)

print file_list
input_file_dict={}

for f in file_list:
    n = re.findall(cfg.combine, f)
    if len(n) != 0:
        if n[0] in input_file_dict:
            input_file_dict[n[0]].append(f)
        else:
            input_file_dict[n[0]]=[f]
if len(input_file_dict) == 0:
  print "Empty file list. Do nothing."
  sys.exit(0)


print "Output dir: ", cfg.output_dir

if not os.path.exists(cfg.output_dir):
    os.mkdir(cfg.output_dir)
else:
    c = raw_input("Output dir exists. Owerwrite content? (y/N):")
    if c != 'y': 
        print "exit"
        sys.exit(0)

submit_filename = cfg.output_dir+'/submit.csh'
submit_filename2 = cfg.output_dir+'/submit.sh';

submit_file = open(submit_filename,"w");
submit_file2 = open(submit_filename2,"w");

if cfg.config == '':
    for key, flist in input_file_dict.items():
        files=""
        for f in flist:
            files=files+'"'+os.path.abspath(f)+'",\n'
        files=files[:-2]
        W = 1.77686*2
        if cfg.W:
          W=float(key)
        else:
          run = int(key)
          if not (55115 <= run <= 55361) : continue
          if 55115 <= run <= 55155: W=3.539068
          if 55157 <= run <= 55161: W=3.550872
          if 55162 <= run <= 55199: W=3.552865
          if 55200 <= run <= 55231: W=3.553934
          if 55232 <= run <= 55239: W=3.560356
          if 55240 <= run <= 55257: W=3.599572
          if 55347 <= run <= 55361: W=3.601510
        cfg_file = cfg.output_dir+'/'+cfg.prefix+key+".cfg"
        output_file = cfg.output_dir+'/'+cfg.prefix + key + ".root"
        f = open(cfg_file,'w')
        print "Creating ", cfg_file, "..."
        config = template % (files, cfg.N, output_file, W)
        f.write(config)
        submit_file.write("boss.condor "+cfg_file+"\n")
        submit_file2.write("boss.condor "+cfg_file+"\n")
else:
    Scan = ScanPoint.read_scan_point_table(cfg.config)
    for key, flist in input_file_dict.items():
        files=""
        for f in flist:
            files=files+'"'+os.path.realpath(os.path.abspath(f))+'",\n'
        files=files[:-2]
        run = int(key)
        b = filter( lambda point: point.runlist.count(run)>0, Scan)
        if len(b) == 0:  continue
        W = b[0].W
        cfg_file = cfg.output_dir+'/'+cfg.prefix+str(run)+ ".cfg"
        output_file = cfg.output_dir+'/'+cfg.prefix + str(run) + ".root"
        f = open(cfg_file,'w')
        print "Creating ", cfg_file, "..."
        config = template % (files, cfg.N, output_file, W)
        f.write(config)
        submit_file.write("boss.condor "+cfg_file+"\n")
        submit_file2.write("boss.condor "+cfg_file+"\n")

submit_file.close()
submit_file2.close()
os.chmod(submit_filename,0722)
os.chmod(submit_filename2,0722)

print "To run signle file: boss.exe <filename.cfg>"
print "To run all file: source submit.csh"
