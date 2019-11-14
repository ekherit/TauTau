#!/bin/python
import re
import copy
# ScanPoint class hold information about scan point: name, energy, energy spread, luminocity and
# runlist
class ScanPoint:
    title=''
    W = 0   # c.m. energy, GeV
    dW = 0  # energy error, GeV
    Sw = 0  # energy spread, MeV
    dSw = 0 # energy spread error, MeV
    L  = 0  # luminosity value, pb^-1
    dL = 0  # luminosity error, pb^-1
    runlist = []  #list of the runs

    def head(self):
        return "%5s %10s %10s %10s %10s %10s %10s %3s %-10s" % ("#name","W,GeV", "dW,GeV", "Sw,GeV", "dSw,GeV", "L,pb^-1", "dL,pb^-1", "", "runs")


#make string with folder runlist
    def fold(self, rl):
        n = len(rl)
        s = ""
        for i in range(0,n):
            if i == 0:  s+=str(rl[i]);continue #begin run
            if i == (n-1):
                if rl[i]==rl[i-1]+1: s+="-"+str(rl[i]); continue
                s+=str(rl[i]);
            if rl[i] != rl[i-1]+1: s+= "-" + str(rl[i-1]) + " "+ str(rl[i])
        return s

#print point ifno
    def str(self):
        s  = "%5s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %3s %-10s" % (self.title, self.W, self.dW, self.Sw, self.dSw, self.L, self.dL, "", self.fold(self.runlist))
        return s

    def print_head(self):
        print self.head()

    def prn(self):
        print self.str()



def read_scan_point_table(cfg_file_name):
    cfg_file = open(cfg_file_name, "r")
    Scan = []
    empty_re = re.compile("^\s*$")
    comment_re = re.compile("^\s*#.*")
    run_range_re=re.compile("^\s*(\d+)\s*-\s*(\d+)")
    single_run_re=re.compile("^\s*\d+")
    for line in cfg_file:
        if re.match(empty_re,line)   : continue
        if re.match(comment_re,line) : continue
        cols = line.split()
        sp = ScanPoint()
        sp.title = cols[0];
        sp.W =   float(cols[1]);
        sp.dW =  float(cols[2]);
        sp.Sw =  float(cols[3]);
        sp.dSw = float(cols[4]);
        sp.L =   float(cols[5]);
        sp.dL =   float(cols[6]);
        runlist = []
        for l in range(7, len(cols)):
            r = re.match(run_range_re,cols[l])
            if r:
                for run in range(int(r.group(1)), int(r.group(2))+1):
                    runlist.append(run)
            r = re.match(single_run_re, cols[l])
            if r:
                runlist.append(int(r.group(0)))
        sp.runlist = sorted(list(set(runlist)))
        Scan.append(sp)
    return Scan

def test():
    Scan = read_scan_point_table("scan_points.txt")
    Scan[0].print_head();
    for sp in Scan:
        sp.prn()




class RunInfo:
    title = ''
    run = 0
    L = 0.0


def read_run_info(filename):
    f = open(filename, "r")
    RI = []
    empty_re = re.compile("^\s*$")
    comment_re = re.compile("^\s*#.*")
    not_number_re = re.compile("^\s*[^\d].*")
    for line in f:
        if re.match(empty_re,line)   : continue
        if re.match(comment_re,line) : continue
        if re.match(not_number_re,line) : continue
        cols = line.split()
        ri = RunInfo()
        #print cols[0], cols[1]
        ri.run = int(cols[0])
        ri.L   = float(cols[1])*1e-3
        RI.append(ri)
    return RI


def calculate_luminosity_from_runinfo(S, RI):
    result = []
    for point in S:
        p = copy.copy(point)
        p.L = 0
        p.dL = point.dL
        for run in point.runlist:
            b = filter(lambda x: x.run == run, RI)
            if len(b) >1:
                print "WARNING: duplicate runinfo for run: ",run
            if len(b) == 1:
                p.L += b[0].L
        result.append(p)
    return result



def make_runinfo_lum(scanpoint_file, runinfo_file):
    Scan = read_scan_point_table(scanpoint_file)
    print "Original file: "
    Scan[0].print_head();
    for sp in Scan:
        sp.prn()
    RI = read_run_info(runinfo_file)
    ScanNew = calculate_luminosity_from_runinfo(Scan, RI)
    print "New file: "
    Scan[0].print_head();
    for sp in ScanNew:
        sp.prn()

    print ""
    print "Difference: "
    for i in range(0,len(Scan)):
        print "%10s %10.6f %10.6f %10.6f %10.6f" % (Scan[i].title, Scan[i].W, Scan[i].L, ScanNew[i].L,  (Scan[i].L-ScanNew[i].L))
    return ScanNew

def write_to_file(Scan, filename):
    f = open(filename,"w")
    if len(Scan)>0:
        f.write(Scan[0].head()+"\n")
    for point in Scan:
        f.write(point.str()+"\n")
    f.close()



#ScanNew = make_runinfo_lum("all_scan_points_ems3.txt", "../share/tau_scan_run_info_55060_55361.txt");
#write_to_file(ScanNew,"new.txt")





