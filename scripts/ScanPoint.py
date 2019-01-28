#!/bin/python
import re

# ScanPoint class hold information about scan point: name, energy, energy spread, luminocity and
# runlist
class ScanPoint:
    title=''
    W = 0   # c.m. energy, GeV
    dW = 0  # energy error, GeV
    Sw = 0  # energy spread, MeV
    dSw = 0 # energy spread error, MeV
    L = 0   # luminocity
    runlist = []  #list of the runs

    def print_head(self):
        s = "%5s %10s %10s %10s %10s %10s %s" % ("#name","W,GeV", "dW", "Sw", "dSw", "L,pb^-1", "   runs")
        print s

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
    def prn(self):
        s = "%5s %10.6f %10.6f %10.3f %10.3f %10.3f" % (self.title, self.W, self.dW, self.Sw, self.dSw, self.L)
        print s,"  ",self.fold(self.runlist)


def read_scan_point_table(cfg_file_name):
    cfg_file = open(cfg_file_name, "r")
    Scan = []
    for line in cfg_file:
        cols = line.split()
        if cols[0] == '#': continue
        sp = ScanPoint()
        sp.title = cols[0];
        sp.W =   float(cols[1]);
        sp.dW =  float(cols[2]);
        sp.Sw =  float(cols[3]);
        sp.dSw = float(cols[4]);
        sp.L =   float(cols[5]);
        run_range_re=re.compile("^\s*(\d+)\s*-\s*(\d+)")
        single_run_re=re.compile("^\s*\d+")
        runlist = []
        for l in range(6, len(cols)):
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

