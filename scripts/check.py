#!/usr/bin/python3
import os
import re
dir = "./"
count = 1
#re_file_pattern = ".+.cfg.bosslog"
re_file_pattern = ".+.cfg$"
r = re.compile("Event proceed:")
okfile = open("OK.txt", "w")
failfile = open("FAIL.txt", "w")
for subdir, dirs, files in os.walk(dir):
    for f in files:
        count += 1
        if re.match(re_file_pattern, f):
            print(f)
            fname = os.path.join(subdir, f)
            fname_boss = fname+".bosslog"
            print(fname_boss)
            match = False
            try:
                file = open(fname_boss,"r");
                for line in file.readlines():
                    if(re.match(r, line)):
                        match = True
                        okfile.write("%-50s %10s\n" % (fname, "OK"))
                        break
                file.close()
            except:
                match = False
            if not match:
                failfile.write("%-50s %10s\n" % (fname, "FAIL"))
            if(count>100): break
    if(count>100): break

