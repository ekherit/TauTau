#!/usr/bin/python3
import os
import re
dir = "./"
count = 1
#re_file_pattern = ".+.cfg.bosslog"
re_file_pattern = ".+.cfg$"
r = re.compile("Event proceed:")
okfile = open("check.txt", "w")
#failfile = open("FAIL.txt", "w")
fmt = "%-50s %50s %50s %10.3f Mb %10s\n" 
for subdir, dirs, files in os.walk(dir):
    for f in files:
        count += 1
        if re.match(re_file_pattern, f):
            fname = os.path.join(subdir, f)
            fname_boss = fname+".bosslog"
            base_fname = os.path.splitext(fname)[0]  
            root_fname = base_fname + ".root"
            #print("Checking %s  %s" % (fname, fname_boss))
            ok = False
            file_size = 0
            try :
                st = os.stat(root_fname)
                file_size =  st.st_size/1024./1024. #in Mb
                file = open(fname_boss,"r");
                for line in file.readlines():
                    if(re.match(r, line)):
                        ok = file_size > 0.01
                        break
                file.close()
            except:
                ok = False
            result = fmt % (fname, fname_boss, root_fname, file_size, "OK" if ok else "FAIL")
            okfile.write(result)
            print(result)

