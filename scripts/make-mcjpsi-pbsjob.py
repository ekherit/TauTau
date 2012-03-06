#!/usr/bin/python
import os
import string
import sys
import fileinput
import re

ENERGY=1843
RANDOM_SEED=0
EVENT_NUMBER=1
JOB_NAME="test"
#PBS_QUEUE="publicq@torqsrv"
PBS_QUEUE="besq@torqsrv"
ENERGY_SPREAD=1.12

def configure(source_file_name,  target_file_name,  JOB_NAME,  ENERGY,  RANDOM_SEED,  EVENT_NUMBER,  RUN_NUMBER):
#Open source template files with simulation configuration
	source_file = open(source_file_name, 'r')
	target_file = open(target_file_name, 'w')
	for line in source_file:
		line = re.sub("TEMPLATE_NAME", JOB_NAME, line)
		line = re.sub("TEMPLATE_BEAM_ENERGY_SPREAD", str(ENERGY_SPREAD/1.0e3), line)
		line = re.sub("TEMPLATE_BEAM_ENERGY", str(2*ENERGY/1e3), line)
		line = re.sub("TEMPLATE_RANDOM_SEED", str(RANDOM_SEED), line)
		line = re.sub("TEMPLATE_EVENT_NUMBER", str(EVENT_NUMBER), line)
		line = re.sub("TEMPLATE_PBS_QUEUE", PBS_QUEUE, line)
		line = re.sub("TEMPLATE_RUN_NUMBER", RUN_NUMBER, line)
		target_file.write(line)
	source_file.close()
	target_file.close()

def do_mc(energy,  event_number,  job_number,  run_number):
	random_seed=int(energy*100)/10+job_number
	name="mcpsip_"+str(energy)+"_"+str(job_number)
	template_dir="/bes3fs/groups/tauqcd/tauqcdgroup/nikolaev/mc/mcpsip-template"
	work_dir = name
	if os.path.exists(work_dir):
		print "Directory "+work_dir+" exist! Do nothing!!"
		return
	else:
		print "Creating directory "+work_dir
		os.mkdir(work_dir)
	JOB_NAME=name
	RANDOM_SEED=random_seed
	EVENT_NUMBER=event_number
	ENERGY=energy
	RUN_NUMBER=run_number

#Do substitution for required value
	configure(template_dir+"/template_sim.cfg", work_dir+"/"+name+"_sim.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/template_rec.cfg", work_dir+"/"+name+"_rec.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/template_ana.cfg", work_dir+"/"+name+"_ana.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/common_sim.cfg", work_dir+"/common_sim.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/pbsjobs", work_dir+"/pbsjobs", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)

	print "Starting the "+PBS_QUEUE+" job " + name
	os.system("qsub "+work_dir+"/pbsjobs")
#end of do_mc function

N=50000
jobs = 20
for job in range(1,jobs+1):
	do_mc(1777, N, job,  "-20334, -20335, -20339, -20336")
