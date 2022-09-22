from os import listdir
from os.path import isfile, join
import os
import sys
from string import ascii_uppercase
from sklearn.exceptions import NotFittedError
import random
import numpy as np
import argparse
import math
from sklearn.cluster import KMeans
import pandas as pd
from sklearn import metrics
import config
import time
import datetime
import subprocess

start_date = datetime.datetime.now()


text = "IQTREE partitioning with Kmeans Partitioning"

parser = argparse.ArgumentParser(description = text)  

parser.add_argument("-f", "--filex", help="The Phylip")
parser.add_argument("-o", "--output", help="The output folder")
args = parser.parse_args()

filex = ""
if args.filex:
	filex = args.filex

kpath = "/home/lkthu/kmeans/"
iqtree_path = config.iqtree_path

output = "output"
if args.output:
	output = args.output
if(not os.path.isdir(output)):
	os.system("mkdir "+output)

if(not os.path.isdir(output+"/Result")):
	os.system("mkdir "+output+"/Result")

if(not os.path.isdir(output+"/process")):
	os.system("mkdir "+output+"/process")

path = os.getcwd()
print("Path: "+path)
home = "/home/lkthu"

if not output.startswith("/"):
	output = path + "/" +output

print("Output: "+output)

def logging(fname,content):
	now = datetime.datetime.now() # current date and time
	date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
	wlog = open(fname,"a+")
	wlog.write(content+"\t"+str(date_time)+"\n")
	wlog.close()

filename = filex
if "/" in filex:
	filename = filex.split("/")[len(filex.split("/"))-1]
filePath = filex.replace(filename,"")
os.system("cp "+filex+" "+output+"/process/"+filename)

def strIntersection(s1, s2):
	out = ""
	for c in s1:
		if c in s2 and not c in out:
			out += c
	return out
	
def isConst(string1,prot):
	if prot ==0:
		stringa = list(set(string1.upper()))
		#stringa.remove("-")
		bases = "ACGT"
		ot = strIntersection(stringa,bases.upper())
		if len(ot) <= 1:
			return 1
		else:
			return 0
			#print(ot)
	else:
		stringa = list(set(string1.upper()))
		bases = "ARNDCQEGHILKMFPSTWYV" #string.ascii_uppercase[:26]
		result = ""
		result = strIntersection(stringa,bases.upper())
		if len(result) <= 1:
			return 1
		else:
			return 0

##CHECK DNA or AA #prot=0 DNA
prot = 0	
cdna = open(filex,"r")
zz = 0
for ll in cdna:
	if zz == 1:
		ll = ll.upper()
		if " " in ll:
			if len(ll.split(" ",1)[1].strip()) == ll.split(" ",1)[1].strip().count("-") or len(ll.split(" ",1)[1].strip()) == ll.split(" ",1)[1].strip().count("N"):
				zz = zz-1
			else:
				if len(ll.split(" ",1)[1].strip())*8/10 < (ll.split(" ",1)[1].strip().count("-") + ll.split(" ",1)[1].strip().count("A") + ll.split(" ",1)[1].strip().count("T") + ll.split(" ",1)[1].strip().count("G") + ll.split(" ",1)[1].strip().count("C") + ll.split(" ",1)[1].strip().count("N")):
					prot = 0
				else:
					prot = 1 
		elif "\t" in ll:
			if len(ll.split("\t",1)[1].strip()) == ll.split("\t",1)[1].strip().count("-") or len(ll.split("\t",1)[1].strip()) == ll.split("\t",1)[1].strip().count("N"):
				zz = zz-1
			else:
				if len(ll.split("\t",1)[1].strip())*8/10 < (ll.split("\t",1)[1].strip().count("-") + ll.split("\t",1)[1].strip().count("A") + ll.split("\t",1)[1].strip().count("T") + ll.split("\t",1)[1].strip().count("G") + ll.split("\t",1)[1].strip().count("C") + ll.split("\t",1)[1].strip().count("N")):
					prot = 0
				else:
					prot = 1 
	zz += 1
cdna.close()


#####Add Bsub Function#####
def addJob(command,jobid,nthreads):
	bjFile = "job" + str(jobid) + ".bsub"
	BJ = open(bjFile, 'w')
	BJ.write('#BSUB -J job' + str(jobid) + '\n')
	BJ.write('#BSUB -n ' + str(nthreads) + '\n')
	BJ.write('#BSUB -R "span[hosts=1]"\n')
	BJ.write('#BSUB -q normal\n')
	#BJ.write('#BSUB -m "fit01 fit02 fit03 fit04 fit06 fit07 fit08 fit09 fit10 fit11 fit12 fit13 fit14 fit15 fit16"\n')
	BJ.write('#BSUB -e '+path+'/logs/log' + str(jobid) + '.err\n')
	BJ.write('#BSUB -o '+path+'/logs/log' + str(jobid) + '.out\n')
	BJ.write("\n"+command+"\n")
	BJ.close()
	
	cmd = "bsub < "+bjFile
	#os.system(cmd)
	job = subprocess.check_output(cmd, shell=True)
	jID = str(job).split (">")[0].split ("<")[1]
	return jID

##########CHECK JOB AVAILABLE##################
def checkJob(jobId):
	try:
		cmd = "bjobs -u all | grep "+str(jobId)
		job = subprocess.check_output(cmd, shell=True)
		if str(jobId) in str(job).split("\n"):
			return True
		else:
			return False
		return False
	except subprocess.CalledProcessError:
		return False
	
#########CHECK JOBS AVAIABLE##################
def checkJobs(jobIds):
	for jobid in jobIds:
		if(checkJob(jobid)):
			return True
	return False


###TRUE IF 10 JOBS AVAIL####
def checkFullJob(pcmd):
	try:
		cmd = pcmd
		job = subprocess.check_output(cmd, shell=True)
		lineNum = str(job).split("\n")
		#print(lineNum)
		if len(lineNum) > 15:
			return True
		else:
			return False
		return False
	except subprocess.CalledProcessError:
		return False

def getInfoFrom(fromFile, startText, endText):
	text = ""
	with open(fromFile) as infile:
		copy = False
		for line in infile:
			if startText in line.strip():
				text += line.split(startText)[1]
				text = text.strip()
				text = text.strip("\n")
				copy = True
			elif endText in line.strip():
				copy = False
			elif copy:
				text += line
	text = text.strip()
	text = text.strip("\n")
	return text

####GET SEQUENCE
inline = open(filex,"r")
i = 1
tax = 0
numSite = 0
seq = ""
for line in inline:
	if i ==1:
		if(len(line.strip())==0):
			#print "File format's error."
			sys.exit("File format's error.")
		else:
			if " " in line:
				numSite = int(line.split(" ",1)[1].strip())
			elif "\t" in line:
				numSite = int(line.split("\t")[1].strip())
	else:
		if(len(line.strip()) > 0):
			if " " in line:
				if len(line.split(" ",1)[1].strip()) != numSite:
					#print "File format's error."
					sys.exit("File format's error.")
				else:
					seq += line.split(" ",1)[1].strip()+"\n"
			elif "\t" in line:
				if len(line.split("\t")[1].strip()) != numSite:
					#saprint "File format's error"
					sys.exit("File format's error")
				else:
					seq += line.split("\t")[1].strip()+"\n"
			tax += 1
	i += 1
inline.close()

def getAICc(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Corrected Akaike information criterion (AICc) score:" in line.strip():
				text += line.split(" ")[6]
				return text


modelList = ['LG','JTT','WAG']
lhList = [-9999999.0,-9999999.0,-9999999.0]
slh = []
logging(output+"/"+filename+".log","Start kmeans...")
for model in modelList:
	mi = modelList.index(model)
	command = iqtree_path+"iqtree -s "+filex+" -t BIONJ -m "+model+"+G4+I -wsl -fast -safe -nt AUTO -pre "+output+"/process/"+filename+"."+model+"  -seed 0 -redo \n"
	os.system(command)
	logging(output+"/"+filename+".log","Complete iqtree with model "+model)
	if os.path.isfile(output+"/process/"+filename+"."+model+".sitelh"):
		rh = open(output+"/process/"+filename+"."+model+".sitelh")
		sumlh = 0.0
		slht = []
		for linerh in rh:
			if "Site_Lh" in linerh:
				lh = linerh.strip("Site_Lh").strip().split()
				for lhi in lh:
					sumlh += float(lhi)
					slht.append(float(lhi))
		rh.close()
		slh.append(slht)
		lhList[mi] = sumlh

maxsumLh = lhList[0]
maxidx = 0
mi = 0 
while mi < len(modelList):
	if lhList[mi] >= maxsumLh:
		maxsumLh = lhList[mi]
		maxidx = mi
	mi += 1

logging(output+"/"+filename+".log","Max likelihood model "+modelList[maxidx]+" with lh: "+str(maxsumLh))

#SPLIT INVARIANT SITES AND VARIANT SITES
i = 0
invasites = []
varsites = []
parList = []
parNameList = []
waitParList = []
while i < numSite:
	parList.append('1')
	parNameList.append(filename)
	siteContent  = ""
	for line in seq.splitlines():
		siteContent += line[i]
	siteContent = siteContent.strip()
	if isConst(siteContent,prot) == 1:
		invasites.append(i)
		#print(str(i+1)+" | INVA")
	else:
		#print(str(i+1)+" | VAR")
		varsites.append(i)
	i+=1

logging(output+"/"+filename+".log","Start run SRate...")
waitParList.append(filename)
#CREATE MATRIX WITH TIGER RATES AND ENTROPY RATES FOR EACH SITES
os.system("$PYTHON3/python "+kpath+"genMatrix.py -f "+filex+" -o "+output+" -e 2")
logging(output+"/"+filename+".log","Finish calc SRate...")


sRate = []
eRate = []
if os.path.isfile(output+"/"+filename+".rate"):
	w_m = open(output+"/"+filename+".rate","r")
	var = open(output+"/"+filename+".load.rate","w+")
	#var.write("SITE_ID\tSRATE\tENTRO\n")
	var.write("SITE_ID\tSRATE\n")
	# invar = open(output+"/"+filename+".invar.rate","w+")
	for line in w_m:
		if len(line.strip()) > 0 and "SITE_ID" not in line:
			if "\t" in line:
				ll = line.strip().split("\t")
				#var.write(ll[0].strip()+"\t"+ll[1].strip()+"\t"+ll[2].strip()+"\n")
				var.write(ll[0].strip()+"\t"+ll[1].strip()+"\n")
				sRate.append(ll[1].strip())
				eRate.append(ll[2].strip())
	var.close()
	w_m.close()

logging(output+"/"+filename+".log","Loop Kmeans...")
while len(waitParList) > 0:
	for parIdx in waitParList:
		if not os.path.isfile(output+"/"+parIdx+".load.rate"):
			var = open(output+"/"+parIdx+".load.rate","w+")
			#var.write("SITE_ID\tSRATE\tENTRO\n")
			var.write("SITE_ID\tSRATE\n")
			sI = 0
			nsI = 0
			while sI < len(parNameList):
				if parNameList[sI] == parIdx:
					#var.write(str(nsI+1)+"\t"+sRate[sI].strip()+"\t"+eRate[sI].strip()+"\n")
					var.write(str(nsI+1)+"\t"+sRate[sI].strip()+"\n")
					nsI += 1
				sI += 1
			var.close()
		data = pd.read_csv(output+"/"+parIdx+".load.rate",sep="\t")
		data.head()
		        
		data_x = data.drop(['SITE_ID'], axis=1)
		data_x.head()

		data_site = data_x.to_numpy()
		print(data_site)
		n_clusters_ = 0
		try:
			db = KMeans(n_clusters=2, init='k-means++', n_init=10, max_iter=300, algorithm='auto', random_state=0).fit(data_site)
			labels = db.labels_
			n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
			#print(str(xeps)+ " nCluster: "+str(n_clusters_)+" ")
			#print(str(n_clusters_))
			print(labels)
		except ValueError as err:
			print(str(err))
			n_clusters_ = 1
		logging(output+"/"+filename+".log","Complete kmeans clustering (sklearn) "+parIdx)
		if n_clusters_ == 2:
			if not os.path.isfile(output+"/process/"+parIdx+"_AICc.iqtree"):
				command = ""
				if not os.path.isfile(output+"/"+filename+".treefile"):
					command = iqtree_path+"iqtree -s "+output+"/process/"+parIdx+" -m MF -t BIONJ -nt AUTO -mset JTT,WAG,rtREV,LG -fast -safe -pre "+output+"/process/"+parIdx+"_AICc  -seed 0 -redo \n"
				else:
					command = iqtree_path+"iqtree -s "+output+"/process/"+parIdx+" -m MF -nt AUTO -te "+output+"/"+filename+".treefile -mset JTT,WAG,rtREV,LG -fast -safe -pre "+output+"/process/"+parIdx+"_AICc  -seed 0 -redo \n"
				logging(output+"/"+filename+".log","Run iqtree AICC "+parIdx)
				#command = iqtree_path+"iqtree -s "+output+"/process/"+parIdx+" -m MF -mset JTT,WAG,rtREV,LG -fast -safe -pre "+output+"/process/"+parIdx+"_AICc  -seed 0 -redo \n"

				#print(command)
				os.system(command)
				logging(output+"/"+filename+".log","Complete iqtree AICC "+parIdx)
				if not os.path.isfile(output+"/"+filename+".treefile") and os.path.isfile(output+"/process/"+parIdx+"_AICc.treefile"):
					os.system("cp "+output+"/process/"+parIdx+"_AICc.treefile "+output+"/"+filename+".treefile")

			orgAICc = float(getAICc(output+"/process/"+parIdx+"_AICc.iqtree").strip())
			iClus = 0 
			#par1Num = 0
			#par2Num = 0
			parFile = open(output+"/process/par."+parIdx+"","w+")
			parFile.write("#nexus\nbegin sets;\n")
			while iClus < n_clusters_: 
				parFile.write("\tcharset S"+str(iClus+1)+"=")
				iS = 0
				char = ""
				while iS < len(labels):
					if labels[iS] == iClus:
						char += " "+str(iS+1)
					iS += 1
				parFile.write(char+";\n")
				iClus += 1

			parFile.write("end;")
			parFile.close()

			if not os.path.isfile(output+"/process/"+parIdx+"_aAICc.iqtree"):
				logging(output+"/"+filename+".log","Start iqtree after AICC "+parIdx)
				#command = iqtree_path+"iqtree -s "+output+"/process/"+parIdx+" -spp "+output+"/process/par."+parIdx+" -m MF -mset JTT,WAG,rtREV,LG -fast -redo -safe -pre "+output+"/process/"+parIdx+"_aAICc  -seed 0  \n"
				command = iqtree_path+"iqtree -s "+output+"/process/"+parIdx+" -spp "+output+"/process/par."+parIdx+" -m MFP -blfix -te "+output+"/"+filename+".treefile -mset JTT,WAG,rtREV,LG -fast -redo -safe -pre "+output+"/process/"+parIdx+"_aAICc  -seed 0  \n"
				os.system(command)
				logging(output+"/"+filename+".log","Complete iqtree after AICC "+parIdx)
			transAICc = float(getAICc(output+"/process/"+parIdx+"_aAICc.iqtree").strip())
			if(transAICc <= orgAICc):
				os.system("python "+kpath+"splitPartition.py -f "+output+"/process/"+parIdx+" -p "+output+"/process/par."+parIdx)
				if os.path.isfile(output+"/process/"+parIdx+"S1"):
					if list(labels).count(0) >= 50 and list(labels).count(1) >= 50:
						waitParList.append(parIdx+"S1")
						sIdx = 0
						csIdx = 0
						while sIdx < len(parList):
							if parNameList[sIdx] == parIdx:
								if labels[csIdx] == 0:
									parNameList[sIdx] = parIdx+"S1"
									parList[sIdx] = parList[sIdx]+'1'
								csIdx += 1
							sIdx += 1
				if os.path.isfile(output+"/process/"+parIdx+"S2"):
					if list(labels).count(0) >= 50 and list(labels).count(1) >= 50:
						waitParList.append(parIdx+"S2")
						sIdx = 0
						csIdx = 0
						while sIdx < len(parList):
							if parNameList[sIdx] == parIdx:
								if labels[csIdx] == 1:
									parNameList[sIdx] = parIdx+"S2"
									parList[sIdx] = parList[sIdx]+'2'
								csIdx += 1
							sIdx += 1
		waitParList.remove(parIdx)

#print(parList)
logging(output+"/"+filename+".log","Complete loop kmeans")

parFile = open(output+"/Result/par."+filename+"","w+")
parFile.write("#nexus\nbegin sets;\n")
iClus = 0
parClust = list(set(parList))
while iClus < len(parClust): 
	parFile.write("\tcharset S"+str(iClus+1)+"=")
	iS = 0
	char = ""
	while iS < len(parList):
		if parList[iS] == parClust[iClus]:
			char += " "+str(iS+1)
		iS += 1
	parFile.write(char+";\n")
	iClus += 1
parFile.write("end;")
parFile.close()

logging(output+"/"+filename+".log","Complete create partition "+filename)
