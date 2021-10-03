from os import listdir
from os.path import isfile, join
import os
import platform
import sys
import random
import numpy as np
import argparse
import math
import config
import time
np.seterr(divide='ignore', invalid='ignore')

os_type = platform.system()
print("OS Platform: " + os_type+"\n")

text = "Partitioning with SRate"

parser = argparse.ArgumentParser(description = text)  

parser.add_argument("-f", "--filex", help="The Phylip")
parser.add_argument("-o", "--output", help="The output folder")
parser.add_argument("-k", "--kvalue", help="min corelation value")
args = parser.parse_args()

T_Start = time.time()

filex = ""
if args.filex:
	filex = args.filex

cwd = os.getcwd() #get current dir
iqtree_path = config.iqtree_path

output = "output"
if args.output:
	output = args.output
if(not os.path.isdir(output)):
	os.system("mkdir "+output)

if "WINDOW" in os_type.upper():
	if not iqtree_path.endswith("\\") and not iqtree_path.endswith("/"):
		iqtree_path += "\\"
else:
	if not iqtree_path.endswith("/"):
		iqtree_path += "/"

if not output.startswith("/"):
	if "WINDOW" in os_type.upper():
		output = cwd+"\\"+output
	else:
		output = cwd+"/"+output

KVALUE = 0.9995
if args.kvalue:
	KVALUE = float(args.kvalue)
	
def strIntersection(s1, s2):
	out = ""
	for c in s1:
		if c in s2 and not c in out:
			out += c
	return out

def countLines(files):
	if(os.path.isfile(files)):
		with open(files) as f:
			return len(f.readlines())
	else:
		return 0
		
def line_prepender(filename, line):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)

#check aligment file if exists sequences which are gap only
def removeGapsPhylip(filex, dest):
	cfile = open(filex,"r")
	i = 0
	t = 0
	treefn = ""
	if "/" in filex:
		treef = filex.split("/")
		treefn = dest.rstrip("/")+"/"+treef[len(treef)-1].strip()+".RG"
	else:
		treefn = dest.rstrip("/")+"/"+filex+".RG"
	if os.path.isfile(treefn):
		if "WINDOW" in os_type.upper():
			os.system("del "+treefn.replace("/","\\"))
		else: 
			os.system("rm "+treefn)
		
	rfile = open(treefn,"w")
	for line in cfile:
		if i == 0:
			if " " in line:
				tax = int(line.strip().split(" ",1)[0].strip())
				sitex = int(line.strip().split(" ",1)[1].strip())
			elif "\t" in line:
				tax = int(line.strip().split("\t",1)[0].strip())
				sitex = int(line.strip().split("\t",1)[1].strip())
		elif i >= 1:
			if " " in line:
				checkLine = line.split(" ",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
			elif "\t" in line:
				checkLine = line.split("\t",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
		i += 1
	cfile.close()
	rfile.close()
	if t == 0:
		print("The file: "+filex+" is ok.")
		if os.path.isfile(treefn):
			if "WINDOW" in os_type.upper():
				os.system("del "+treefn.replace("/","\\"))
			else:
				os.system("rm "+treefn)
		return 0
		#newfile = filex
	else:
		print("The file: "+filex+" has "+str(t)+" seqs only gaps or missing data.")
		tax = tax - t
		line_prepender(treefn,str(tax)+" "+str(sitex))
		return 1


def modelAccuracy(modelName):
	if "+F{" in modelName:
		tempModel = modelName.split("+F{")
		para = tempModel[1].split("}",1)
		freqs = para[0].split(",")
		total_q = 0.0
		i = 0
		freq = "+F{"
		for par in freqs:
			if i < 3:
				total_q += float(par)
				freq += par+","
				i += 1
		last_par = 1.0 - total_q
		freq += str(last_par)+"}"
		newModel = tempModel[0]+freq
		if len(para) > 1:
			newModel += para[1]
		return newModel
	elif "+FQ+I{" in modelName:
		tempModel = modelName.split("+FQ+I{")
		para = tempModel[1].split("}",1)
		newModel = tempModel[0]+"+FQ+I"
		if len(para) > 1:
			newModel += para[1]
		return newModel
	else:
		return modelName

#copy the content between startText and endText from fromFile
def getInfoFrom(fromFile, startText, endText):
	text = ""
	with open(fromFile) as infile:
		copy = False
		for line in infile:
			if startText in line.strip():
				#text += line.split(startText)[1]
				#text = text.strip()
				#text = text.strip("\n")
				copy = True
			elif endText in line.strip():
				copy = False
			elif copy:
				text += line
	text = text.strip()
	text = text.strip("\n")
	return text

#check site const
def isConst(string1):

	tmpString = "".join(string1)
	stringa = list(set(tmpString.upper()))
	#stringa.remove("-")
	bases = "ACGT"
	ot = strIntersection(stringa,bases.upper())
	if len(ot) <= 1:
		return 1
	else:
		return 0
	#print(ot)

#get model base
def getBase(modelName):
	vBaseModel = modelName
	vBaseModel = modelName.split("+")[0]
	vBaseModel = vBaseModel.split("{")[0].strip()
	return vBaseModel

#get parameter sequence of model
def getParSequence(modelName):
	vreturn = []
	if "{" not in modelName:
		vParCount = modelName.count("+")
		i = 0
		while i < vParCount+1:
			vreturn.append(1.0)
			i += 1
	else:
		modelExt = modelName.split("+")
		i = 0
		while i < len(modelExt):
			if "{" in modelExt[i]:
				ext = modelExt[i].split("{",1)[1]
				ext = ext.split("}")[0]
				pExt = ext.split(",")
				j = 0
				while j < len(pExt):
					vreturn.append(float(pExt[j]))
					j += 1
			else:
				vreturn.append(1.0)
			i += 1
	return vreturn

#get model parameter values
def modelParameter(modelName):
	vreturn = []
	newModel = modelAccuracy(modelName)
	freqs = []
	#get frequences
	if "+F" not in newModel:
		freqs = [0.25,0.25,0.25,0.25]
	elif "+FQ" in newModel:
		freqs = [0.25,0.25,0.25,0.25]
	else:
		par = newModel.split("+F{",1)[1].split("}")[0]
		par = par.split(",")
		total = 1.0
		iz = 0 
		while iz < 3:
			freqs.append(float(par[iz]))
			total = total - float(par[iz])
			iz += 1
		freqs.append(total)
	#declare r: the 4x4 matrix with elements value: 1.0 
	r = []
	q = []
	iz = 0
	while iz < 4:
		r.append([1.0,1.0,1.0,1.0])
		q.append([1.0,1.0,1.0,1.0])
		iz += 1
	baseModel = getBase(modelName)
	if baseModel == "K80" or baseModel == "K2P" or baseModel == "HKY" or baseModel == "HKY85":
		par1 = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par1.split(",")[0])
		r[0][2] = par1
		r[2][0] = par1
		r[1][3] = par1
		r[3][1] = par1
	elif baseModel == "TN" or baseModel == "TN93" or baseModel == "TNe":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		r[0][2] = par1
		r[2][0] = par1
		r[1][3] = par2
		r[3][1] = par2
	elif baseModel == "K81" or baseModel == "K3P" or baseModel == "K81u":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		r[0][2] = par1
		r[2][0] = par1
		r[1][3] = par1
		r[3][1] = par1
		r[0][3] = par2
		r[3][0] = par2
		r[2][1] = par2
		r[1][2] = par2
	elif baseModel == "TPM2" or baseModel == "TPM2u":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		r[0][1] = par1
		r[1][0] = par1
		r[0][3] = par1
		r[3][0] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[2][1] = par2
		r[1][2] = par2
	elif baseModel == "TPM3" or baseModel == "TPM3u":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		r[0][1] = par1
		r[1][0] = par1
		r[1][2] = par1
		r[2][1] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[1][3] = par2
		r[3][1] = par2
	elif baseModel == "TIM" or baseModel == "TIMe":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		par3 = float(par.split(",")[2])
		r[0][2] = par1
		r[2][0] = par1
		r[0][3] = par2
		r[3][0] = par2
		r[1][2] = par2
		r[2][1] = par2
		r[1][3] = par3
		r[3][1] = par3
	elif baseModel == "TIM2" or baseModel == "TIM2e":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		par3 = float(par.split(",")[2])
		r[0][1] = par1
		r[1][0] = par1
		r[0][3] = par1
		r[3][0] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[1][3] = par3
		r[3][1] = par3
	elif baseModel == "TIM3" or baseModel == "TIM3e":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		par3 = float(par.split(",")[2])
		r[0][1] = par1
		r[1][0] = par1
		r[1][2] = par1
		r[2][1] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[1][3] = par3
		r[3][1] = par3
	elif baseModel == "TVM" or baseModel == "TVMe":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		par3 = float(par.split(",")[2])
		par4 = float(par.split(",")[3])
		r[0][1] = par1
		r[1][0] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[1][3] = par2
		r[3][1] = par2
		r[0][3] = par3
		r[3][0] = par3
		r[1][2] = par4
		r[2][1] = par4
	elif baseModel == "SYM" or baseModel == "GTR":
		par = modelName.replace(baseModel,"").split("{")[1].split("}")[0]
		par1 = float(par.split(",")[0])
		par2 = float(par.split(",")[1])
		par3 = float(par.split(",")[2])
		par4 = float(par.split(",")[3])
		par5 = float(par.split(",")[4])
		r[0][1] = par1
		r[1][0] = par1
		r[0][2] = par2
		r[2][0] = par2
		r[0][3] = par3
		r[3][0] = par3
		r[1][2] = par4
		r[2][1] = par4
		r[1][3] = par5
		r[3][1] = par5
	i = 0
	while i < 4:
		qii = 0.0 
		j = 0
		while j < 4:
			if i != j:
				q[i][j]= freqs[i]*r[i][j]
				qii += q[i][j]
			j += 1
		q[i][i] = qii*(-1)
		i += 1
	m = 0.0
	i = 0
	while i < 4:
		m += freqs[i]*q[i][i]
		i += 1
	
	i = 0
	while i < 4:
		j = 0
		while j < 4:
			q[i][j] = q[i][j]/m
			vreturn.append(q[i][j])
			j += 1
		i += 1
	return vreturn

filename = filex
if "/" in filex:
	filename = filex.split("/")[len(filex.split("/"))-1]

def listIntersection(Input1, Input2):
	d = 0
	for index, elem in enumerate(Input2):
		if elem == Input1[index]:
			d+=1
	return d

####GET SEQUENCE
inline = open(filex,"r")
i = 1
tax = 0
numSite = 0
seq = []
for line in inline:
	if i ==1:
		#print(line)
		if(len(line.strip())==0):
			#print "File format's error."
			sys.exit("File format's error.")
		else:
			if " " in line:
				numSite = int(line.split(" ",1)[1].strip())
			elif "\t" in line:
				numSite = int(line.split("\t")[1].strip())
		i+=1
	else:
		if (len(line.strip()) > 0):
			if " " in line:
				if len(line.split(" ", 1)[1].strip()) != numSite:
					# print "File format's error."
					sys.exit("File format's error.")
				else:
					seq.append(line.split(" ", 1)[1].strip())
			elif "\t" in line:
				if len(line.split("\t")[1].strip()) != numSite:
					# saprint "File format's error"
					sys.exit("File format's error")
				else:
					seq.append(line.split("\t")[1].strip())
			tax += 1
		i += 1
inline.close()

site_list = list(map(list, zip(*seq)))
#print("Complete to get the tax and numOfSite:\t"+str(tax)+"\t"+str(numSite)+"\t"+str(len(site_list)))

#SPLIT INVARIANT SITES AND VARIANT SITES
i = 0
siteType = []
while i < numSite:
	siteContent  = ""
	siteContent = site_list[i]
	if isConst(siteContent) == 1:
		siteType.append(1)
	else:
		siteType.append(0)
	i+=1

T_SRstart = time.time()
sRate = []
if not os.path.isfile(output+"/"+filename+".srate"):
	#calculate f(i,j) with f(i,j) = sum(seq[i] == seq[j])
	tIdx = 0
	f = []
	while tIdx < tax - 1:
		sameSite = []
		tI = tIdx + 1
		seqi = list(seq[tIdx].strip())
		while tI < tax:
			sIdx = 0
			total = 0
			seqj = list(seq[tI].strip())
			equalSite = listIntersection(seqi, seqj)

			sameSite.append(equalSite)
			tI += 1
		f.append(sameSite)
		tIdx += 1
	print("Complete to calculate the f(i,j)")

	#total f(i,j)
	i = 0
	total_fij = 0
	while i < tax - 1:
		j = 0
		while j < len(f[i]):
			total_fij += f[i][j]
			j += 1
		i += 1

	#calc site rate
	sIdx = 0
	while sIdx < numSite:
		same_fij = 0
		siteContent = site_list[sIdx]
		totalset = 0

		genGroup = []
		for gen in list(set(siteContent)):
			gIdx = 0
			tempgroup = []
			tempgroup.append(gen)
			while gIdx < len(siteContent):
				if siteContent[gIdx] == gen:
					tempgroup.append(gIdx)
				gIdx += 1
			genGroup.append(tempgroup)

		for gR in genGroup:
			totalset += ((len(gR)-1)*(len(gR)-2))/2 #n*(n-1)/2 remove element 0 is the gene
			ix = 1
			while ix < len(gR)-1:
				yx = ix + 1
				while yx < len(gR):
					same_fij += f[gR[ix]][gR[yx]-gR[ix]-1]
					yx += 1
				ix += 1

		if totalset == 0:
			sRate.append(0.0)
		else:
			sRate.append(float(same_fij))
		if(round(int(sIdx)*100/numSite,1) % 10 == 0.0):
			print("Complete "+str(float(sIdx)*100/numSite)+"%")
		sIdx += 1

	sIdx = 0
	slog = open(output+"/"+filename+".srate","w+")
	while sIdx < len(sRate):
		sRate[sIdx] = float(sRate[sIdx])/total_fij
		slog.write(str(sRate[sIdx])+"\n")
		sIdx += 1
	slog.close()
	print("Complete to calculate sRate:")	
	print(sRate)
else:
	slog = open(output+"/"+filename+".srate","r")
	for sr in slog:
		sRate.append(float(sr.strip()))
	print("Complete to load sRate:")
	print(sRate)
	slog.close()
T_SRfin = time.time()
##################################################
nInv = siteType.count(1) #number of Invariant sites
nVar = siteType.count(0) #number of variant sites
if(2*tax-3 >= 100):
	eps = float(2*tax-3)/nVar
else:
	eps = float(100)/nVar
xeps = 0.0
subsetIdx = []
tempSubset = []
while xeps <= 1:
	i = 0 
	while i < numSite:
		if sRate[i] >= xeps and sRate[i] < xeps+eps and siteType[i] == 0:
			tempSubset.append(i)
		i += 1
	if(len(tempSubset)>=50):
		tempSubset.sort()
		subsetIdx.append(tempSubset)
		tempSubset=[]
	xeps += eps

if(len(tempSubset) > 0):
	tmp = []
	for x in tempSubset:
		tmp.append(x)
	for y in subsetIdx[len(subsetIdx)-1]:
		tmp.append(y)
	tmp.sort()
	subsetIdx[len(subsetIdx)-1] = tmp
	
i = 0
while i < numSite:
	if siteType[i] == 1:
		random.seed()
		ranInx = random.randint(0,len(subsetIdx)-1)
		subsetIdx[ranInx].append(i)
	i += 1

zI = 0
while zI < len(subsetIdx):
	temp = subsetIdx[zI]
	temp.sort() 
	subsetIdx[zI] = temp
	zI += 1

rs = open(output+"/"+filename+".nex","w+")
rs.write("#nexus\nbegin sets;\n");
subsetId = 1
for ss in subsetIdx:
	rs.write("\tcharset Subset"+str(subsetId)+" =")
	for z in ss:
		rs.write(" "+str(z+1))
	rs.write(";\n")
	subsetId += 1
rs.write("end;")	
rs.close()
print("Complete to create the partition by sRate")

#using iqtree to get the best models for each subsets
if not os.path.isfile(output+"/"+filename+".S1.iqtree"):
	command = iqtree_path+"iqtree -s "+filex+" -m MFP -fast -spp "+output+"/"+filename+".nex -safe -pre "+output+"/"+filename+".S1 -seed 0 -redo"
	os.system(command)

#best model list
bmodel = []
if(os.path.isfile(output+"/"+filename+".S1.iqtree")):
	bestModel = getInfoFrom(output+"/"+filename+".S1.iqtree", "Edge-linked-proportional partition model", "MAXIMUM LIKELIHOOD TREE")
	for bm in bestModel.splitlines():
		if len(bm.strip()) > 0:
			if "Model" not in bm:
				if modelAccuracy(bm.strip().split()[3]) not in bmodel:
					bmodel.append(modelAccuracy(bm.strip().split()[3]))

if(len(bmodel)>0):
	if "WINDOW" in os_type.upper():
		os.system("del "+output+"\\"+filename+".S1.*")
	else:
		os.system("rm "+output+"/"+filename+".S1.*")

bm = open(output+"/"+filename+".bestmodel","w+")
bIdx = 1
for bM in bmodel:
	bm.write(str(bIdx)+"\t"+bM+"\n")
	bIdx += 1
bm.write("\n==========================\n")

#MERGE the same base model
mergeModel = []
bIdx = 0 
while bIdx < len(bmodel)-1:
	#ck = 0
	temp_idx = []
	meridx = []
	mIdx = 0
	inMerge = -1
	while mIdx < len(mergeModel):
		x = mergeModel[mIdx]
		meridx = meridx + x
		if x.count(bIdx) >= 1:
			temp_idx = x
			inMerge = mIdx
			mIdx = len(mergeModel)
		mIdx += 1
	
	if inMerge == -1:
		cIdx = bIdx + 1
		if len(temp_idx) == 0:
			temp_idx.append(bIdx)
		while cIdx < len(bmodel):
			if meridx.count(cIdx) == 0:
				list1 = modelParameter(bmodel[bIdx])
				list2 = modelParameter(bmodel[cIdx])
				if list1 == list2:
					corelValue = 1.0
				else:
					corelValue = np.corrcoef(list1, list2)[0, 1]
				print(corelValue)
				if corelValue >= KVALUE:
					temp_idx.append(cIdx)
					bm.write(bmodel[bIdx]+" vs "+bmodel[cIdx]+"\t"+str(corelValue)+"\n")
			cIdx += 1
		if len(temp_idx) > 1:
			if inMerge == -1:
				mergeModel.append(temp_idx)
			else:
				mergeModel[inMerge] = temp_idx
	bIdx += 1

for x in mergeModel:
	bm.write(str(x)+"\n")

bm.write("\n============================\n")
#RELOAD bestModel
newBestModel = []
bIdx = 0
while bIdx < len(bmodel):
	ck = 0
	for x in mergeModel:
		if bIdx in x:
			if bIdx == x[0]:
				newBestModel.append(bmodel[bIdx])
			ck += 1
	
	if ck == 0:
		newBestModel.append(bmodel[bIdx])
	bIdx += 1

bIdx = 1
for bM in newBestModel:
	bm.write(str(bIdx)+"\t"+bM+"\n")
	bIdx += 1
bm.close()

#check and merge small subsets
needRm = []
for mergeData in mergeModel:
	merged = 1
	tmp = subsetIdx[mergeData[0]]
	while merged < len(mergeData):
		for da in subsetIdx[mergeData[merged]]:
			tmp.append(da)
		needRm.append(mergeData[merged])
		merged += 1
	tmp.sort()
	subsetIdx[mergeData[0]] = tmp

ssid = 0
tmpSS = []
while ssid < len(subsetIdx):
	if needRm.count(ssid) == 0:
		tmpSS.append(subsetIdx[ssid])
	ssid += 1

subsetIdx = []
subsetIdx = tmpSS

#re-create partition scheme file
rs = open(output+"/"+filename+".merge.nex","w+")
rs.write("#nexus\nbegin sets;\n");
subsetId = 1
for ss in subsetIdx:
	rs.write("\tcharset Subset"+str(subsetId)+" =")
	for z in ss:
		rs.write(" "+str(z+1))
	rs.write(";\n")
	subsetId += 1
rs.write("end;")	
rs.close()
print("\n================================\nComplete The First Step....Build the 1st Partition file "+output+"/"+filename+".merge.nex\n===================\n")

#adjust partition-id of sites
parId = []
sIdx = 0
while sIdx < numSite:
	parId.append(0)
	sIdx += 1

ssIdx = 0
while ssIdx < len(subsetIdx):
	for x in subsetIdx[ssIdx]:
		parId[int(x)] = ssIdx
	ssIdx += 1

subsetId = 0
ssLh = open(output+"/"+filename+".sitelh","w+")
while subsetId < len(newBestModel):
	if(removeGapsPhylip(filex, output) == 0):
		if not os.path.isfile(output+"/"+filename+".treefile"):
			command = iqtree_path+"iqtree -s "+filex+" -m \""+newBestModel[subsetId].replace("+ASC","")+"\" -fast -te BIONJ -safe -pre "+output+"/"+filename+"."+str(subsetId)+" -wsl -blfix -seed 0 -redo"
		else:
			command = iqtree_path+"iqtree -s "+filex+" -m \""+newBestModel[subsetId].replace("+ASC","")+"\" -fast -te "+output+"/"+filename+".treefile -safe -pre "+output+"/"+filename+"."+str(subsetId)+" -wsl -blfix -seed 0 -redo"
	else:
		if not os.path.isfile(output+"/"+filename+".treefile"):
			command = iqtree_path+"iqtree -s "+output+"/"+filename+".RG -m \""+newBestModel[subsetId].replace("+ASC","")+"\" -fast -te BIONJ -safe -pre "+output+"/"+filename+"."+str(subsetId)+" -wsl -blfix -seed 0 -redo"
		else:
			command = iqtree_path+"iqtree -s "+output+"/"+filename+".RG -m \""+newBestModel[subsetId].replace("+ASC","")+"\" -fast -te "+output+"/"+filename+".treefile -safe -pre "+output+"/"+filename+"."+str(subsetId)+" -wsl -blfix -seed 0 -redo"
	os.system(command)
	if not os.path.isfile(output+"/"+filename+".treefile") and os.path.isfile(output+"/"+filename+"."+str(subsetId)+".treefile"):
		if "WINDOW" in os_type.upper():
			os.system("move "+output+"\\"+filename+"."+str(subsetId)+".treefile "+output+"\\"+filename+".treefile")
		else:
			os.system("mv "+output+"/"+filename+"."+str(subsetId)+".treefile "+output+"/"+filename+".treefile")
	
	ck = 0
	if os.path.isfile(output+"/"+filename+"."+str(subsetId)+".sitelh"):
		getLh = open(output+"/"+filename+"."+str(subsetId)+".sitelh")
		for lh in getLh:
			if "Site_Lh" in lh:
				ssLh.write(str(subsetId)+" "+lh.split(" ",1)[1].strip()+"\n")
				ck += 1
		getLh.close()
		if "WINDOW" in os_type.upper():
			os.system("del "+output+"\\"+filename+"."+str(subsetId)+".*")
		else:
			os.system("rm "+output+"/"+filename+"."+str(subsetId)+".*")
	subsetId += 1
ssLh.close()
T_SLH = time.time()
#SPLIT THE SITE IN SUBSET
if(os.stat(output+"/"+filename+".sitelh").st_size > 0):
	siteID = 0
	sitetemp = open(output + "/" + filename + ".sitelh")
	lhArr = []
	for rl in sitetemp:
		rl = rl.split(" ",1)[1]
		slh = rl.split()
		lhArr.append(slh)
	sitetemp.close()
	nOP = len(lhArr)  # no Of Par

	while siteID < numSite:
		# VARIANT
		if siteType[siteID] == 0:
			maxLh = -10000.0
			bestBM = 0
			bIdx = 0
			while bIdx < nOP:
				if bIdx == 0:
					maxLh = float(lhArr[bIdx][siteID])
				if float(lhArr[bIdx][siteID]) > float(maxLh):
					maxLh = float(lhArr[bIdx][siteID])
					bestBM = bIdx
				bIdx+=1
			parId[int(siteID)] = int(bestBM)
		# INVARIANT BY PROB
		else:
			tmp_value = 0.0
			a_value = []
			tempLh = []
			for bIdx in range(0, nOP - 1):
				tempLh.append(lhArr[bIdx][siteID])
			for xlh in tempLh:
				tmp_value += math.exp(float(xlh))
				a_value.append(tmp_value)
			xz = 0
			while xz < len(a_value):
				a_value[xz] = float(a_value[xz]) / tmp_value
				xz += 1
			random.seed()
			prob = random.uniform(0.0, 1.0)
			xz = 0
			# maxvalue = 0.0
			g_data = 0
			while xz < len(a_value):
				if a_value[xz] >= prob:
					g_data = xz
					xz = len(a_value)
				else:
					xz += 1
			parId[int(siteID)] = int(g_data)
		siteID += 1
else:
	print("Can not split this subset: "+str(subsetId))	
T_cur = time.time()
#MERGE TINY SUBSET INTO NEAREST SUBSET
print("Re-divide partition")
minSubsetCount = 50
needMerge = []
for bSS in list(set(parId)):
	if parId.count(bSS) < minSubsetCount:
		needMerge.append(bSS)

for x in needMerge:
	i = 0
	x_avg = 0.0
	min_distance = 9999999999999.99
	nearIdx = 0
	for sI in parId:
		if sI == x:
			x_avg += sRate[i]
		i += 1
	if parId.count(x) > 0:
		x_avg = float(x_avg)/parId.count(x)
		for bSS in list(set(parId)):
			if bSS not in needMerge:
				y_avg = 0.0
				j = 0
				for sI in parId:
					if sI == bSS:
						y_avg += sRate[j]
					j += 1
				y_avg = float(y_avg)/parId.count(bSS)
				if x_avg > y_avg:
					tmp_dis = x_avg - y_avg
				else:
					tmp_dis = y_avg - x_avg
				if(tmp_dis < min_distance):
					min_distance = tmp_dis
					nearIdx = bSS
		i = 0
		for sI in parId:
			if sI == x:
				parId[i] = nearIdx
			i += 1

T_cur = time.time()
#COMPLETE THE FINAL NEX FILE
rs = open(output+"/"+filename+".FINAL.nex","w+")
rs.write("#nexus\nbegin sets;\n");
sIdx = 1
tmpSubset = list(set(parId))
tmpSubset.sort() 
for bSS in tmpSubset:
	i = 1
	ssId = str(bSS).split("_",1)
	if len(ssId) > 1:
		rs.write("\tcharset Subset"+str(int(ssId[0])+1)+""+ssId[1]+" =")
	else:
		rs.write("\tcharset Subset"+str(int(bSS)+1)+" =")
	while i <= len(parId):
		if(parId[i-1] == bSS):
			rs.write(" "+str(i))
		i+=1
	rs.write(";\n")
	sIdx += 1
rs.write("end;")	
rs.close()
T_Fin = time.time()
print("RUNTIME: Start:\t"+str(T_Start)+"\tSRstart:\t"+str(T_SRstart)+"\tSTfinish:\t"+str(T_SRfin)+"\tLikelihood:\t"+str(T_SLH)+"\tFinish:\t"+str(T_Fin))
print("Mission completed.")
