#Author: Thulek@gmail.com
#Setting up 
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time

os_type = platform.system()

if os.path.isfile("config.py"):
	print("File config.py is exists.\n")
	rconfig = open("config.py","r")
	print("==================CONFIG.PY FILE=================")
	for line in rconfig:
		print(line.strip()+"\n")
	rconfig.close()
	reinstall = str(input("Do you want to overwrite? [Y/N]: "))
	if reinstall == "Y" or reinstall == "y":
		if "WINDOW" in os_type.upper():
			os.system("del config.py")
		else:
			os.system("rm config.py")
		rconfig = open("config.py","w")
		#iqtree_path = "$IQTREE/" #IQ-TREE path including splash
		#tiger_path = "/home/lkthu/tiger_original/"
		iq_path = input("Input the IQ-TREE path including splash: [$IQTREE]:")
		rconfig.write("iqtree_path = \""+str(iq_path)+"\" \n")
		rconfig.close()
		if os.path.isfile("config.py"):
			print("Setup completed.")
		else:
			print("Setup has an error.")
else:
	rconfig = open("config.py","w")
	#iqtree_path = "$IQTREE/" #IQ-TREE path including splash
	#tiger_path = "/home/lkthu/tiger_original/"
	iq_path = input("Input the IQ-TREE path including splash: [$IQTREE]:")
	#tiger_path = raw_input("Input the TIGER path including splash: [tiger_original/]:")
	rconfig.write("iqtree_path = \""+str(iq_path)+"\" \n")
	#rconfig.write("tiger_path = \""+str(tiger_path)+"\" \n")
	rconfig.close()
	if os.path.isfile("config.py"):
		print("Setup completed.")
	else:
		print("Setup has an error.")
		