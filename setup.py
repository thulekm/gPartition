import sys
import os
#print (sys.version)
python_version = sys.version
if python_version.startswith("3"):
	os.system("python setup_v3.py")
else:
	os.system("python setup_v2.py")