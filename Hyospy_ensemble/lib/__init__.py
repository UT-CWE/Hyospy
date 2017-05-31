import sys
import os

base_path = os.path.dirname(__file__)
sys.path.insert(0, base_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','DataDownload'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','DataIO'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','GIS'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','SUNTANS'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','UNTRIM'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS','Utils'))
sys.path.append(lib_path)


