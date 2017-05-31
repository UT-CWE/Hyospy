import sys
import os

base_path = os.path.dirname(__file__)
sys.path.insert(0, base_path)

lib_path = os.path.abspath(os.path.join(base_path,'DataDownload'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'DataIO'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'GIS'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'SUNTANS'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'UNTRIM'))
sys.path.append(lib_path)

lib_path = os.path.abspath(os.path.join(base_path,'Utils'))
sys.path.append(lib_path)
