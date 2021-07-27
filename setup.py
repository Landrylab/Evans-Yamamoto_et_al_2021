import os
from setuptools import setup

PATH = os.path.abspath(".")
requirementPath = PATH + '/requirements.txt'

if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        require_list = list(f.read().splitlines())
       
setup(
    name="BFG",
    version="1.0.0",
    install_requires= require_list,
    extras_require={}
      )
