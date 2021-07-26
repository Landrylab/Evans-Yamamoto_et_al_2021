import os

PATH = os.path.abspath(".")
requirementPath = PATH + '/requirements.txt'

install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()
setup(name="BFG", install_requires=install_requires, [...])