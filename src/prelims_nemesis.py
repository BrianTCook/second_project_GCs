from amuse.lab import*
from nemesis import Nemesis, HierarchicalParticles
from main import parent_worker, sub_worker, py_worker

Mgal, Rgal = 5e6|units.MSun, 10.|units.parsec

nemesis = Nemesis(parent_worker, sub_worker, py_worker)

print('hello world!')
