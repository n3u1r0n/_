from nt import *
from poly import *
from helper import *

try:
  from IPython.display import Latex

  def latex(object):
    return Latex('$' + str(object) + '$')
except:
  pass

sieve_all(1_000_000)