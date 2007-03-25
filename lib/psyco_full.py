"""
Attempt to call psyco.full, but ignore any errors.
"""

import sys

try:
    import psyco
    psyco.full()
except:
	pass
    #print >> sys.stderr, "Psyco not found, continuing without it"


