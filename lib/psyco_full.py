"""
Attempt to call psyco.full, but ignore any errors.
"""

try:
    import psyco
    psyco.full()
except Exception:
    pass
