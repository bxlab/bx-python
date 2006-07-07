cdef extern from "pwm_utils.h":
    int pattern_match( char* string, char* pattern, int n)

def c_match_consensus( sequence, pattern, size ):
    return pattern_match( sequence, pattern, size )
