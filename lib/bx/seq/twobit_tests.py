import sys
import twobit
import random

def quick_fasta_iter( f ):
    current_header = None
    current_sequence = []
    for line in f:
        if line.startswith( "#" ):
            continue
        if line.startswith( ">" ):
            if current_sequence:
                ## print current_header, "".join( current_sequence )
                yield current_header, "".join( current_sequence )
                current_sequence = []
            current_header = line.strip()[1:]
        else:
            current_sequence.append( "".join( line.split() ) )
    if current_sequence:
        yield current_header, "".join( current_sequence )
        current_sequence = []

def test():
    """
    Nose test generator
    """
    for t in "test", "testN", "testMask":
        test_fa = "test_data/seq_tests/%s.fa" % t
        test_twobit = "test_data/seq_tests/%s.2bit" % t
        yield check_random_subseq_matches, test_fa, test_twobit

def check_random_subseq_matches( test_fa, test_twobit ):
    # Load Fasta data
    expected = {}
    for h, s in quick_fasta_iter( open( test_fa ) ):
        expected[h] = s
    # Open 2bit
    t = twobit.TwoBitFile( open( test_twobit ) )
    for k, s in expected.iteritems():
        assert k in t.index
        # assert t.index[k].size == len(s)
        length = len(s)
        for i in range( 100 ):
            start = random.randint( 0, length-2 )
            end = random.randint( start+1, length )
            assert t[k].get(start, end) == s[start:end]
            assert t[k][start:end] == s[start:end], \
                "seq: %s, start: %d, end: %d\nExpected:\n%s\nActual:\n%s\n" % ( k, start, end, s[start:end], t.get( k, start, end ) )
