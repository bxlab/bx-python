import unittest

import base


class Test(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_select.py ${features}"
    input_features = base.TestFile("""0
                                       0
                                       0
                                       0
                                       0
                                       0
                                       0
                                       1""")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm8_chr7_tiny_no_index.maf")
    output_stdout = base.TestFile("""##maf version=1
a score=8132.0
s mm8.chr7  80083009 147 + 145134094 TAGGGAGGTTGGCATTGGTGCTGGAACTTTCCTTGGCCCCCCAATTTATCGAAGTACTAAGGGTTGGAAGTCTCTGGAGCTGCAGGAGTT--GAGTTTGAGAAAAGGCTCTTGGTGGTTTAAAGAGA----------------GGTTTCAACTGC--------------------------CTCTGGCCTC
s rn4.chr1 136012452 190 + 267910886 TAGGGAGATTGGGATTGGTACTGGAACTTTCCTTGGCCTCCCAGTGTATT-CAGTACTAAGGGTTGGAAGTCTCGGGTGCTACAAGAATTAAGAGTTTGAGAAGAGGCTCTTGGTAGTTTAGAAAGAGAGAAGGACATCTTTGGGTTTCGACTACCTGTGGTGGCAGTGTCAGAATTCAGGCTCTGGCCTC

""")


class TestWithE(base.BaseScriptTest, unittest.TestCase):
    command_line = "./scripts/maf_select.py ${features}"
    input_features = base.TestFile("""0
                                       1
                                       0
                                       0
                                       0
                                       0
                                       0
                                       0""")
    input_stdin = base.TestFile(filename="./test_data/maf_tests/mm10_chr12_lessspe.maf")
    output_stdout = base.TestFile("""##maf version=1
a score=-247111.0
s mm10.chr12                       56694986  40 + 120129022 CTCTTTAG---TCTGG--------TTTTTTAATTTTTTTTTCC-T------CA-----CTGCA
s hetGla2.JH602151                  6912245  40 -   7060556 AATTTCAGCCCCCCCG--------ATGCCTAGGTTTCC---CC-G------CA-----CTGGA
i hetGla2.JH602151                                          C 0 C 0
s micMur1.scaffold_1897              109315  37 +    248211 CGCTTCAG--CTGCAC--------GTGTTTAAATTCCC---GGCA------CA-----TGG--
i micMur1.scaffold_1897                                     I 1 C 0
s tupBel1.scaffold_149545.1-136892     2307  42 +    136892 AGCTTTCG--CTGCGT--------GTGGTTACTTTCTC---CGCA------CACCGCGCAG--
i tupBel1.scaffold_149545.1-136892                          I 1 C 0
s pteVam1.scaffold_182                42258  37 +    455609 AGCTTTAG--CTGCAA--------GTGGTTACATTCTC---TGCA------CA-----CCT--
i pteVam1.scaffold_182                                      I 1 I 2
s eriEur1.scaffold_370206             11243  36 +     65867 ACTTTAAG--TAACTTAAAGAAACTCGGCTACACACTC-------------------------
i eriEur1.scaffold_370206                                   I 1 I 2
s sorAra1.scaffold_233549               581  38 +     65803 GCTTTTAG--CCACTCAAGG----TCGGTTATGTCCAC---TGCA------CA----------
i sorAra1.scaffold_233549                                   I 1 I 2
s loxAfr3.scaffold_9               12283164  24 +  83325590 ------------GCTT--------GTAGTTAAGTTCTC---GGTA------CA----------
i loxAfr3.scaffold_9                                        C 0 C 0
e bosTau7.chr21                    47705583 182 +  69078422 I
""")
