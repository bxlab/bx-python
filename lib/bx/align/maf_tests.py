"""
Tests for `bx.align.maf`.
"""

from io import StringIO

import bx.align as align
import bx.align.maf as maf

# A simple MAF from the rat paper days
test_maf = """##maf version=1 scoring=humor.v4
# humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet300/chr1.maf
# /cluster/data/hg15/bed/blastz.rn3/axtNet300/chr1.maf

a score=0.128
s human_hoxa 100  8 + 100257 ACA-TTACT
s horse_hoxa 120  9 -  98892 ACAATTGCT
s fugu_hoxa   88  7  + 90788 ACA--TGCT


a score=0.071
s human_unc 9077 8 + 10998 ACAGTATT
# Comment
s horse_unc 4555 6 -  5099 ACA--ATT
s fugu_unc  4000 4 +  4038 AC----TT
"""

# A more complicated MAF with synteny annotation and such
test_maf_2 = """##maf version=1 scoring=autoMZ.v1
a score=3656.000000
s hg17.chr1                   2005 34 + 245522847 TGTAACTTAATACCACAACCAGGCATAGGGG--AAA-------------
s rheMac2.chr11            9625228 31 + 134511895 TGTAACCTCTTACTGCAACAAGGCACAGGGG------------------
i rheMac2.chr11           C 0 I 1678
s panTro1.chr1                2014 34 + 229575298 TGTAACTTAATACCACAACCAGGCATGGGGG--AAA-------------
i panTro1.chr1            C 0 C 0
s bosTau2.chr5            64972365 47 +  76426644 TCCAGCCATGTGTTGTGATCAG--CCAGGGGCTAAAGCCATGGCGGTAG
i bosTau2.chr5            C 0 I 1462
s canFam2.chr27           45129665 31 +  48908698 TTTGACTCTGTGCTCTTATCAGGCCCAAGGG------------------
i canFam2.chr27           C 0 I 1664
e danRer3.chr18            2360867 428 +  50308305 I
e oryCun1.scaffold_139397      643 1271 -      4771 I
e loxAfr1.scaffold_5603      58454 1915 +     68791 I
e echTel1.scaffold_212365     4641 1430 +      9822 I
e echTel1.scaffold_212365     4641 1430 +      9822 I
e rn3.chr4                29161032 1524 - 187371129 I
e mm7.chr6                28091695 3290 - 149646834 I

"""

# A MAF to test slicing upon
test_maf_3 = """##maf version=1 scoring=none
a score=0
s apple  34 64 + 110 AGGGA---GTTCGTCACT------GTCGTAAGGGTTCAGA--CTGTCTATGTATACACAAGTTGTGTTGCA--ACCG
s orange 19 61 - 100 AGGGATGCGTT--TCACTGCTATCGTCGTA----TTCAGACTTCG-CTATCT------GAGTTGT---GCATTACCG

"""

complex_maf = align.Alignment()
complex_maf.score = "7009"
complex_maf.components.append(
    align.Component(src="human_hoxa", start=100, size=8, strand="+", src_size=100257, text="ACA-TTACT")
)
complex_maf.components.append(
    align.Component(src="horse_hoxa", start=120, size=9, strand="-", src_size=98892, text="ACAATTGCT")
)
complex_maf.components[-1].synteny_left = (maf.MAF_NEW_STATUS, 0)
complex_maf.components[-1].synteny_right = (maf.MAF_CONTIG_STATUS, 0)
complex_maf.components.append(
    align.Component(src="unknown_1", start=150, size=3, strand="-", src_size=98892, text="---ATT---")
)
complex_maf.components.append(
    align.Component(src="unknown_2", start=12, size=1000, strand="+", src_size=1200, text=None)
)
complex_maf.components[-1].empty = True
complex_maf.components[-1].synteny_empty = maf.MAF_INSERT_STATUS
complex_maf.text_size = 9


def test_reader():
    reader = maf.Reader(StringIO(test_maf))
    assert reader.attributes["version"] == "1"
    assert reader.attributes["scoring"] == "humor.v4"

    a = next(reader)
    assert a.score == 0.128
    assert len(a.components) == 3
    check_component(a.components[0], "human_hoxa", 100, 8, "+", 100257, "ACA-TTACT")
    check_component(a.components[1], "horse_hoxa", 120, 9, "-", 98892, "ACAATTGCT")
    check_component(a.components[2], "fugu_hoxa", 88, 7, "+", 90788, "ACA--TGCT")

    a = next(reader)
    assert a.score == 0.071
    assert len(a.components) == 3
    check_component(a.components[0], "human_unc", 9077, 8, "+", 10998, "ACAGTATT")
    check_component(a.components[1], "horse_unc", 4555, 6, "-", 5099, "ACA--ATT")
    check_component(a.components[2], "fugu_unc", 4000, 4, "+", 4038, "AC----TT")

    a = next(reader)
    assert a is None

    reader.close()


def test_writer():
    val = StringIO()
    writer = maf.Writer(val, {"scoring": "foobar"})

    a = align.Alignment()
    a.score = 7009

    a.components.append(
        align.Component(src="human_hoxa", start=100, size=9, strand="+", src_size=1000257, text="ACA-TTACT")
    )
    a.components.append(
        align.Component(src="horse_hoxa", start=120, size=10, strand="-", src_size=98892, text="ACAATTGCT")
    )

    check_component(a.components[0], "human_hoxa", 100, 9, "+", 1000257, "ACA-TTACT")
    check_component(a.components[1], "horse_hoxa", 120, 10, "-", 98892, "ACAATTGCT")

    writer.write(a)

    assert (
        val.getvalue()
        == """##maf version=1 scoring=foobar
a score=7009
s human_hoxa 100  9 + 1000257 ACA-TTACT 
s horse_hoxa 120 10 -   98892 ACAATTGCT 

"""  # noqa: W291
    )


def test_slice():
    b = complex_maf.slice_by_component(0, 101, 105)

    check_component(b.components[0], src="human_hoxa", start=101, size=4, strand="+", src_size=100257, text="CA-TT")
    check_component(b.components[1], src="horse_hoxa", start=121, size=5, strand="-", src_size=98892, text="CAATT")
    check_component(b.components[2], src="unknown_1", start=150, size=3, strand="-", src_size=98892, text="--ATT")
    check_component(b.components[3], src="unknown_2", start=12, size=1000, strand="+", src_size=1200, text=None)
    assert b.components[3].empty
    assert b.components[3].synteny_empty == maf.MAF_INSERT_STATUS

    # test slicing with + strand src
    reader = maf.Reader(StringIO(test_maf_3))
    a = next(reader)
    b = a.slice_by_component(0, 40, 62)
    check_component(
        b.components[0], src="apple", start=40, size=22, strand="+", src_size=110, text="TTCGTCACT------GTCGTAAGGGTTC"
    )
    check_component(
        b.components[1], src="orange", start=28, size=22, strand="-", src_size=100, text="TT--TCACTGCTATCGTCGTA----TTC"
    )

    # test slicing with - strand src
    b = a.slice_by_component(1, 30, 68)
    check_component(
        b.components[0],
        src="apple",
        start=46,
        size=41,
        strand="+",
        src_size=110,
        text="ACT------GTCGTAAGGGTTCAGA--CTGTCTATGTATACACAAGTTG",
    )
    check_component(
        b.components[1],
        src="orange",
        start=32,
        size=38,
        strand="-",
        src_size=100,
        text="ACTGCTATCGTCGTA----TTCAGACTTCG-CTATCT------GAGTTG",
    )

    a = next(reader)
    assert a is None


def test_reverse_complement():
    b = complex_maf.reverse_complement()

    check_component(
        b.components[0], src="human_hoxa", start=100257 - 100 - 8, size=8, strand="-", src_size=100257, text="AGTAA-TGT"
    )
    check_component(
        b.components[1], src="horse_hoxa", start=98892 - 120 - 9, size=9, strand="+", src_size=98892, text="AGCAATTGT"
    )
    assert b.components[1].synteny_right == (maf.MAF_NEW_STATUS, 0)
    assert b.components[1].synteny_left == (maf.MAF_CONTIG_STATUS, 0)
    check_component(
        b.components[2], src="unknown_1", start=98892 - 150 - 3, size=3, strand="+", src_size=98892, text="---AAT---"
    )
    check_component(
        b.components[3], src="unknown_2", start=1200 - 12 - 1000, size=1000, strand="-", src_size=1200, text=None
    )
    assert b.components[3].empty
    assert b.components[3].synteny_empty == maf.MAF_INSERT_STATUS


def test_column_iter():
    expected = [
        ["A", "A", "-"],
        ["C", "C", "-"],
        ["A", "A", "-"],
        ["-", "A", "A"],
        ["T", "T", "T"],
        ["T", "T", "T"],
        ["A", "G", "-"],
        ["C", "C", "-"],
        ["T", "T", "-"],
    ]
    for i, c in enumerate(complex_maf.column_iter()):
        assert c == expected[i]


def test_remove_all_gap_column():
    complex_maf_gap = align.Alignment()
    complex_maf_gap.score = "7009"
    complex_maf_gap.components.append(
        align.Component(src="human_hoxa", start=100, size=8, strand="+", src_size=100257, text="-ACA--TTACT")
    )
    complex_maf_gap.components.append(
        align.Component(src="horse_hoxa", start=120, size=9, strand="-", src_size=98892, text="-ACA-ATTGCT")
    )
    complex_maf_gap.components[-1].synteny_left = (maf.MAF_NEW_STATUS, 0)
    complex_maf_gap.components[-1].synteny_right = (maf.MAF_CONTIG_STATUS, 0)
    complex_maf_gap.components.append(
        align.Component(src="unknown_1", start=150, size=3, strand="-", src_size=98892, text="-----ATT---")
    )
    complex_maf_gap.components.append(
        align.Component(src="unknown_2", start=12, size=1000, strand="+", src_size=1200, text=None)
    )
    complex_maf_gap.components[-1].empty = True
    complex_maf_gap.components[-1].synteny_empty = maf.MAF_INSERT_STATUS
    complex_maf_gap.text_size = 11
    complex_maf_gap.remove_all_gap_columns()
    assert complex_maf_gap == complex_maf


def test_read_with_synteny():
    reader = maf.Reader(StringIO(test_maf_2), parse_e_rows=True)

    a = next(reader)
    check_component(
        a.components[0], "hg17.chr1", 2005, 34, "+", 245522847, "TGTAACTTAATACCACAACCAGGCATAGGGG--AAA-------------"
    )
    check_component(
        a.components[1],
        "rheMac2.chr11",
        9625228,
        31,
        "+",
        134511895,
        "TGTAACCTCTTACTGCAACAAGGCACAGGGG------------------",
    )
    print(a.components[1].synteny_left)
    assert a.components[1].synteny_left == (maf.MAF_CONTIG_STATUS, 0)
    assert a.components[1].synteny_right == (maf.MAF_INSERT_STATUS, 1678)

    rat = a.get_component_by_src_start("rn3.")
    check_component(rat, "rn3.chr4", 29161032, 1524, "-", 187371129, None)
    assert rat.synteny_empty == maf.MAF_INSERT_STATUS


def test_write_with_synteny():
    reader = maf.Reader(StringIO(test_maf_2), parse_e_rows=True)
    a = next(reader)
    val = StringIO()
    writer = maf.Writer(val, {"scoring": "foobar"})
    writer.write(a)
    actual = val.getvalue()
    expected = """##maf version=1 scoring=foobar
a score=3656.0
s hg17.chr1                   2005   34 + 245522847 TGTAACTTAATACCACAACCAGGCATAGGGG--AAA------------- 
s rheMac2.chr11            9625228   31 + 134511895 TGTAACCTCTTACTGCAACAAGGCACAGGGG------------------ 
i rheMac2.chr11                                     C 0 I 1678                                        
s panTro1.chr1                2014   34 + 229575298 TGTAACTTAATACCACAACCAGGCATGGGGG--AAA------------- 
i panTro1.chr1                                      C 0 C 0                                           
s bosTau2.chr5            64972365   47 +  76426644 TCCAGCCATGTGTTGTGATCAG--CCAGGGGCTAAAGCCATGGCGGTAG 
i bosTau2.chr5                                      C 0 I 1462                                        
s canFam2.chr27           45129665   31 +  48908698 TTTGACTCTGTGCTCTTATCAGGCCCAAGGG------------------ 
i canFam2.chr27                                     C 0 I 1664                                        
e danRer3.chr18            2360867  428 +  50308305 I                                                 
e oryCun1.scaffold_139397      643 1271 -      4771 I                                                 
e loxAfr1.scaffold_5603      58454 1915 +     68791 I                                                 
e echTel1.scaffold_212365     4641 1430 +      9822 I                                                 
e echTel1.scaffold_212365     4641 1430 +      9822 I                                                 
e rn3.chr4                29161032 1524 - 187371129 I                                                 
e mm7.chr6                28091695 3290 - 149646834 I                                                 

"""  # noqa: W291
    print(actual)
    print("---")
    print(expected)
    assert actual == expected


def check_component(c, src, start, size, strand, src_size, text):
    assert c.src == src
    assert c.start == start
    assert c.size == size
    assert c.strand == strand
    assert c.src_size == src_size
    assert c.text == text
