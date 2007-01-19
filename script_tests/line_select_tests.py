import base
import unittest

class Test( base.BaseScriptTest, unittest.TestCase ):
    command_line = "./scripts/line_select.py ${features}"
    input_features = base.TestFile( """0
                                       1
                                       1
                                       0
                                       1
                                       0""" )
    input_stdin = base.TestFile( """a
                                    b

                                    d
                                    e
                                    f""" )
    output_stdout = base.TestFile( """b

                                      d""" )