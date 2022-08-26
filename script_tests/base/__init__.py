import filecmp
import os
import string
import subprocess
import tempfile
import unittest
from io import StringIO


class TestFile:
    def __init__(self, text=None, filename=None):
        assert text is None or filename is None, "Cannot specify both text and filename for input"
        self.text = text
        self.filename = filename
        if self.filename is None:
            _, tf_name = tempfile.mkstemp()
            tf = open(tf_name, "w")
            sio = StringIO(self.text)
            for line in sio:
                print(line.lstrip(" ").rstrip("\r\n"), file=tf)
            tf.close()
            self.tempfile = True
            self.filename = tf_name
        else:
            self.tempfile = False

    def check(self, other_fname):
        assert filecmp.cmp(self.filename, other_fname), f"Files do not match ({self.filename}, {other_fname})"

    def __del__(self):
        if self.tempfile:
            os.remove(self.filename)


class BaseScriptTest:
    """
    Helper class for testing a command line tool
    """

    def test_script(self):
        # Accumulate parameters
        input_files = dict()
        output_files = dict()
        out_dir = None
        stdin = stdout = stderr = None
        for key in dir(self):
            if key == "command_line":
                command_line = getattr(self, key)
            elif key.startswith("input_"):
                value = getattr(self, key)
                assert isinstance(value, TestFile)
                arg_name = key[6:]
                input_files[arg_name] = value
            elif key.startswith("output_"):
                value = getattr(self, key)
                assert isinstance(value, TestFile)
                arg_name = key[7:]
                output_files[arg_name] = value
            elif key == "out_dir":
                out_dir = getattr(self, key)
                assert os.path.isdir(out_dir)
        # Build the command line
        input_fnames = dict()
        output_fnames = dict()
        all_fnames = dict()
        for key, value in input_files.items():
            input_fnames[key] = value.filename
            all_fnames[key] = input_fnames[key]
            if key == "stdin":
                stdin = open(input_fnames[key])
        for key in output_files.keys():
            _, tf_name = tempfile.mkstemp()
            output_fnames[key] = tf_name
            all_fnames[key] = output_fnames[key]
            if key == "stdout":
                stdout = open(output_fnames[key], "w")
                stdout.flush()
            if key == "stderr":
                stderr = open(output_fnames[key], "w")
                stdout.flush()
        if out_dir is not None:
            temp_out_dir = tempfile.mkdtemp()
            all_fnames["out_dir"] = temp_out_dir
            for root, _, files in os.walk(out_dir):
                for file in files:
                    output_files[os.path.join(root, file)] = TestFile(filename=os.path.join(root, file))
                    output_fnames[os.path.join(root, file)] = os.path.join(temp_out_dir, file)
        real_command = string.Template(command_line).substitute(all_fnames)
        # Augment PYTHONPATH, bit of a HACK here! need to suck this data from setuptools or something?
        env = dict(os.environ)
        if "PYTHONPATH" in env:
            env["PYTHONPATH"] = "./lib:" + env["PYTHONPATH"]
        else:
            env["PYTHONPATH"] = "./lib"
        # Run the command
        subprocess.check_call(real_command, stdin=stdin, stdout=stdout, stderr=stderr, shell=True, env=env)
        # Check the outputs
        for key, value in output_files.items():
            value.check(output_fnames[key])
        # Cleanup
        for value in output_fnames.values():
            os.remove(value)


class TestTest(BaseScriptTest, unittest.TestCase):
    input_in1 = TestFile("""Foo\nBar\nBaz""")
    output_stdout = TestFile("""Foo""")
    command_line = "/usr/bin/head -1 ${in1}"


class TestTest2(BaseScriptTest, unittest.TestCase):
    input_in1 = TestFile("/etc/passwd")
    output_stdout = TestFile("/etc/passwd")
    command_line = "cat ${in1}"


if __name__ == "__main__":
    unittest.main()
