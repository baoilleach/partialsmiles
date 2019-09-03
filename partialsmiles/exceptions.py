# MIT License
# 
# Copyright (c) 2019 Noel O'Boyle
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

class Error(Exception):
    """Base class for partialsmiles exceptions."""
    def __init__(self, message, smi, idx):
        self.smi = smi
        self.idx = idx
        self.message = message
    def __str__(self):
        error = [self.message, "  " + self.smi, " "*(self.idx + 2)+"^"]
        return "\n".join(error)

class SMILESSyntaxError(Error):
    """Exception raised for syntax errors

    Attributes:
        message -- explanation of the error
        smi -- input SMILES in which the error occurred
        idx -- index of the character in the SMILES when the error occurred
    """

class KekulizationFailure(Error):
    """Exception raised when an aromatic system cannot be kekulized

    Attributes:
        message -- explanation of the error
        smi -- input SMILES in which the error occurred
        idx -- index of the character in the SMILES when the error occurred
    """

class ValenceError(Error):
    """Exception raised when the valence is not consistent with allowed values

    Attributes:
        message -- explanation of the error
        smi -- input SMILES in which the error occurred
        idx -- index of the character in the SMILES when the error occurred
    """
