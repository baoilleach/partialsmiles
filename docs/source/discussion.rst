Discussion
==========

Testing
-------

How can you be sure that the library does what I say it does?

To begin with there is a test suite, admittedly fairly small, but which provides high code coverage and is useful for catching regressions.

Really, most of the work winkling out bugs is done by:

1. Testing against an existing codebase (Open Babel) to look for cases where there is disagreement in whether a particular complete SMILES is accepted/rejected.

2. Testing the invariant, that *for any complete SMILES string, if any prefix is rejected then any longer prefix (including the complete SMILES) should also be rejected*. The corrollary is that if a complete SMILES is accepted, then any prefix of that should also be accepted.

These are tested using a large number of SMILES strings which I construct on-the-fly from ChEMBL:

1. Read two SMILES strings at a time, A and B

2. Construct new SMILES strings by combining every prefix of A with every suffix of B, and *vice versa*

So how can you be sure that the library does what I say it does? You can't. Not until you run the test scripts, which are all included::

  >>> python suite.py
  ...
  OK
  >>> python test_invariant.py chembl.smi
  0
  100
  ...
  
Alternate tokens
----------------

In the :ref:`introduction`, I mention that one application of the library might be to flag a problem so that an alternate token could be suggested instead. Is this actually feasible? And even if so, should one do it?

...

