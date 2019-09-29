.. _introduction:

Introduction
============

**partialsmiles** is a Python library that provides a validating SMILES parser with support for partial SMILES. The parser checks syntax, ability to kekulize aromatic system, and checks whether atoms' valences appear on a list of allowed valences.

The main use of the library is to accept or reject SMILES strings while they are being generated character-by-character by a machine-learning method. Previously the only way to check such strings was after the entire string was generated, a process that could take 10s of seconds. Using **partialsmiles**, SMILES strings with invalid syntax or problematic semantics can be detected while they are being generated, allowing them to be discarded or for alternate tokens to be suggested.

Some other toolkits (such as OEChem and the CDK) can read partial SMILES; should you use these instead? Certainly, if they meet your needs. Their goal is to return a molecule object that captures bonding information so far (useful for 'depiction as you type', for example), perhaps converting partial aromatic rings to single bonds and unmatched ring connections to dangling bonds. This is in contrast to the goal of **partialsmiles** which is to work out whether the partial SMILES string could be the prefix of a complete SMILES string that would be accepted by the parser. For example, ``CC(C)(C)(`` may be read as isobutane by other readers, but **partialsmiles** will reject it with a valence error as completion of the string will lead to a five or more valent carbon (see :ref:`valence_errors`).

A secondary goal of this project is to make available a permissively-licensed SMILES parser that correctly reads all aromatic systems in the SMILES reading benchmark. This is with the hope to encourage improvements in existing parsers.

.. note::

    The literature in this area is full of terms such as "valid/invalid SMILES", or worse still "valid/invalid molecules", vague terms which are never defined and have no intrinsic meaning except in the context of syntax. Here I use the terms accepted or rejected. The criteria for acceptance or rejection are described in the section :ref:`errors`. Users may wish to use phrases such as "95% of SMILES were accepted by the partialsmiles library" or "I didn't use the partialsmiles library because it made a fuss about terminology".
