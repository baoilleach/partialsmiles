# partialsmiles

A validating SMILES parser, with support for partial SMILES (e.g. as generated character by character by a generative model).

It validates syntax and ability to kekulize aromatic systems. It also checks the valence against a set of allowed values (which the user should edit if necessary - see valence.py). It does not check stereochemistry.

This is still in development and not yet ready for production use, but feedback and bug reports welcome.

```
>>> import partialsmiles as ps
>>> mol = ps.ParseSmiles("CC(")
Traceback (most recent call last):
...
partialsmiles.exceptions.SMILESSyntaxError: 1 branch has not been closed
  CC(
     ^
>>> mol = ps.ParseSmiles("CC(", partial=True)
>>> mol.debug()
Atom: 0 6 charge 0 implh 3 expdeg 1 iso 0 arom 0
Atom: 1 6 charge 0 implh 2 expdeg 1 iso 0 arom 0
Bond: 0 0->1 bo 1 arom 0
```

TODO or WIP:
1. Test that shortened versions of valid strings are also valid
2. Test that shortened versions of invalid strings are sometimes invalid
3. Rethink handling of aromatic bonds
