Installation and Use
====================

**partialsmiles** is a Python library. It can be installed with :command:`pip`::

  pip install -U partialsmiles

The library provides a single API function, ``ParseSmiles``.

.. py:function:: ParseSmiles(smi, partial=False)

   Parse a SMILES string, by default treating it as a complete SMILES string.
   Note that even where ``partial`` is set to ``True``, appending a space or
   a tab causes it to be treated as a complete string.

To test, try the following at the Python prompt.::

        >>> import partialsmiles as ps
        >>> ps.__version__
        1.0
        >>> mol = ps.ParseSmiles("CC(", partial=False)
        Traceback (most recent call last):
        ...
        partialsmiles.exceptions.SMILESSyntaxError: 1 branch has not been closed
          CC(
             ^
        >>> mol = ps.ParseSmiles("CC(", partial=True)
        >>> # no error message
     
You should now read the rest of the documentation (there isn't that much) to avoid any surprises!

If you have any problems or find a bug, please file an issue at https://github.com/baoilleach/partialsmiles/issues.

If this library is used to obtain results for publication, please acknowledge its use, including the version and any modifications made. One way to cite it would be ``N.M. O'Boyle. partialsmiles, version 1.0. Available from https://github.com/baoilleach/partialsmiles``.
