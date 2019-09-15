import sys
import partialsmiles as ps

def help_text():
    print("""validate.py - validate a full or partial SMILES string

Usage: validate.py SMILES

Note that the SMILES string should be enclosed in quotation marks if it
contains characters that your shell would interpret (e.g. asterisks).

The SMILES string will first be interpreted as a full SMILES (1) and any
errors reported. It will then be interpreted as a partial SMILES (2) and
any errors reported.

Finally, (3) each substring of the SMILES string (all starting from the first 
character) will be tested until the first (if any) error is found.

If (1) passes but either of (2) or (3) fails, or (2) passes but (3) fails, please file a bug (you will be prompted to do so).

As an example, consider the following SMILES strings and their
associated failures:

1. "CN(C)" passes tests (1), (2) and (3)
2. "CN(C)(" fails (1), passes (2) and (3) - unmatched parenthesis (fine in partial SMILES as it may be matched later on)
3. "CN(C)(=" fails (1) and (2), passes (3) - a nitrogen with four neighbours causes a valence error (stricly speaking the valence of 5 is accepted, it's the degree of 4 that triggers the error)
4. "CN(C)(=O" fails (1), (2) and (3) - the longest substring is the same as in example 3, and so test (3) fails with the same error
""")    

PASS, FAIL = True, False

if __name__ == "__main__":

    if len(sys.argv) != 2:
        help_text()
        sys.exit(1)

    smi = sys.argv[1]

    print("Input: %s" % smi)


    print("\n1. Parse as full SMILES:\n")
    full_status = PASS
    try:
        ps.ParseSmiles(smi)
    except ps.Error as e:
        print("FAIL\n<<<\n" + str(e) + "\n>>>")
        full_status = FAIL
    except Exception:
        raise
    else:
        print("OK")

    print("\n2. Parse as partial SMILES:\n")
    partial_status = PASS
    try:
        ps.ParseSmiles(smi, partial=True)
    except ps.Error as e:
        print("FAIL\n<<<\n" + str(e) + "\n>>>")
        partial_status = FAIL
    except Exception:
        raise
    else:
        print("OK")

    print("\n3. Parse substrings as partial SMILES:\n")
    substring_status = PASS
    try:
        for i in range(1, len(smi)): 
            ps.ParseSmiles(smi[:i], partial=True)
    except ps.Error as e:
        print("FAIL at %s" % smi[:i])
        print("<<<\n" + str(e) + "\n>>>")
        substring_status = FAIL
    except Exception:
        raise
    else:
        print("OK")

    if (full_status == PASS and (partial_status == FAIL or substring_status == FAIL)) or (partial_status == PASS and substring_status == FAIL):
        print("\n*********************** BUG *************************\n")
        print("Congratulations - you've found a bug in the code\n")
        print("Please file a bug report at:")
        print("   https://github.com/baoilleach/partialsmiles/issues")
        print()
        print("Include the following information:")
        print("  version: %s" % ps.__version__)
        print("  input: %s" % smi)
        print("\n*********************** BUG *************************\n")
