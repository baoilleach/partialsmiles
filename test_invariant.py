import sys
import random
random.seed(1)
import partialsmiles as ps

# We want to find examples where:
# 1. Full parsing passes. There is a substring where it fails.
# 2. Full parsing fails. For a substring where it fails, there is a longer
#    substring where it passes.

INCLUDE_ORIGINAL = False

def two_molecules(fname):
    with open(fname) as inp:
        data = []
        for idx, line in enumerate(inp):
            if idx % 10 == 0:
                print(idx, file=sys.stderr)
            smi = line.split()[0]
            if len(smi) > 30 or len(smi) < 4: continue
            data.append(smi)
            if len(data) == 2:
                yield data
                data = []

def mutants(fname):
    for smiA, smiB in two_molecules(fname):
            if INCLUDE_ORIGINAL:
                yield smiA
                yield smiB
            for i in range(1, len(smiA) -1):
                for j in range(1, len(smiB) - 1):
                    yield smiA[:i] + smiB[j:]
                    yield smiB[:j] + smiA[i:]

if __name__ == "__main__":
    fname = sys.argv[1]
    for smi in mutants(sys.argv[1]):
        try:
            # Does it fail for a substring, but then pass for a
            # longer string? (Or the whole treated as a complete string?)
            substring_failed = False
            for i in range(1, len(smi)+1):
                sub = smi[:i]
                try:
                    ps.ParseSmiles(sub, partial=True)
                except ps.Error:
                    substring_failed = True
                else:
                    if substring_failed:
                        print("ERROR!")
                        print(smi)
                        continue

            try:
                ps.ParseSmiles(smi, partial=False)
            except ps.Error:
                pass
            else:
                if substring_failed:
                    print("ERROR!")
                    print(smi)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            print("Exception!")
            print(smi)
