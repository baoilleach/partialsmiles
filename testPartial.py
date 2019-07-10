from smiparser import ParseSmiles

def testPartial(smi):
    for i in range(1, len(smi)+1):
        try:
            ParseSmiles(smi[:i], partial=True)
        except:
            print(smi)
            print(smi[:i])
            fd
        ParseSmiles(smi, partial=False)

def testChembl():
    with open(r"C:\Tools\LargeData\sortedbylength.smi") as inp:
        for line in inp:
            smi = line.split()[0]
            testPartial(smi)

if __name__ == "__main__":
    testPartial("[12H]")
    testChembl()
