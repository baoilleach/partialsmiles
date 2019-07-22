from smiparser import ParseSmiles
from valence import HasCommonValence

def HasUnusualValence(mol):
    for atom in mol.atoms:
        # print(atom, atom.getExplicitValence(), atom.implh)
        if not HasCommonValence(atom):
            return True
    return False
        
def main():
    with open(r"C:\Tools\LargeData\sortedbylength.smi") as inp:
        for line in inp:
            smi = line.split()[0]
            mol = ParseSmiles(smi, partial=False)
            unusual = HasUnusualValence(mol)
            if unusual:
                print(smi)

if __name__ == "__main__":
    main()
    fd
    smi = "Cc1ccccc1"
    smi = "c1ccoc1"
    smi = "Cc1ccccn1"
    mol = ParseSmiles(smi, partial=False)
    unusual = HasUnusualValence(mol)
    print(unusual)
