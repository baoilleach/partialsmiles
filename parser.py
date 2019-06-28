import sys
from elements import elements

bondchars = "-=#$\\/"
bondorders = [1, 2, 3, 4, 1, 1]
CLOSE, OPEN = range(2)

def ToBondOrder(bondchar):
    idx = bondchars.index(bondchar)
    return bondorders[idx] # returns 1 for empty string

class Bond:
    def __init__(self, beg, end, bo):
        self.beg = beg
        self.end = end
        self.bo = bo
        self.arom = False
        self.idx = -1

class Atom:
    def __init__(self, element, charge=0):
        self.element = element
        self.charge = charge
        self.implh = -1
        self.arom = False
        self.idx = -1
        self.bonds = []
        self.isotope = 0
    def __repr__(self):
        return "Atom(elem=%d)" % self.element
    def getExplicitDegree(self):
        return len(self.bonds)
    def getExplicitValence(self):
        return sum(bond.bo for bond in self.bonds)

class Molecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
    def addAtom(self, symbol, prev, bondchar):
        element = ToElement(symbol)
        atom = Atom(element)
        atom.arom = symbol[0].islower()
        atom.idx = len(self.atoms)
        self.atoms.append(atom)
        if prev is not None:
            bo = ToBondOrder(bondchar)
            bond = self.addBond(prev, atom.idx, bo)
            if atom.arom and self.atoms[prev].arom and not bondchar:
                bond.arom = True
        return atom
    def addBond(self, beg, end, bo):
        bond = Bond(beg, end, bo)
        bond.idx = len(self.bonds)
        self.bonds.append(bond)
        self.atoms[beg].bonds.append(bond)
        self.atoms[end].bonds.append(bond)
        return bond
    def debug(self):
        for i, atom in enumerate(self.atoms):
            print("Atom: %d %s charge %d implh %d expdeg %d iso %d arom %d" % (i, atom.element, atom.charge, atom.implh, atom.getExplicitDegree(), atom.isotope, atom.arom))
        for i, bond in enumerate(self.bonds):
            print("Bond: %d %d->%d bo %d arom %d" % (i, bond.beg, bond.end, bond.bo, bond.arom))
        for bcsymbol, (atomidx, bondchar) in self.openbonds.items():
            print("Open bond: %d->? bc '%s'" % (atomidx, bondchar))
    def __repr__(self):
        return "Molecule(atoms=%s)" % [str(x) for x in self.atoms]

class SmilesParser:
    def __init__(self):
        pass
    def parse(self, smi):
        self.smi = smi
        self.N = len(smi)
        self.mol = Molecule()
        self.openbonds = {} # dict from symbol -> (atom idx, bondchar)
        self.hcount = []
        self.idx = 0
        self.prev = [None]
        while self.idx < self.N:
            x = smi[self.idx]
            if x == '.':
                self.prev = [None]
                self.idx += 1
            elif x == ')':
                self.prev.pop()
                self.idx += 1
            elif x == '(':
                self.prev.append(self.prev[-1])
                self.idx += 1
            elif x in bondchars:
                self.idx += 1
            elif x.isdigit() or x=='%':
                self.handleBCSymbol()
            elif x in "[*BCNPOSFIbcnos":
                self.parseAtom()
            else:
                raise Exception
        self.setImplicitHydrogenCount()
        self.mol.openbonds = dict(self.openbonds)
        return self.mol

    def setImplicitHydrogenCount(self):
        for i, atom in enumerate(self.mol.atoms):
            hcount = self.hcount[i]
            if hcount > -1:
                atom.implh = hcount
            else:
                explicitvalence = atom.getExplicitValence()
                implicitvalence = SmilesValence(atom.element, explicitvalence)
                implh = implicitvalence - explicitvalence
                if implh > 0 and atom.arom:
                    implh -= 1
                atom.implh = implh

    def handleBCSymbol(self):
        x = self.smi[self.idx]
        bondchar = self.smi[self.idx-1] if self.idx > 0 and self.smi[self.idx-1] in bondchars else ""
        if x == "%":
            bcsymbol = smi[self.idx:self.idx+3]
            self.idx += 2
        else:
            bcsymbol = x
        if bcsymbol in self.openbonds:
            opening, openbc = self.openbonds[bcsymbol]
            openbo = 0 if not openbc else ToBondOrder(openbc)
            closebo = 0 if not bondchar else ToBondOrder(bondchar)
            if closebo and openbo and closebo != openbo:
                raise Exception
            bo = closebo if closebo else openbo if openbo else 1
            self.mol.bonds.append( Bond(opening, self.prev[-1], bo) )
        else:
            self.openbonds[bcsymbol] = (self.prev[-1], bondchar)
        self.idx += 1

    def notAtEnd(self):
        return self.idx < self.N

    def consume(self, char):
        if self.notAtEnd() and self.smi[self.idx] == char:
            self.idx += 1
            return True
        return False

    def parseAtom(self):
        x = self.smi[self.idx]
        bondchar = self.smi[self.idx-1] if self.idx > 0 and self.smi[self.idx-1] in bondchars else ""

        if x == '[':
            self.idx += 1
            if self.notAtEnd() and self.smi[self.idx].isdigit():
                isotope = int(self.smi[self.idx])
                self.idx += 1
                if self.notAtEnd() and self.smi[self.idx].isdigit():
                    isotope = isotope*10 + int(self.smi[self.idx])
                    self.idx += 1
                    if self.notAtEnd() and self.smi[self.idx].isdigit():
                        isotope = isotope*10 + int(self.smi[self.idx])
                        self.idx += 1
            else:
                isotope = 0
            if self.notAtEnd() and self.smi[self.idx].upper() + self.smi[self.idx+1] in elements:
                symbol = self.smi[self.idx:self.idx+2]
                self.idx += 2
            else:
                symbol = self.smi[self.idx]
                self.idx += 1
            self.consume("@")
            self.consume("@")
            hasH = self.consume("H")
            if hasH:
                if self.notAtEnd() and self.smi[self.idx].isdigit():
                    hcount = int(self.smi[self.idx])
                    self.idx += 1
                else:
                    hcount = 1
            else:
                hcount = 0
            hasPlus = self.consume("+")
            if hasPlus:
                if self.notAtEnd() and self.smi[self.idx].isdigit():
                    charge = int(self.smi[self.idx])
                    self.idx += 1
                else:
                    charge = 1
            else:
                hasMinus = self.consume("-")
                if hasMinus:
                    if self.notAtEnd() and self.smi[self.idx].isdigit():
                        charge = -int(self.smi[self.idx])
                        self.idx += 1
                    else:
                        charge = -1
                else:
                    charge = 0
            if self.notAtEnd() and self.smi[self.idx] == ']':
                self.idx += 1
        elif self.notAtEnd() and self.smi[self.idx:self.idx+2] in ["Cl", "Br"]:
            symbol = self.smi[self.idx:self.idx+2]
            self.idx += 2
            hcount = -1
            charge = 0
            isotope = 0
        else:
            self.idx += 1
            symbol = x
            hcount = -1
            charge = 0
            isotope = 0

        atom = self.mol.addAtom(symbol, self.prev[-1], bondchar)
        atom.charge = charge
        atom.isotope = isotope
        self.prev[-1] = atom.idx
        self.hcount.append(hcount)

valencemodel = {
        5: [3], 6: [4], 7: [3, 5], 15: [3, 5],
        8: [2], 16: [2, 4, 6],
        9: [1], 17: [1], 35: [1], 53: [1]
        }
def SmilesValence(elem, val):
    implvalences = valencemodel[elem]
    for implvalence in implvalences:
        if val <= implvalence:
            return implvalence
    return val

def ToElement(symbol):
    if symbol[0].islower():
        symbol = symbol[0].upper() + symbol[1:]
    return elements[symbol]

def ParseSmiles(smi):
    sp = SmilesParser();
    mol = sp.parse(smi)
    return mol

def main():
    text = sys.argv[1]
    mol = ParseSmiles(text)
    mol.debug()

def testChembl():
    with open(r"C:\Tools\LargeData\sortedbylength.smi") as inp:
        for line in inp:
            print(line)
            sp = SmilesParser()
            mol = sp.parse(line.split()[0])

if __name__ == "__main__":
    # testChembl()
    main()
