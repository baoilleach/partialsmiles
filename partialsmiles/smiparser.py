import sys
from .elements import elements
from . import valence
from . import kekulize
from .exceptions import *

bondchars = "-=#$\\/:"
bondorders = [1, 2, 3, 4, 1, 1, 1]
CLOSE, OPEN = range(2)

aromatic_list = ["c", "n", "p", "o", "s", "te", "se", "b"] # more?
aromatic_set = set(aromatic_list)
firstLetterOfElements = set(x[0] for x in list(elements.keys()) + aromatic_list)

def ToBondOrder(bondchar):
    idx = bondchars.index(bondchar)
    return bondorders[idx] # returns 1 for empty string

class Bond:
    def __init__(self, beg, end, order):
        self.beg = beg
        self.end = end
        self.order = order
        self.arom = False
        self.idx = -1

    def getNbr(self, atom):
        return self.beg if self.end==atom else self.end

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
        return "Atom(elem=%d,chg=%d,idx=%d)" % (self.element, self.charge, self.idx)

    def getExplicitDegree(self):
        return len(self.bonds)

    def getExplicitValence(self):
        return sum(bond.order for bond in self.bonds)

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
            bond = self.addBond(prev, atom, bo)
            if bondchar == ':' or (atom.arom and prev.arom and not bondchar):
                bond.arom = True
        return atom

    def addBond(self, beg, end, bo):
        bond = Bond(beg, end, bo)
        bond.idx = len(self.bonds)
        self.bonds.append(bond)
        beg.bonds.append(bond)
        end.bonds.append(bond)
        return bond

    def getBond(self, atomA, atomB):
        for nbr_bond in atomA.bonds:
            nbr = nbr_bond.end if nbr_bond.beg == atomA else nbr_bond.beg
            if nbr == atomB:
                return nbr_bond
        return None

    def debug(self):
        for i, atom in enumerate(self.atoms):
            print("Atom: %d %s charge %d implh %d expdeg %d iso %d arom %d" % (i, atom.element, atom.charge, atom.implh, atom.getExplicitDegree(), atom.isotope, atom.arom))
        for i, bond in enumerate(self.bonds):
            print("Bond: %d %d->%d bo %d arom %d" % (i, bond.beg.idx, bond.end.idx, bond.order, bond.arom))
        for bcsymbol, (atom, bondchar) in self.openbonds.items():
            print("Open bond: %d->? bc '%s'" % (atom.idx, bondchar))

    def __repr__(self):
        return "Molecule(atoms=%s)" % [str(x) for x in self.atoms]

class SmilesParser:

    def __init__(self, partial, rulesToIgnore):
        self.partial = partial
        self.rulesToIgnore = rulesToIgnore

    def parse(self, smi):
        self.smi = smi
        self.N = len(smi)
        self.mol = Molecule()
        self.openbonds = {} # dict from symbol -> (atom idx, bondchar)
        self.hcount = []
        self.idx = 0
        self.prev = [None]
        self.smiidx = []
        self.reaction_part = 0
        self.bondchar = None
        self.parsing_atom = False
        while self.idx < self.N:
            x = smi[self.idx]
            if x in "CcONon[BPSFIbps*":
                self.handleError(SMILESSyntaxError, self.parseAtom())
            elif x in '. \t>':
                self.handleComponent()

                self.prev[-1] = None
                if x == '.':
                    self.handleError(SMILESSyntaxError, self.validateSyntax(dot=True))
                    self.idx += 1
                elif x == '>':
                    self.reaction_part += 1
                    if self.reaction_part == 3:
                        self.handleError(SMILESSyntaxError, "Reactions only have three parts")
                    self.handleError(SMILESSyntaxError, self.validateSyntax(dot=True))
                    self.idx += 1
                else: # regard whitespace as the stop token
                    self.partial = False
                    break
            elif x == ')':
                if not self.rulesToIgnore & 2 and (self.idx > 1 and smi[self.idx-1]=='('):
                    self.handleError(SMILESSyntaxError, "Empty branches are not allowed")
                if not self.rulesToIgnore & 4 and self.idx > 1 and smi[self.idx-1]==')':
                    self.handleError(SMILESSyntaxError, ("The final branch should not be within parentheses", self.idx-1))
                if self.bondchar:
                    self.handleError(SMILESSyntaxError, "An atom must follow a bond symbol")
                self.prev.pop()
                if not self.prev:
                    self.handleError(SMILESSyntaxError, "Unmatched close parenthesis")
                self.idx += 1
            elif x == '(':
                if self.prev[-1] is None or smi[self.idx-1]=='(':
                    self.handleError(SMILESSyntaxError, "An atom must precede an open parenthesis")
                if self.bondchar:
                    self.handleError(SMILESSyntaxError, "A bond symbol should not precede an open parenthesis")
                self.prev.append(self.prev[-1])
                self.idx += 1
            elif x in bondchars:
                if self.prev[-1] is None:
                    self.handleError(SMILESSyntaxError, "An atom must precede a bond symbol")
                if self.bondchar:
                    self.handleError(SMILESSyntaxError, "Only a single bond symbol should be used")
                if not self.rulesToIgnore & 64 and x == ':':
                    self.handleError(SMILESSyntaxError, "Aromatic bond symbols are rejected by default")
                self.bondchar = x
                self.idx += 1
            elif x.isdigit() or x=='%':
                if self.prev[-1] is None:
                    self.handleError(SMILESSyntaxError, "An atom must precede a bond closure symbol")
                if not self.rulesToIgnore & 32:
                    precedingtext = self.smi[self.smiidx[self.prev[-1].idx]+1:self.idx]
                    if ")" in precedingtext:
                        self.handleError(SMILESSyntaxError, "Ring closure symbols must immediately follow an atom")
                    if "(" in precedingtext:
                        self.handleError(SMILESSyntaxError, "Ring closure symbols should not be in parentheses")
                self.handleError(SMILESSyntaxError, self.handleBCSymbol())
            else:
                self.handleError(SMILESSyntaxError, "Illegal character")

        if not self.partial:
            self.handleComponent()

        self.handleError(SMILESSyntaxError, self.validateSyntax())

        self.setImplicitHydrogenCount()
        self.incompleteAtoms = set(x[0] for x in self.openbonds.values())
        if self.partial:
            self.incompleteAtoms.update(x for x in self.prev if x is not None)

        self.handleError(ValenceError, self.validateValence())
        self.handleError(KekulizationFailure, self.validateKekulization())
        self.mol.openbonds = dict(self.openbonds)
        return self.mol

    def handleError(self, errtype, msg):
        if msg:
            if type(msg) == type(()):
                raise errtype(msg[0], self.smi, msg[1])
            else:
                raise errtype(msg, self.smi, self.idx)

    def handleComponent(self):
        # Check some conditions at the end of parsing a component
        if not self.rulesToIgnore & 1:
            if self.idx == 0 or self.smi[self.idx-1] in "." or (self.smi[self.idx-1] == '>' and self.reaction_part == 2):
                self.handleError(SMILESSyntaxError, "Empty molecules are not allowed")
        if self.bondchar:
            self.handleError(SMILESSyntaxError, "An atom must follow a bond symbol")
        if not self.rulesToIgnore & 4 and self.idx > 1 and self.smi[self.idx-1]==')':
            self.handleError(SMILESSyntaxError, ("The final branch should not be within parentheses", self.idx-1))

    def validateValence(self):
        # ----- Check for unusual valence -------
        for atom in self.mol.atoms:
            if not self.hasCommonValence(atom):
                return ("Uncommon valence or charge state", self.smiidx[atom.idx])

    def validateKekulization(self):
        # Try to kekulize all aromatic ring systems that are 'complete'
        # i.e. where no member has an unclosed ring opening or
        #      is the current atom being added to by the parser

        # To do this:
        # 1. Identify incomplete atoms
        # 2. Identify any aromatic ring system that contains an incomplete
        #    atom.
        # 3. Run the normal kekulization procedure after stripping aromaticity
        #    from the atoms in an incomplete system.

        incompleteAromaticAtoms = set(x for x in self.incompleteAtoms if x.arom)

        # Find aromatic systems (networks of aromatic bonds)
        seen = [0]*len(self.mol.atoms)
        bond_marker = [0]*len(self.mol.bonds)
        arom_system = 0
        incomplete_systems = set()
        for atom in self.mol.atoms:
            if seen[atom.idx] or not any(bond.arom for bond in atom.bonds): continue
            arom_system += 1
            stack = [atom]
            while len(stack):
                curr = stack.pop()
                if seen[curr.idx]: continue
                seen[curr.idx] = arom_system
                if curr in incompleteAromaticAtoms:
                    incomplete_systems.add(arom_system)
                for bond in curr.bonds:
                    if bond.arom:
                        bond_marker[bond.idx] = arom_system
                        nbr = bond.getNbr(curr)
                        if not seen[nbr.idx]:
                            stack.append(nbr)

        # Strip aromaticity
        for idx, bond in enumerate(self.mol.bonds):
            if bond_marker[idx] in incomplete_systems:
                bond.arom = False

        # Kekulize
        result = kekulize.Kekulize(self.mol)
        if result:
            return ("Aromatic system cannot be kekulized", self.smiidx[result])

        # Reset aromaticity
        if True: # set To False for a minor speed-up if you are not reusing
                  # the molecule afterwards.
            for idx, bond in enumerate(self.mol.bonds):
                if bond_marker[idx] in incomplete_systems:
                    bond.arom = True

        return None

    def validateSyntax(self, dot=False):
        # ----- Check syntax -------
        # Check that all ring bonds have been closed
        # Check that all brackets have been closed
        if self.partial and not dot:
            return None
        if (not dot or not self.rulesToIgnore & 16) and self.openbonds:
            text = "s have" if len(self.openbonds) > 1 else " has"
            return "{} ring opening{} not been closed".format(len(self.openbonds), text)
        if not self.rulesToIgnore & 8 and len(self.prev) > 1:
            text = "branches have" if len(self.prev) > 2 else "branch has"
            return "{} {} not been closed".format(len(self.prev)-1, text)
        return None

    def setImplicitHydrogenCount(self):
        for i, atom in enumerate(self.mol.atoms):
            hcount = self.hcount[i]
            if hcount > -1:
                atom.implh = hcount
            else:
                explicitvalence = self.getAdjustedExplicitValence(atom)
                implicitvalence = SmilesValence(atom.element, explicitvalence)
                implh = implicitvalence - explicitvalence
                if implh > 0 and (atom.arom or any(bond.arom for bond in atom.bonds)):
                    implh -= 1
                atom.implh = implh

    def handleBCSymbol(self):
        x = self.smi[self.idx]
        if x == "%":
            bcsymbol = self.smi[self.idx:self.idx+3]
            self.idx += 2
        else:
            bcsymbol = x
        if bcsymbol in self.openbonds:
            opening, openbc = self.openbonds.pop(bcsymbol)
            if opening == self.prev[-1]:
                return "Cannot have a bond opening and closing on the same atom"
            if self.mol.getBond(opening, self.prev[-1]):
                return "Cannot have a second bond between the same atoms"
            openbo = 0 if not openbc else ToBondOrder(openbc)
            closebo = 0 if not self.bondchar else ToBondOrder(self.bondchar)
            if closebo and openbo and closebo != openbo:
                return "Inconsistent bond orders"
            bo = closebo if closebo else openbo
            if self.bondchar == ':':
                arom = True
            else:
                arom = False
                if not bo:
                    if opening.arom and self.prev[-1].arom:
                        arom = True
                    bo = 1
            bond = self.mol.addBond(opening, self.prev[-1], bo)
            bond.arom = arom
        else:
            self.openbonds[bcsymbol] = (self.prev[-1], self.bondchar)
        self.idx += 1
        self.bondchar = None
        return None

    def getAdjustedExplicitValence(self, atom):
        """Adjust the explicit valence for the attachee in partial SMILES
        
        For example: "C(" has explicit valence of 0, but adjusted to 2
        For example: "C=" has explicit valence of 0, but adjusted to 2
        "C" -> 0
        "C(" -> 2
        
        """
        ans = atom.getExplicitValence()
        if not self.partial:
            return ans
        if atom == self.prev[-1]:
            if self.bondchar:
                ans += ToBondOrder(self.bondchar)
            elif self.parsing_atom or self.smi[-1] == '(' or (self.openbonds and len(self.prev) == 1 and self.smi[-1] != ')'):
                ans += 1
            if not self.rulesToIgnore & 4:
                if len(self.prev) > 1 and self.prev[-1] == self.prev[-2]:
                    ans += 1 # just opened a new branch
                if self.smi[-1] == ')':
                    ans += 1 # just closed a branch, but there must be another branch
        elif atom in self.prev:
            if not self.rulesToIgnore & 4:
                ans += 1
        for bcsymbol, (beg, symbol) in self.openbonds.items():
            if beg == atom:
                ans += 1 if not symbol else ToBondOrder(symbol)
        return ans

    def hasCommonValence(self, atom):
        data = valence.common_valencies.get(atom.element, None)

        # How to handle elements not in the list?
        if data is None:
            if atom.element == 0:
                return True # How to handle asterisk
            return False # Alternatively, you may wish to return True

        allowed = data.get(atom.charge, None)
        if allowed is None:
            return False # unusual charge state

        explval = self.getAdjustedExplicitValence(atom)
        # adjust valence for aromatic atoms (erring on the side of caution for incompleteAtoms)
        if valence.NeedsDblBond(atom) and not (self.partial and atom in self.incompleteAtoms):
            explval += 1

        totalvalence = explval + atom.implh
        if totalvalence in allowed:
            return True

        if self.partial and atom in self.incompleteAtoms:
            if totalvalence <= max(allowed):
                return True # still possibly normal valence

        # Note to reader: you could comment out the following line if
        # you don't need to support TEMPO-like (stable) oxygen radicals
        if atom.element==8 and atom.charge==0 and explval==1 and atom.implh==0 and valence.IsAttachedToNitrogen(atom): # TEMPO-like
            return True
        return False

    def notAtEnd(self):
        return self.idx < self.N

    def atEnd(self):
        return self.idx == self.N

    def incrementAndTestForEnd(self):
        self.idx += 1
        if self.atEnd():
            if self.partial:
                return True, None
            else:
                return True, "An open square brackets is present without the corresponding close square brackets"
        return False, None

    def parseAtom(self):
        x = self.smi[self.idx]
        self.parsing_atom = True

        if x == '[':
            end, msg = self.incrementAndTestForEnd()
            if end:
                return msg

            # Handle isotope
            if self.smi[self.idx].isdigit():
                isotope = int(self.smi[self.idx])
                if isotope == 0:
                    return "Isotope value of 0 not allowed"
                end, msg = self.incrementAndTestForEnd()
                if end:
                    return msg
                if self.smi[self.idx].isdigit():
                    isotope = isotope*10 + int(self.smi[self.idx])
                    end, msg = self.incrementAndTestForEnd()
                    if end:
                        return msg
                    if self.smi[self.idx].isdigit():
                        isotope = isotope*10 + int(self.smi[self.idx])
                        end, msg = self.incrementAndTestForEnd()
                        if end:
                            return msg
            else:
                isotope = 0

            # Handle element
            if self.smi[self.idx] not in firstLetterOfElements:
                return "An element symbol is required"
            if self.idx+1 < self.N and self.smi[self.idx].upper() + self.smi[self.idx+1] in elements:
                symbol = self.smi[self.idx:self.idx+2]
                self.idx += 1
            else:
                symbol = self.smi[self.idx]
            end, msg = self.incrementAndTestForEnd()
            if end:
                return msg

            # Handle tet stereo
            if self.smi[self.idx] == '@':
                end, msg = self.incrementAndTestForEnd()
                if end:
                    return msg
                if self.smi[self.idx] == '@':
                    end, msg = self.incrementAndTestForEnd()
                    if end:
                        return msg

            # Handle H count
            if self.smi[self.idx] == 'H':
                end, msg = self.incrementAndTestForEnd()
                if end:
                    return msg
                if self.smi[self.idx].isdigit():
                    hcount = int(self.smi[self.idx])
                    end, msg = self.incrementAndTestForEnd()
                    if end:
                        return msg
                else:
                    hcount = 1
            else:
                hcount = 0

            # Handle charge
            if self.smi[self.idx] in "+-":
                charge = 1 if self.smi[self.idx] == '+' else -1
                end, msg = self.incrementAndTestForEnd()
                if end:
                    return msg
                if self.smi[self.idx].isdigit():
                    numcharge = int(self.smi[self.idx])
                    charge *= numcharge
                    end, msg = self.incrementAndTestForEnd()
                    if end:
                        return msg
                elif self.smi[self.idx] in "+-":
                    numcharge = 1
                    while self.smi[self.idx] == self.smi[self.idx-1]:
                        numcharge += 1
                        end, msg = self.incrementAndTestForEnd()
                        if end:
                            return msg
                    charge *= numcharge
            else:
                charge = 0

            # Handle close-bracket
            if self.smi[self.idx] == ']':
                self.idx += 1
            else:
                return "Missing the close bracket"

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

        atom = self.mol.addAtom(symbol, self.prev[-1], self.bondchar if self.bondchar else "")
        atom.charge = charge
        atom.isotope = isotope
        self.prev[-1] = atom
        self.hcount.append(hcount)
        self.smiidx.append(self.idx-1)

        # Reset
        self.bondchar = None
        self.parsing_atom = False

valencemodel = {
        5: [3], 6: [4], 7: [3, 5], 15: [3, 5],
        8: [2], 16: [2, 4, 6],
        9: [1], 17: [1], 35: [1], 53: [1]
        }

def SmilesValence(elem, val):
    if elem == 0:
        return val
    implvalences = valencemodel[elem]
    for implvalence in implvalences:
        if val <= implvalence:
            return implvalence
    return val

def ToElement(symbol):
    if symbol[0].islower():
        symbol = symbol[0].upper() + symbol[1:]
    return elements[symbol]

###################### PUBLIC API ########################

def ParseSmiles(smi, partial=False, rulesToIgnore=0):
    sp = SmilesParser(partial, rulesToIgnore);
    mol = sp.parse(smi)
    return mol
