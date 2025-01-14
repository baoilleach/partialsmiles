import sys
import copy
from .elements import elements
from . import valence
from . import kekulize
from .exceptions import *

atomchars = set("CcONon[BPSFIbps*")
bondchars_list = "-=#$\\/:"
bondchars = set(bondchars_list)
bondorders = [1, 2, 3, 4, 1, 1, 1]
CLOSE, OPEN = range(2)

aromatic_list = ["c", "n", "p", "o", "s", "te", "se", "b"] # more?
aromatic_set = set(aromatic_list)
firstLetterOfElements = set(x[0] for x in list(elements.keys()) + aromatic_list)

def ToBondOrder(bondchar):
    idx = bondchars_list.index(bondchar)
    return bondorders[idx] # returns 1 for empty string

class Bond:
    __slots__ = ('beg', 'end', 'order', 'arom', 'idx')
    def __init__(self, beg, end, order):
        self.beg = beg
        self.end = end
        self.order = order
        self.arom = False
        self.idx = -1

    def getNbr(self, atom):
        return self.beg if self.end==atom else self.end

    def __repr__(self):
        return "Bond(beg={},end={},idx={},order={},arom={}".format(self.beg, self.end, self.idx, self.order, self.arom)

class Atom:
    __slots__ = ('element', 'charge', 'implh', 'arom', 'idx', 'bonds', 'isotope')
    def __init__(self, element, charge=0):
        self.element = element
        self.charge = charge
        self.implh = -1
        self.arom = False
        self.idx = -1
        self.bonds = []
        self.isotope = 0

    def __repr__(self):
        return "Atom(elem={},chg={},idx={},implh={},arom={})".format(self.element, self.charge, self.idx, self.implh, self.arom)

    def getExplicitDegree(self):
        return len(self.bonds)

    def getExplicitValence(self):
        return sum(bond.order for bond in self.bonds)

class Molecule:
    __slots__ = ('atoms', 'bonds', 'openbonds')
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

    def __repr__(self):
        return "Molecule(atoms={},bonds={})".format([str(x) for x in self.atoms], [str(x) for x in self.bonds])

class State:
    """The state of it"""
    __slots__ = ('mol', 'openbonds', 'hcount', 'idx', 'prev', 'smiidx', 'reaction_part', 'bondchar', 'parsing_atom')
    def __init__(self):
        self.mol = Molecule()
        self.openbonds = {} # dict from symbol -> (atom idx, bondchar)
        self.hcount = []
        self.idx = 0
        self.prev = [None]
        self.smiidx = []
        self.reaction_part = 0
        self.bondchar = None
        self.parsing_atom = False

    def __repr__(self):
        """Look at the state of it"""
        out = f"{self.mol}\n"
        out += f"{self.openbonds=}\n{self.hcount=}\n{self.idx=}\n{self.prev=}\n{self.smiidx=}\n{self.reaction_part=}\n{self.bondchar=}\n{self.parsing_atom=}"
        return out

class SmilesParser:

    def __init__(self, partial, rulesToIgnore):
        self.partial = partial
        self.rulesToIgnore = rulesToIgnore

    def parse(self, smi):
        self.smi = smi
        self.N = len(smi)
        state = State()

        while state.idx < self.N:
            self.parse_token(state)

        if not self.partial:
            self.handleComponent(state)

        self.handleError(SMILESSyntaxError, (self.validateSyntax(state), state.idx))

        self.setImplicitHydrogenCount(state)
        self.incompleteAtoms = set(x[0] for x in state.openbonds.values())
        if self.partial:
            self.incompleteAtoms.update(x for x in state.prev if x is not None)

        self.handleError(ValenceError, self.validateValence(state))
        self.handleError(KekulizationFailure, self.validateKekulization(state.mol, state))
        
        return state

    def parse_token(self, state):
        x = self.smi[state.idx]
        if x in atomchars:
            self.handleError(SMILESSyntaxError, (self.parseAtom(state), state.idx))
        elif x == ')':
            if not self.rulesToIgnore & 2 and (state.idx > 1 and self.smi[state.idx-1]=='('):
                self.handleError(SMILESSyntaxError, ("Empty branches are not allowed", state.idx))
            if not self.rulesToIgnore & 4 and state.idx > 1 and self.smi[state.idx-1]==')':
                self.handleError(SMILESSyntaxError, ("The final branch should not be within parentheses", state.idx-1))
            if state.bondchar:
                self.handleError(SMILESSyntaxError, ("An atom must follow a bond symbol", state.idx))
            state.prev.pop()
            if not state.prev:
                self.handleError(SMILESSyntaxError, ("Unmatched close parenthesis", state.idx))
            state.idx += 1
        elif x == '(':
            if state.prev[-1] is None or self.smi[state.idx-1]=='(':
                self.handleError(SMILESSyntaxError, ("An atom must precede an open parenthesis", state.idx))
            if state.bondchar:
                self.handleError(SMILESSyntaxError, ("A bond symbol should not precede an open parenthesis", state.idx))
            state.prev.append(state.prev[-1])
            state.idx += 1
        elif x in bondchars:
            if state.prev[-1] is None:
                self.handleError(SMILESSyntaxError, ("An atom must precede a bond symbol", state.idx))
            if state.bondchar:
                self.handleError(SMILESSyntaxError, ("Only a single bond symbol should be used", state.idx))
            if not self.rulesToIgnore & 64 and x == ':':
                self.handleError(SMILESSyntaxError, ("Aromatic bond symbols are rejected by default", state.idx))
            state.bondchar = x
            state.idx += 1
        elif x.isdigit() or x=='%':
            if state.prev[-1] is None:
                self.handleError(SMILESSyntaxError, ("An atom must precede a bond closure symbol", state.idx))
            if not self.rulesToIgnore & 32:
                precedingtext = self.smi[state.smiidx[state.prev[-1].idx]+1:state.idx]
                if ")" in precedingtext:
                    self.handleError(SMILESSyntaxError, ("Ring closure symbols must immediately follow an atom", state.idx))
                if "(" in precedingtext:
                    self.handleError(SMILESSyntaxError, ("Ring closure symbols should not be in parentheses", state.idx))
            self.handleError(SMILESSyntaxError, (self.handleBCSymbol(state), state.idx))
        elif x in '.>':
            self.handleComponent(state)

            state.prev[-1] = None
            if x == '.':
                self.handleError(SMILESSyntaxError, (self.validateSyntax(state, dot=True), state.idx))
                state.idx += 1
            elif x == '>':
                state.reaction_part += 1
                if state.reaction_part == 3:
                    self.handleError(SMILESSyntaxError, ("Reactions only have three parts", state.idx))
                self.handleError(SMILESSyntaxError, (self.validateSyntax(state, dot=True), state.idx))
                state.idx += 1
        else:
            self.handleError(SMILESSyntaxError, ("Illegal character", state.idx))

    def handleError(self, errtype, errinfo):
        """Do nothing if no errorinfo supplied, otherwise raise an Exception"""
        if errinfo and errinfo[0]:
            raise errtype(errinfo[0], self.smi, errinfo[1])

    def handleComponent(self, state):
        # Check some conditions at the end of parsing a component
        if not self.rulesToIgnore & 1:
            if state.idx == 0 or self.smi[state.idx-1] in "." or (self.smi[state.idx-1] == '>' and state.reaction_part == 2):
                self.handleError(SMILESSyntaxError, ("Empty molecules are not allowed", state.idx))
        if state.bondchar:
            self.handleError(SMILESSyntaxError, ("An atom must follow a bond symbol", state.idx))
        if not self.rulesToIgnore & 4 and state.idx > 1 and self.smi[state.idx-1]==')':
            self.handleError(SMILESSyntaxError, ("The final branch should not be within parentheses", state.idx-1))

    def validateValence(self, state):
        # ----- Check for unusual valence -------
        for atom in state.mol.atoms:
            if not self.hasCommonValence(atom, state):
                return ("Uncommon valence or charge state", state.smiidx[atom.idx])

    def validateKekulization(self, mol, state):
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
        seen = [0]*len(mol.atoms)
        bond_marker = [0]*len(mol.bonds)
        arom_system = 0
        incomplete_systems = set()
        for atom in mol.atoms:
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
        for idx, bond in enumerate(mol.bonds):
            if bond_marker[idx] in incomplete_systems:
                bond.arom = False

        # Kekulize
        result = kekulize.Kekulize(mol)
        if result is not None:
            return ("Aromatic system cannot be kekulized", state.smiidx[result])

        # Reset aromaticity
        if True: # set To False for a minor speed-up if you are not reusing
                  # the molecule afterwards.
            for idx, bond in enumerate(mol.bonds):
                if bond_marker[idx] in incomplete_systems:
                    bond.arom = True

        return None

    def validateSyntax(self, state, dot=False):
        # ----- Check syntax -------
        # Check that all ring bonds have been closed
        # Check that all brackets have been closed
        if self.partial and not dot:
            return None
        if (not dot or not self.rulesToIgnore & 16) and state.openbonds:
            text = "s have" if len(state.openbonds) > 1 else " has"
            return "{} ring opening{} not been closed".format(len(state.openbonds), text)
        if not self.rulesToIgnore & 8 and len(state.prev) > 1:
            text = "branches have" if len(state.prev) > 2 else "branch has"
            return "{} {} not been closed".format(len(state.prev)-1, text)
        return None

    def setImplicitHydrogenCount(self, state):
        for i, atom in enumerate(state.mol.atoms):
            hcount = state.hcount[i]
            if hcount > -1:
                atom.implh = hcount
            else:
                explicitvalence = self.getAdjustedExplicitValence(atom, state)
                implicitvalence = SmilesValence(atom.element, explicitvalence)
                implh = implicitvalence - explicitvalence
                if implh > 0 and (atom.arom or any(bond.arom for bond in atom.bonds)):
                    implh -= 1
                atom.implh = implh

    def handleBCSymbol(self, state):
        x = self.smi[state.idx]
        if x == "%":
            bcsymbol = self.smi[state.idx:state.idx+3]
            state.idx += 2
        else:
            bcsymbol = x
        if bcsymbol in state.openbonds:
            opening, openbc = state.openbonds.pop(bcsymbol)
            if opening == state.prev[-1]:
                return "Cannot have a bond opening and closing on the same atom"
            if state.mol.getBond(opening, state.prev[-1]):
                return "Cannot have a second bond between the same atoms"
            openbo = 0 if not openbc else ToBondOrder(openbc)
            closebo = 0 if not state.bondchar else ToBondOrder(state.bondchar)
            if closebo and openbo and closebo != openbo:
                return "Inconsistent bond orders"
            bo = closebo if closebo else openbo
            if state.bondchar == ':':
                arom = True
            else:
                arom = False
                if not bo:
                    if opening.arom and state.prev[-1].arom:
                        arom = True
                    bo = 1
            bond = state.mol.addBond(opening, state.prev[-1], bo)
            bond.arom = arom
        else:
            state.openbonds[bcsymbol] = (state.prev[-1], state.bondchar)
        state.idx += 1
        state.bondchar = None
        return None

    def getAdjustedExplicitValence(self, atom, state):
        """Adjust the explicit valence for the attachee in partial SMILES
        
        For example: "C(" has explicit valence of 0, but adjusted to 2
        For example: "C=" has explicit valence of 0, but adjusted to 2
        "C" -> 0
        "C(" -> 2
        
        """
        ans = atom.getExplicitValence()
        if not self.partial:
            return ans
        if atom == state.prev[-1]:
            if state.bondchar:
                ans += ToBondOrder(state.bondchar)
            elif state.parsing_atom or self.smi[-1] == '(' or (state.openbonds and len(state.prev) == 1 and self.smi[-1] != ')'):
                ans += 1
            if not self.rulesToIgnore & 4:
                if len(state.prev) > 1 and state.prev[-1] == state.prev[-2]:
                    ans += 1 # just opened a new branch
                if self.smi[-1] == ')':
                    ans += 1 # just closed a branch, but there must be another branch
        elif atom in state.prev:
            if not self.rulesToIgnore & 4:
                ans += 1
        for bcsymbol, (beg, symbol) in state.openbonds.items():
            if beg == atom:
                ans += 1 if not symbol else ToBondOrder(symbol)
        return ans

    def hasCommonValence(self, atom, state):
        data = valence.common_valencies.get(atom.element, None)

        # How to handle elements not in the list?
        if data is None:
            if atom.element == 0:
                return True # How to handle asterisk
            return False # Alternatively, you may wish to return True

        allowed = data.get(atom.charge, None)
        if allowed is None:
            return False # unusual charge state

        explval = self.getAdjustedExplicitValence(atom, state)
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

    def notAtEnd(self, state):
        return state.idx < self.N

    def atEnd(self, state):
        return state.idx == self.N

    def incrementAndTestForEnd(self, state):
        state.idx += 1
        if self.atEnd(state):
            if self.partial:
                return True, None
            else:
                return True, "An open square brackets is present without the corresponding close square brackets"
        return False, None

    def parseAtom(self, state):
        x = self.smi[state.idx]
        state.parsing_atom = True

        if x == '[':
            end, msg = self.incrementAndTestForEnd(state)
            if end:
                return msg

            # Handle isotope
            if self.smi[state.idx].isdigit():
                isotope = int(self.smi[state.idx])
                if isotope == 0:
                    return "Isotope value of 0 not allowed"
                end, msg = self.incrementAndTestForEnd(state)
                if end:
                    return msg
                if self.smi[state.idx].isdigit():
                    isotope = isotope*10 + int(self.smi[state.idx])
                    end, msg = self.incrementAndTestForEnd(state)
                    if end:
                        return msg
                    if self.smi[state.idx].isdigit():
                        isotope = isotope*10 + int(self.smi[state.idx])
                        end, msg = self.incrementAndTestForEnd(state)
                        if end:
                            return msg
            else:
                isotope = 0

            # Handle element
            if self.smi[state.idx] not in firstLetterOfElements:
                return "An element symbol is required"
            if state.idx+1 < self.N and self.smi[state.idx].upper() + self.smi[state.idx+1] in elements:
                symbol = self.smi[state.idx:state.idx+2]
                state.idx += 1
            else:
                symbol = self.smi[state.idx]
                if symbol.upper() not in elements:
                    return "An element symbol is required"
            end, msg = self.incrementAndTestForEnd(state)
            if end:
                return msg

            # Handle tet stereo
            if self.smi[state.idx] == '@':
                end, msg = self.incrementAndTestForEnd(state)
                if end:
                    return msg
                if self.smi[state.idx] == '@':
                    end, msg = self.incrementAndTestForEnd(state)
                    if end:
                        return msg

            # Handle H count
            if self.smi[state.idx] == 'H':
                end, msg = self.incrementAndTestForEnd(state)
                if end:
                    return msg
                if self.smi[state.idx].isdigit():
                    hcount = int(self.smi[state.idx])
                    end, msg = self.incrementAndTestForEnd(state)
                    if end:
                        return msg
                else:
                    hcount = 1
            else:
                hcount = 0

            # Handle charge
            if self.smi[state.idx] in "+-":
                charge = 1 if self.smi[state.idx] == '+' else -1
                end, msg = self.incrementAndTestForEnd(state)
                if end:
                    return msg
                if self.smi[state.idx].isdigit():
                    numcharge = int(self.smi[state.idx])
                    charge *= numcharge
                    end, msg = self.incrementAndTestForEnd(state)
                    if end:
                        return msg
                elif self.smi[state.idx] in "+-":
                    numcharge = 1
                    while self.smi[state.idx] == self.smi[state.idx-1]:
                        numcharge += 1
                        end, msg = self.incrementAndTestForEnd(state)
                        if end:
                            return msg
                    charge *= numcharge
            else:
                charge = 0

            # Handle close-bracket
            if self.smi[state.idx] == ']':
                state.idx += 1
            else:
                return "Missing the close bracket"

        elif self.notAtEnd(state) and self.smi[state.idx:state.idx+2] in ["Cl", "Br"]:
            symbol = self.smi[state.idx:state.idx+2]
            state.idx += 2
            hcount = -1
            charge = 0
            isotope = 0
        else:
            state.idx += 1
            symbol = x
            hcount = -1
            charge = 0
            isotope = 0

        atom = state.mol.addAtom(symbol, state.prev[-1], state.bondchar if state.bondchar else "")
        atom.charge = charge
        atom.isotope = isotope
        state.prev[-1] = atom
        state.hcount.append(hcount)
        state.smiidx.append(state.idx-1)

        # Reset
        state.bondchar = None
        state.parsing_atom = False

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
    state = sp.parse(smi)
    return state
