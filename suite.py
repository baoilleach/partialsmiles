import unittest
from partialsmiles.smiparser import SmilesParser, IndexOfLastToken
from partialsmiles.exceptions import *

class ValenceTests(unittest.TestCase):

    def testAdjustedValence(self):
        data = [
                ("C", 0), # might turn out to be the final char
                ("C[", 1),
                ("C(", 2),
                ("C1(", 3),
                ("C=[", 2),
                ("C=", 2),
                ("C.", 0),
                ("C1", 2),
                ("C=1", 3),
                ("C(C)", 2), # rule 3
                ]
        for smi, val in data:
            sp = SmilesParser()
            state = sp.parse(smi, True)
            mol = state.mol
            atom = mol.atoms[0]
            self.assertEqual(sp.getAdjustedExplicitValence(atom, state), val)
        data = [
                ("C(", 2),
                ("C(C)(", 3),
                ("S(=C)(=", 5),
                ]
        for smi, val in data:
            sp = SmilesParser()
            state = sp.parse(smi, True)
            mol = state.mol
            atom = mol.atoms[0]
            self.assertEqual(sp.getAdjustedExplicitValence(atom, state), val)

    def testBranched(self):
        smi = "C(C(C)("
        val = [2, 4, 1]
        sp = SmilesParser()
        state = sp.parse(smi, True)
        mol = state.mol
        for v, atom in zip(val, mol.atoms):
            self.assertEqual(sp.getAdjustedExplicitValence(atom, state), v)

        smi = "C(C(C)"
        val = [2, 3, 1]
        sp = SmilesParser()
        state = sp.parse(smi, True)
        mol = state.mol
        for v, atom in zip(val, mol.atoms):
            self.assertEqual(sp.getAdjustedExplicitValence(atom, state), v)

    def testParser(self):
        smis = ["C(C)(C)(C)(C)C", "[N+5]"]
        sp = SmilesParser()
        for smi in smis:
            self.assertRaises(ValenceError, sp.parse, smi, False)
            self.assertRaises(ValenceError, sp.parse, smi, True)

    def testInvalid(self):
        sp = SmilesParser()
        state = sp.parse("[CH3]C", True)
        self.assertRaises(ValenceError, sp.parse, "[CH2]C", True)
        self.assertRaises(ValenceError, sp.parse, "[CH3]1C", True)
        self.assertRaises(ValenceError, sp.parse, "[CH2]=1C", True)

    def testCoverage(self):
        sp = SmilesParser()
        for smi in ["[Eu]C", "[N](=C)(C)(C)C"]:
            self.assertRaises(ValenceError, sp.parse, smi, True)
        # TEMPO
        state = sp.parse("C1(C)(C)N([O])C(C)(C)CCC1", partial=False)
        state = sp.parse("FC[CH2]", partial=True)

    def testBugs(self):
        sp = SmilesParser()
        # Don't raise valence errors for the following
        good = ["c1O", "c1cc(sc1="]
        for smi in good:
            state = sp.parse(smi, partial=True)

        self.assertRaises(ValenceError, sp.parse, "CN(C)=", True)
        # The following is 5-valent at least, and we only allow 3-valent
        # nitrogens.
        self.assertRaises(ValenceError, sp.parse, "CCN(=O)(", True)
        # The second bond closure was not contributing to the
        # adjusted valence due to a 'break' statement
        self.assertRaises(ValenceError, sp.parse, "C=C3=1", True)

    def testImpliedValence(self):
        sp = SmilesParser()
        sp_lenient = SmilesParser(4)
        # Given Rule 3, an open parenthesis implies at least two
        # nbrs.
        self.assertRaises(ValenceError, sp.parse, "C1CN1(", True)
        state = sp_lenient.parse("C1CN1(", True)
        self.assertRaises(ValenceError, sp.parse, "CN(=", True)
        state = sp_lenient.parse("CN(=", True)

    def testBug(self):
        sp = SmilesParser()
        # The latter was failing with a valence error, but the former was
        # not
        smis = ["ncC(N)(=O)[O", "ncC(N)(=O)["]
        for smi in smis:
            self.assertRaises(ValenceError, sp.parse, smi, True)

    def testEffectOfOpenRing(self):
        # If there's a ring opening, there must be an additional
        # valence (if at the outermost branching level).
        sp = SmilesParser()
        sp.parse("C1C=N", True)
        self.assertRaises(ValenceError, sp.parse, "C1C=N2", True)
        sp.parse("C1C=NC", True)
        sp.parse("c1(=O", True)
        sp.parse("c1(=O)", True)

class MolTest(unittest.TestCase):

    def testBonds(self):
        sp = SmilesParser()
        state = sp.parse("CCO", False)
        mol = state.mol
        self.assertTrue(mol.getBond(mol.atoms[0], mol.atoms[1]))
        self.assertTrue(mol.getBond(mol.atoms[1], mol.atoms[2]))
        self.assertFalse(mol.getBond(mol.atoms[0], mol.atoms[2]))

    def testCharge(self):
        data = [
                ("[NH4+]", 1),
                ("[CH3-]", -1),
                ("[CH2--]", -2),
               ]
        sp = SmilesParser()
        for smi, charge in data:
            state = sp.parse(smi, False)
            mol = state.mol
            self.assertTrue(mol.atoms[0].charge, charge)

    def testAsterisk(self):
        sp = SmilesParser()
        state = sp.parse("*", True)
        mol = state.mol
        self.assertEqual(mol.atoms[0].element, 0)
        self.assertEqual(mol.atoms[0].implh, 0)
        state = sp.parse("[*H]", True)
        mol = state.mol
        self.assertEqual(mol.atoms[0].element, 0)
        self.assertEqual(mol.atoms[0].implh, 1)

class RuleTests(unittest.TestCase):

    def testRuleOne(self):
        """Whether to allow empty molecules, like C..C"""
        sp = SmilesParser()
        sp_lenient = SmilesParser(1)
        smi = "C..C"
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
        sp_lenient.parse(smi, False)

    def testRuleTwo(self):
        "Whether to allow an empty branch like C()C"
        sp = SmilesParser()
        sp_lenient = SmilesParser(2)
        smi = "C()C"
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
        sp_lenient.parse(smi, False)

    def testRuleThree(self):
        """Whether to require that the final branch be unparenthesised"""
        sp = SmilesParser()
        sp_lenient = SmilesParser(4)
        smis = ["C(O(C))", "C(O).C"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
            sp_lenient.parse(smi, False)

        smis = ["C(O)"] # these are fine as partial but not otherwise
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
            sp.parse(smi, True)
            sp_lenient.parse(smi, False)

    def testRuleFour(self):
        """Whether to allow a dot within parentheses"""
        sp = SmilesParser()
        sp_lenient = SmilesParser(8)
        smis = ["C(C.C)C"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
            sp_lenient.parse(smi, False)

    def testRuleFive(self):
        """Whether to allow bond closures across disconnected components"""
        sp = SmilesParser()
        sp_lenient = SmilesParser(16)
        smis = ["C1.C1"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
            sp_lenient.parse(smi, False)

    def testRuleSix(self):
        """Whether to require ring closure symbols immediately after atoms"""
        sp = SmilesParser()
        sp_lenient = SmilesParser(32)
        smis = ["C(1)CCC1", "C(CCC1)1"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
            sp_lenient.parse(smi, False)

    def testRuleSeven(self):
        """Whether to reject aromatic bond symbols"""
        # Aromatic bond symbols are not needed in a SMILES string
        # - an aromatic bond connecting two lowercase atoms is implicit
        # - the meaning of an aromatic bond attached to an uppercase atom
        #   is a matter of intense debate
        sp = SmilesParser()
        sp_lenient = SmilesParser(64)
        smi = "c1:c:c:c:c:c1"
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)
        self.assertRaises(SMILESSyntaxError, sp.parse, smi, True)
        sp_lenient.parse(smi, False)

class ParserTests(unittest.TestCase):

    def check(self, good, bad):
        sp = SmilesParser()
        for smi in good:
            sp.parse(smi, partial=False)
        for smi in bad:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, partial=False)

    def testIllegalChar(self):
        self.check(["C"], ["!", "CC[L]"])

    def testParentheses(self):
        good = ["C(=O)Cl"]
        bad = ["C(", "C(C(", "C)", "CCC.)", "(C)", "C((C))",
               "C.(C)", ")C", "C(C))C"]
        self.check(good, bad)

    def testBondClosures(self):
        good = ["C1CCC1", "C%23CCC%23", "C-1OC1", "C-1OC-1"]
        bad = ["C1CC", "C11C", "1C", "%12C", "C.1C", "C-1OC=1", "C1C1"]
        self.check(good, bad)

    def testDots(self):
        good = ["C.C"]
        bad = [".C", "C..C", "C-.", "C.", "C. "]
        self.check(good, bad)

    def testBondChar(self):
        good = ["C-C#C", "C/C=C/Cl"]
        bad = ["-C", "C.-C", "C-=C", "C--C", "C-(C)", "C=", "C= "]
        self.check(good, bad)

    def testIsotope(self):
        good = ["[12CH4]"]
        bad = ["[0CH4]", "[1", "[12", "[123", "[v", "[C@", "[C@@",
               "[CH", "[CH4", "[C+", "[C+2", "[C++", "[C++ "]
        self.check(good, bad)
        sp = SmilesParser()
        self.assertEqual(sp.parse("[12CH4]").mol.atoms[0].isotope, 12)
        self.assertEqual(sp.parse("[123CH4]").mol.atoms[0].isotope, 123)

    def testBrackets(self):
        good = ["[16CH4]", "[NH3++]", "C[C@@H](C)C"]
        bad = ["[16C)", "C(C-)Cl", "[C"]
        self.check(good, bad)

    def testCharge(self):
        good = ["[NH4+]", "[NH3++]", "[CH3-]"]
        bad = ["[N4+H]", "[NH+-]"]
        self.check(good, bad)

    def testTripTests(self):
        # Roger's from https://www.nextmovesoftware.com/talks/Sayle_SmilesPlus_InChI_201908.pdf
        smis = """C1.C1
# C%00CC%00 # we treat as 0
C(C.C)C
C(C)1CC1
C(.C)
C()
(CO)=O
(C)
.C
C..C
C.
C=(O)C
C((C))
C.(C)
C1CC(=1)
C1CC(1)
C(C.)
C==C""".split("\n")
        sp = SmilesParser()
        for smi in smis:
            if smi.startswith("#"): continue
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, False)

class KekulizationTests(unittest.TestCase):

    def testBasic(self):
        sp = SmilesParser()
        state = sp.parse("c1ccccc1", False) # Will be kekulized
        mol = state.mol
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 3)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

        state = sp.parse("c1ccccc1.", partial=True) # Will be kekulized
        mol = state.mol
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 3)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

        state = sp.parse("c1ccccc1", True) # Won't be kekulized
        mol = state.mol
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 6)
        self.assertEqual(sum(bond.arom for bond in mol.bonds), 6)

        state = sp.parse("c1ccccc1C", True) # Will be kekulized
        mol = state.mol
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 4)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

    def testHardCases(self):
        smis = ["c12c3c4c5c1c6c7c8c9c6c%10c5c%11c%12c4c%13c%14c%15c-%16c%17c%13c%12c%18c%19c%17c%20c%16c%21c%22=c%23c%24c(c8c%25C%24%26C%27(C%26)c%22c%20c%28c%19c%29c%30c(c%27%28)-c%25c9c%30c%10c%11c%29%18)c%31c7c2c%32-c%31c%23C%21%33C%15(c%32c%143)C%33",
                 "c12c3c4c5c6c1c-7c8c9c2c%10c3c%11c%12c%13c%14c%15c%16c%17c%13c%11c4c%18c%17c-%19c%20C%21%22c%23c%24-c%25c%26C%21(C%22)c(c%19%16)c-%27c%15C%28%29c%14c-%30c%12c%10c%31c9c%32c-%33c(c%31%30)C%28(C%29)c(c%26%27)c%33c%25C%34(C%328C%34)c%24c7c6c%23c%20c5%18",
                 "c12c3c4c5c6c7c4c8c9c%10c%11c8c%12c7c%13c-%14c6c%15c%16c%14c%17c%18c%19c%20c%21=c%22c%23c%24c(c1c(c93)c%25c%24c%26c%22c%27c%21c%18c%28C%17%29C%13(c%12c%30-c%28c%27c(c%11%30)c%26c%10%25)C%31C%29C(=N)CCC%31)c%32c-%33c2c5c%15c%33c(c%19%16)C%20C%32%23",
                 "c12c3c4c5c6c7c8c5c9c%10c4c1c%11c%12c%13c2c%14c%15c3c6c%16c%15c%17c%18c%19c%16c7c%20c%19c%21c%22c%18c%23c%17c%14c%13c%24c%23c%25c%22c%26c%27c%28c%29c(c%12c%24c%28%25)-c%11c%10C%30%31C%29(c%27c-%32c(c%26%21)c%20c8c%32c%309)CC(=N)CC%31"
                 ]
        sp = SmilesParser()
        for smi in smis:
            state = sp.parse(smi, False)

    def testFailures(self):
        smis = [
                "n1c[n]nn1",
                "c12ncncc1c(=C)[n]n2C",
                "c1snn2c1nc1ccccc21",
                "s1cc2nc3c(n2n1)cccc3"
                ]
        sp = SmilesParser()
        for smi in smis:
            self.assertRaises(KekulizationFailure, sp.parse, smi, False)

    
    def testCoverage(self):
        smis = [
                "Cs1(=O)ccccn1", # S-oxides
                "c1c#cccc1", # benzyne
                "n1cnc[n-]1", # from SMILES benchmark
                "[nH]1cnnn1", # ditto
                "n1cnco1",    # ditto
                "n1(cno[cH+]1)C", # ditto
                "O[p+]1(O)cccc1", # invented molecule
                "[s-]1(O)cccc1"   # invented molecule
                ]
        sp = SmilesParser()
        for smi in smis:
            state = sp.parse(smi, partial=False)

class AromaticBondTest(unittest.TestCase):
    
    def testBasic(self):
        smis = ["cc", "c:c", "C:C"] # What about c=c? Is this [C]=[C]?
        sp = SmilesParser(rulesToIgnore=64)
        for smi in smis:
            state = sp.parse(smi)
            mol = state.mol
            self.assertEqual(2, mol.bonds[0].order)
            self.assertTrue(mol.bonds[0].arom)

    def testRingClosures(self):
        smis = ["c1.c1", "c1.c:1", "C1.C:1"]
        sp = SmilesParser(rulesToIgnore=16|64)
        for smi in smis:
            state = sp.parse(smi)
            mol = state.mol
            self.assertEqual(2, mol.bonds[0].order)
            self.assertTrue(mol.bonds[0].arom)

class ReactionTests(unittest.TestCase):

    def testBasic(self):
        
        good = ["C>O>N", "C>>N"]
        empty = ["C>>", ">>N", ">O>", "C.>>.N"]
        bad = ["C>O>N>S"]
        sp = SmilesParser()
        for smi in good:
            state = sp.parse(smi, partial=False)
        for smi in bad + empty:
            self.assertRaises(SMILESSyntaxError, sp.parse, smi, partial=False)
        sp = SmilesParser(rulesToIgnore=1)
        for smi in empty:
            state = sp.parse(smi, partial=False)

class LastToken(unittest.TestCase):

    def testBasic(self):
        data = [("CC", 1),
                ("C", 0),
                ("r", -1),
                ("", -1),
                ("CCl", 1),
                ("C[fdjkdfk]", 1),
                ("Cfdjkdfk]", -1),
                ("C1", 1),
                ("C12", 2),
                ("C%12", 1),
                ("C%", 1),
                ("C=", 1),
                ("C(", 1),
                ("1", 0),
                ("12", 1),
                ("%12", 0),
                ]
        for smi, idx in data:
            ans = IndexOfLastToken(smi)
            self.assertEqual(idx, ans)

class CacheTests(unittest.TestCase):

    def testBasic(self):
        parser = SmilesParser()

        state = parser.parse("CO", partial=True)
        self.assertFalse(state.from_cache)
        state = parser.parse("CN", partial=True)
        self.assertTrue(state.from_cache)
        state = parser.parse("ON", partial=True)
        self.assertFalse(state.from_cache)

if __name__ == "__main__":
    unittest.main()
