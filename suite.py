import unittest
from partialsmiles.smiparser import SmilesParser, ParseSmiles
from partialsmiles.exceptions import *

class ValenceTests(unittest.TestCase):

    def testAdjustedValence(self):
        data = [
                ("C", 0), # might turn out to be the final char
                ("C[", 1),
                ("C1(", 2),
                ("C=[", 2),
                ("C=", 2),
                ("C.", 0),
                ("C1", 1),
                ("C=1", 2),
                ("C(C)", 1),
                ]
        for smi, val in data:
            sp = SmilesParser(True, 0)
            mol = sp.parse(smi)
            atom = mol.atoms[-1]
            self.assertEqual(sp.getAdjustedExplicitValence(atom), val)
        data = [
                ("C(", 1),
                ("C(C)(", 2),
                ("C(=C)(=", 4),
                ]
        for smi, val in data:
            sp = SmilesParser(True, 0)
            mol = sp.parse(smi)
            atom = mol.atoms[0]
            self.assertEqual(sp.getAdjustedExplicitValence(atom), val)

    def testParser(self):
        smis = ["C(C)(C)(C)(C)C", "[N+5]"]
        for smi in smis:
            self.assertRaises(ValenceError, ParseSmiles, smi, False)
            self.assertRaises(ValenceError, ParseSmiles, smi, True)

    def testInvalid(self):
        mol = ParseSmiles("[CH3]C", True)
        self.assertRaises(ValenceError, ParseSmiles, "[CH2]C", True)
        self.assertRaises(ValenceError, ParseSmiles, "[CH3]1C", True)
        self.assertRaises(ValenceError, ParseSmiles, "[CH2]=1C", True)

    def testCoverage(self):
        for smi in ["[Eu]C", "[N](=C)(C)(C)C"]:
            self.assertRaises(ValenceError, ParseSmiles, smi, True)
        # TEMPO
        mol = ParseSmiles("C1(C)(C)N([O])C(C)(C)CCC1", partial=False)
        mol = ParseSmiles("FC[CH2]", partial=True)

    def testBugs(self):
        # Don't raise valence errors for the following
        mol = ParseSmiles("c1O", partial=True)
        mol = ParseSmiles("c1cc(sc1=", partial=True)

class MolTest(unittest.TestCase):

    def testBonds(self):
        mol = ParseSmiles("CCO", False)
        self.assertTrue(mol.getBond(mol.atoms[0], mol.atoms[1]))
        self.assertTrue(mol.getBond(mol.atoms[1], mol.atoms[2]))
        self.assertFalse(mol.getBond(mol.atoms[0], mol.atoms[2]))

    def testCharge(self):
        data = [
                ("[NH4+]", 1),
                ("[CH3-]", -1),
                ("[CH2--]", -2),
               ]
        for smi, charge in data:
            mol = ParseSmiles(smi, False)
            self.assertTrue(mol.atoms[0].charge, charge)

    def testAsterisk(self):
        mol = ParseSmiles("*", True)
        self.assertEqual(mol.atoms[0].element, 0)
        self.assertEqual(mol.atoms[0].implh, 0)
        mol = ParseSmiles("[*H]", True)
        self.assertEqual(mol.atoms[0].element, 0)
        self.assertEqual(mol.atoms[0].implh, 1)

class RuleTests(unittest.TestCase):

    def testRuleOne(self):
        """Whether to allow empty molecules, like C..C"""
        smi = "C..C"
        self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, False)
        self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, True)
        ParseSmiles(smi, False, 1)

    def testRuleTwo(self):
        "Whether to allow an empty branch like C()C"
        smi = "C()C"
        self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, False)
        self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, True)
        ParseSmiles(smi, False, 2)

    def testRuleThree(self):
        """Whether to allow an open parenthesis without a preceding
        atom like (CC)"""
        smis = ["(CC)", "C.(CC)"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, False)
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, True)
            ParseSmiles(smi, False, 4)

    def testRuleFour(self):
        """Whether to allow a dot within parentheses"""
        smis = ["C(C.C)C"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, False)
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, True)
            ParseSmiles(smi, False, 8)

    def testRuleFive(self):
        """Whether to allow bond closures across disconnected components"""
        smis = ["C1.C1"]
        for smi in smis:
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, False)
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, True)
            ParseSmiles(smi, False, 16)

class ParserTests(unittest.TestCase):

    def check(self, good, bad):
        for smi in good:
            ParseSmiles(smi, partial=False)
        for smi in bad:
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, partial=False)

    def testIllegalChar(self):
        self.check(["C"], ["!"])

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

    def testBrackets(self):
        good = ["[16CH4]", "[NH3++]", "C[C@@H](C)C"]
        bad = ["[16C)", "C(C-)Cl", "[C"]
        self.check(good, bad)

    def testCharge(self):
        good = ["[NH4+]", "[NH3++]", "[CH3-]"]
        bad = ["[N4+H]", "[NH+-]"]
        self.check(good, bad)

class KekulizationTests(unittest.TestCase):

    def testBasic(self):
        mol = ParseSmiles("c1ccccc1", False) # Will be kekulized
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 3)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

        mol = ParseSmiles("c1ccccc1 ", False) # Will be kekulized
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 3)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

        mol = ParseSmiles("c1ccccc1.", partial=True) # Will be kekulized
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 3)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

        mol = ParseSmiles("c1ccccc1", True) # Won't be kekulized
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 6)
        self.assertEqual(sum(bond.arom for bond in mol.bonds), 6)

        mol = ParseSmiles("c1ccccc1C", True) # Will be kekulized
        bondorders = [bond.order for bond in mol.bonds]
        self.assertEqual(sum(x==1 for x in bondorders), 4)
        self.assertEqual(sum(x==2 for x in bondorders), 3)

    def testHardCases(self):
        smis = ["c12c3c4c5c1c6c7c8c9c6c%10c5c%11c%12c4c%13c%14c%15c-%16c%17c%13c%12c%18c%19c%17c%20c%16c%21c%22=c%23c%24c(c8c%25C%24%26C%27(C%26)c%22c%20c%28c%19c%29c%30c(c%27%28)-c%25c9c%30c%10c%11c%29%18)c%31c7c2c%32-c%31c%23C%21%33C%15(c%32c%143)C%33",
                 "c12c3c4c5c6c1c-7c8c9c2c%10c3c%11c%12c%13c%14c%15c%16c%17c%13c%11c4c%18c%17c-%19c%20C%21%22c%23c%24-c%25c%26C%21(C%22)c(c%19%16)c-%27c%15C%28%29c%14c-%30c%12c%10c%31c9c%32c-%33c(c%31%30)C%28(C%29)c(c%26%27)c%33c%25C%34(C%328C%34)c%24c7c6c%23c%20c5%18",
                 "c12c3c4c5c6c7c4c8c9c%10c%11c8c%12c7c%13c-%14c6c%15c%16c%14c%17c%18c%19c%20c%21=c%22c%23c%24c(c1c(c93)c%25c%24c%26c%22c%27c%21c%18c%28C%17%29C%13(c%12c%30-c%28c%27c(c%11%30)c%26c%10%25)C%31C%29C(=N)CCC%31)c%32c-%33c2c5c%15c%33c(c%19%16)C%20C%32%23",
                 "c12c3c4c5c6c7c8c5c9c%10c4c1c%11c%12c%13c2c%14c%15c3c6c%16c%15c%17c%18c%19c%16c7c%20c%19c%21c%22c%18c%23c%17c%14c%13c%24c%23c%25c%22c%26c%27c%28c%29c(c%12c%24c%28%25)-c%11c%10C%30%31C%29(c%27c-%32c(c%26%21)c%20c8c%32c%309)CC(=N)CC%31"
                 ]
        for smi in smis:
            mol = ParseSmiles(smi, False)

    def testFailures(self):
        smis = [
                "n1c[n]nn1",
                "c12ncncc1c(=C)[n]n2C",
                "c1snn2c1nc1ccccc21",
                "s1cc2nc3c(n2n1)cccc3"
                ]
        for smi in smis:
            self.assertRaises(KekulizationFailure, ParseSmiles, smi, False)

    
    def testCoverage(self):
        smis = ["c1ccn(=O)cc1", "Cs1(=O)ccccn1", # N/S-oxides
                "c1c#cccc1", # benzyne
                "n1cnc[n-]1", # from SMILES benchmark
                "[nH]1cnnn1", # ditto
                "n1cnco1",    # ditto
                "n1(cno[cH+]1)C", # ditto
                "O[p+]1(O)cccc1", # invented molecule
                "[s-]1(O)cccc1"   # invented molecule
                ]
        for smi in smis:
            mol = ParseSmiles(smi, partial=False)

class AromaticBondTest(unittest.TestCase):
    
    def testBasic(self):
        smis = ["cc", "c:c", "C:C"] # What about c=c? Is this [C]=[C]?
        for smi in smis:
            mol = ParseSmiles(smi)
            self.assertEqual(2, mol.bonds[0].order)
            self.assertTrue(mol.bonds[0].arom)

    def testRingClosures(self):
        smis = ["c1.c1", "c1.c:1", "C1.C:1"]
        for smi in smis:
            mol = ParseSmiles(smi, rulesToIgnore=16)
            self.assertEqual(2, mol.bonds[0].order)
            self.assertTrue(mol.bonds[0].arom)

class ReactionTests(unittest.TestCase):

    def testBasic(self):
        
        good = ["C>O>N", "C>>N"]
        empty = ["C>>", ">>N", ">O>", "C.>>.N"]
        bad = ["C>O>N>S"]
        for smi in good:
            mol = ParseSmiles(smi, partial=False)
        for smi in bad + empty:
            self.assertRaises(SMILESSyntaxError, ParseSmiles, smi, partial=False)
        for smi in empty:
            mol = ParseSmiles(smi, partial=False, rulesToIgnore=1)

if __name__ == "__main__":
    unittest.main()
