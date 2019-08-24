import unittest
from smiparser import SmilesParser, ParseSmiles
from exceptions import *

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
            self.assertRaises(SMILESValenceError, ParseSmiles, smi, False)
            self.assertRaises(SMILESValenceError, ParseSmiles, smi, True)

    def testInvalid(self):
        mol = ParseSmiles("[CH3]C", True)
        self.assertRaises(SMILESValenceError, ParseSmiles, "[CH2]C", True)
        self.assertRaises(SMILESValenceError, ParseSmiles, "[CH3]1C", True)
        self.assertRaises(SMILESValenceError, ParseSmiles, "[CH2]=1C", True)

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
        self.assertRaises(Exception, ParseSmiles, smi, False)
        ParseSmiles(smi, False, 1)

    def testRuleTwo(self):
        "Whether to allow an empty branch like C()C"
        smi = "C()C"
        self.assertRaises(Exception, ParseSmiles, smi, False)
        ParseSmiles(smi, False, 2)

    def testRuleThree(self):
        """Whether to allow an open parenthesis without a preceding
        atom like (CC)"""
        smis = ["(CC)", "C.(CC)"]
        for smi in smis:
            self.assertRaises(Exception, ParseSmiles, smi, False)
            ParseSmiles(smi, False, 4)

    def testRuleFour(self):
        """Whether to allow a dot within parentheses"""
        smis = ["C(C.C)C"]
        for smi in smis:
            self.assertRaises(Exception, ParseSmiles, smi, False)
            ParseSmiles(smi, False, 8)

    def testRuleFive(self):
        """Whether to allow bond closures across disconnected components"""
        smis = ["C1.C1"]
        for smi in smis:
            self.assertRaises(Exception, ParseSmiles, smi, False)
            ParseSmiles(smi, False, 16)

class ParserTests(unittest.TestCase):

    def check(self, good, bad):
        for smi in good:
            ParseSmiles(smi, False)
        for smi in bad:
            self.assertRaises(Exception, ParseSmiles, smi, False)

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
        bad = [".C", "C..C", "C-."]
        self.check(good, bad)

    def testBondChar(self):
        good = ["C-C#C", "C/C=C/Cl"]
        bad = ["-C", "C.-C", "C-=C", "C--C", "C-(C)"]
        self.check(good, bad)

    def testIsotope(self):
        good = ["[12CH4]"]
        bad = ["[0CH4]"]
        self.check(good, bad)

    def testBrackets(self):
        good = ["[16CH4]", "[NH3++]", "C[C@@H](C)C"]
        bad = ["[16C)", "C(C-)Cl"]
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

        mol = ParseSmiles("c1ccccc1.", False) # Will be kekulized
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

if __name__ == "__main__":
    unittest.main()
