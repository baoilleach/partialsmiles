import unittest
from smiparser import SmilesParser, ParseSmiles

class PartialTest(unittest.TestCase):

    def testPartialValence(self):
        data = [
                ("C[", 1),
                ("C=[", 2),
                ("C=", 2),
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

if __name__ == "__main__":
    unittest.main()
