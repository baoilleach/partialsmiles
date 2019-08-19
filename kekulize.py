import valence

class Kekulizer:
    def __init__(self, mol):
        self.mol = mol
    def GreedyMatch(self):
        # What atoms need a double bond? The job of kekulization is
        # to give all of these atoms a single double bond
        self.needs_dbl_bond = [valence.NeedsDblBond(atom) for atom in self.mol.atoms]
        # Make a copy of needs_dbl_bond, to restrict the traversal
        # in BackTrack()
        self.kekule_system = self.needs_dbl_bond[:]
        # Create lookup of degrees
        degrees = []
        degreeOneAtoms = []
        for needs_dbl_bond, atom in zip(self.needs_dbl_bond, self.mol.atoms):
            if not needs_dbl_bond:
                degrees.append(0)
                continue
            mdeg = 0
            for bond in atom.bonds:
                if not bond.arom: continue
                nbr = bond.getNbr(atom)
                if self.needs_dbl_bond[nbr.idx]:
                    mdeg += 1
            degrees.append(mdeg)
            if mdeg == 1:
                degreeOneAtoms.append(atom)
         
        # Location of assigned double bonds
        self.doubleBonds = [False]*len(mol.bonds)
         
        finished = False
        while(True): # Main loop

            # Complete all of the degree one nodes
            while len(degreeOneAtoms):
                atom = degreeOneAtoms.pop()
                # some nodes may already have been handled
                if not self.needs_dbl_bond[atom.idx]: continue
                for bond in atom.bonds:
                    if not bond.arom: continue
                    nbr = bond.getNbr(atom)
                    if not self.needs_dbl_bond[nbr.idx]: continue
                    # create a double bond from atom -> nbr
                    self.doubleBonds[bond.idx] = True
                    self.needs_dbl_bond[atom.idx] = False
                    self.needs_dbl_bond[nbr.idx] = False
                    # now update degree information for nbr's neighbors
                    for nbrbond in nbr.bonds:
                        if nbrbond == bond or not nbrbond.arom: continue
                        nbrnbr = nbrbond.getNbr(nbr)
                        if not self.needs_dbl_bond[nbrnbr.idx]: continue
                        degrees[nbrnbr.idx] -= 1
                        if degrees[nbrnbr.idx] == 1:
                            degreeOneAtoms.append(nbrnbr)
                    # only a single double bond can be made to atom so we
                    # can break here
                    break

            if all(x==False for x in self.needs_dbl_bond):
                finished = True
                break
                


    def BackTrack(self):
        pass
    def AssignDoubleBonds(self):
        pass

if __name__ == "__main__":
    import smiparser as sp
    mol = sp.ParseSmiles("c1[nH]ccc1", False)
    kekulizer = Kekulizer(mol)
    kekulizer.GreedyMatch()
    print(kekulizer.doubleBonds)


