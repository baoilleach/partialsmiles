import valence

def nodeIterator(degrees):
    # Iterate over degree 2 nodes first
    for idx, degree in enumerate(degrees):
        if degree == 2:
            yield idx
    for idx, degree in enumerate(degrees):
        if degree > 2:
            yield idx

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

            # Now handle any remaining degree 2 (first) or 3 nodes.
            # Once a double-bond is added that generates more degree 1 nodes
            # exit the iterator.
            change = False
            for atomIdx in nodeIterator(degrees):
                if not self.needs_dbl_bond[atomIdx]: continue
                # The following is almost identical to the code above for deg 1
                # atoms except for handling the variable 'change'
                for bond in atom.bonds:
                    if not bond.arom: continue
                    nbr = bond.getNbr(atom)
                    if not self.needs_dbl_bond[nbr.idx]: continue
                    # create a double bond from atom -> nbr
                    self.doubleBonds[bond.idx] = True
                    self.needs_dbl_bond[atom.idx] = False
                    self.needs_dbl_bond[nbr.idx] = False
                    # now update degree information for both atom's and nbr's neighbors
                    for ref in [atom, nbr]:
                        for nbrbond in ref.bonds:
                            if nbrbond == bond or not nbrbond.arom: continue
                            nbrnbr = nbrbond.getNbr(nbr)
                            if not self.needs_dbl_bond[nbrnbr.idx]: continue
                            degrees[nbrnbr.idx] -= 1
                            if degrees[nbrnbr.idx] == 1:
                                degreeOneAtoms.append(nbrnbr)
                                change = True
                    # only a single double bond can be made to atom so we
                    # can break here
                    break
                if change:
                    break # exit the iterator once we have actually set a double bond

            # We exit if we are finished or if no degree 2/3 nodes can be set
            if not change:
                break

        return finished

    def BackTrack(self):
        return False

    def AssignDoubleBonds(self):
        pass

if __name__ == "__main__":
    import smiparser as sp
    mol = sp.ParseSmiles("c1cccc1", False)
    kekulizer = Kekulizer(mol)
    success = kekulizer.GreedyMatch()
    if not success:
        success = kekulizer.BackTrack()
    kekulizer.AssignDoubleBonds()
    print(success)
