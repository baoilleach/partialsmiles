from . import valence

# Here are some true examples of molecules that cannot be kekulized:
# n1c[n]nn1
# c12ncncc1c(=C)[n]n2C
# c1snn2c1nc1ccccc21
# n1c2c(ns1)cc1c(c2)nsn1
# s1cc2nc3c(n2n1)cccc3

def nodeIterator(degrees):
    # Iterate over degree 2 nodes first
    for idx, degree in enumerate(degrees):
        if degree == 2:
            yield idx
    for idx, degree in enumerate(degrees):
        if degree > 2:
            yield idx

# Set to True to force the greedy match to fail more often, so that
# the BackTracking search gets called more often. Useful in debugging.
CREATE_GREEDY_MATCH_FAILURES = False

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
        self.doubleBonds = [False]*len(self.mol.bonds)

        finished = False
        firsttime = True
        while(True): # Main loop

            if not firsttime or not CREATE_GREEDY_MATCH_FAILURES:
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

            firsttime = False
            if not any(self.needs_dbl_bond):
                finished = True
                break

            # Now handle any remaining degree 2 (first) or 3 nodes.
            # Once a double-bond is added that generates more degree 1 nodes
            # exit the iterator.
            change = False
            for atomIdx in nodeIterator(degrees):
                if not self.needs_dbl_bond[atomIdx]: continue
                atom = self.mol.atoms[atomIdx]
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

    # The isDoubleBond alternates between double and single, as we need to find
    # an alternating path
    def FindPath(self, atomIdx, isDoubleBond, visited, m_path):
        if self.needs_dbl_bond[atomIdx]:
            return True
        visited[atomIdx] = True
        atom = self.mol.atoms[atomIdx]
        for bond in atom.bonds:
            if not bond.arom: continue
            nbr = bond.getNbr(atom)
            if not self.kekule_system[nbr.idx]: continue
            if self.doubleBonds[bond.idx] == isDoubleBond:
                if visited[nbr.idx]: continue
                found_path = self.FindPath(nbr.idx, not isDoubleBond, visited, m_path)
                if found_path:
                    m_path.append(nbr)
                    return True
        visited[atomIdx] = False
        return False

    def BackTrack(self):
        # With an odd number of bits, it's never going to kekulize fully, but
        # let's fill in as many as we can
        count = sum(self.needs_dbl_bond)

        total_handled = 0
        for idx, needs_dbl_bond in enumerate(self.needs_dbl_bond):
            if not needs_dbl_bond: continue
            total_handled += 1
            # If there is no additional atom available to matched this bit
            # then terminate.
            if total_handled == count:
                return False

            # Our goal is to find an alternating path to another atom
            # that also needs a double bond
            self.needs_dbl_bond[idx] = False # to avoid the trivial null path being found
            visited = [False] * len(self.mol.atoms)
            m_path = []
            found_path = self.FindPath(idx, False, visited, m_path)
            if not found_path: # Implies not kekulizable
                self.needs_dbl_bond[idx] = True # reset
                continue
            total_handled += 1
            m_path.append(self.mol.atoms[idx])
            self.needs_dbl_bond[m_path[0].idx] = False
            # Flip all of the bond orders on the path from double<--->single
            for i in range(0, len(m_path) - 1):
                bond = self.mol.getBond(m_path[i], m_path[i+1])
                self.doubleBonds[bond.idx] = i % 2 == 0

        return not any(self.needs_dbl_bond)

    def AssignDoubleBonds(self):
        for bond, is_dbl in zip(self.mol.bonds, self.doubleBonds):
            if is_dbl:
                bond.order = 2

def Kekulize(mol):
    kekulizer = Kekulizer(mol)
    success = kekulizer.GreedyMatch()
    if not success:
        success = kekulizer.BackTrack()
    kekulizer.AssignDoubleBonds()
    if not success:
        for idx, needs_dbl_bond in enumerate(kekulizer.needs_dbl_bond):
            if needs_dbl_bond:
                return idx # an atom in a system that is unkekulizable
    return None
