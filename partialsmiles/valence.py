common_valencies = {1: {0: [1], 1: [0]},
 2: {0: [0]},
 3: {0: [1], 1: [0]},
 4: {0: [2], 1: [1], 2: [0]},
 5: {-2: [3], -1: [4], 0: [3], 1: [2], 2: [1]},
 6: {-2: [2], -1: [3], 0: [4], 1: [3], 2: [2]},
 7: {-2: [1], -1: [2], 0: [3], 1: [4], 2: [3]}, # exclude 5-valent nitrogens
 8: {-2: [0], -1: [1], 0: [2], 1: [3, 5]},
 9: {-1: [0], 0: [1], 1: [2], 2: [3, 5]},
 10: {0: [0]},
 11: {-1: [0], 0: [1], 1: [0]},
 12: {0: [2], 2: [0]},
 13: {-2: [3, 5], -1: [4], 0: [3], 1: [2], 2: [1], 3: [0]},
 14: {-2: [2], -1: [3, 5], 0: [4], 1: [3], 2: [2]},
 15: {-2: [1, 3, 5, 7], -1: [2, 4, 6], 0: [3, 5], 1: [4], 2: [3]},
 16: {-2: [0], -1: [1, 3, 5, 7], 0: [2, 4, 6], 1: [3, 5], 2: [4]},
 17: {-1: [0], 0: [1, 3, 5, 7], 1: [2, 4, 6], 2: [3, 5]},
 18: {0: [0]},
 19: {-1: [0], 0: [1], 1: [0]},
 20: {0: [2], 1: [1], 2: [0]},
 31: {-2: [3, 5], -1: [4], 0: [3], 1: [0], 2: [1], 3: [0]},
 32: {-2: [2, 4, 6], -1: [3, 5], 0: [4], 1: [3], 4: [0]},
 33: {-3: [0], -2: [1, 3, 5, 7], -1: [2, 4, 6], 0: [3, 5], 1: [4], 2: [3]},
 34: {-2: [0], -1: [1, 3, 5, 7], 0: [2, 4, 6], 1: [3, 5], 2: [4]},
 35: {-1: [0], 0: [1, 3, 5, 7], 1: [2, 4, 6], 2: [3, 5]},
 36: {0: [0, 2]},
 37: {-1: [0], 0: [1], 1: [0]},
 38: {0: [2], 1: [1], 2: [0]},
 49: {-2: [3, 5], -1: [2, 4], 0: [3], 1: [0], 2: [1], 3: [0]},
 50: {-2: [2, 4, 6], -1: [3, 5], 0: [2, 4], 1: [3], 2: [0], 4: [0]},
 51: {-2: [1, 3, 5, 7], -1: [2, 4, 6], 0: [3, 5], 1: [2, 4], 2: [3], 3: [0]},
 52: {-2: [0], -1: [1, 3, 5, 7], 0: [2, 4, 6], 1: [3, 5], 2: [2, 4]},
 53: {-1: [0], 0: [1, 3, 5, 7], 1: [2, 4, 6], 2: [3, 5]},
 54: {0: [0, 2, 4, 6, 8]},
 55: {-1: [0], 0: [1], 1: [0]},
 56: {0: [2], 1: [1], 2: [0]},
 81: {0: [1, 3]},
 82: {-2: [2, 4, 6], -1: [3, 5], 0: [2, 4], 1: [3], 2: [0]},
 83: {-2: [1, 3, 5, 7], -1: [2, 4, 6], 0: [3, 5], 1: [2, 4], 2: [3], 3: [0]},
 84: {0: [2, 4, 6]},
 85: {-1: [0], 0: [1, 3, 5, 7], 1: [2, 4, 6], 2: [3, 5]},
 86: {0: [0, 2, 4, 6, 8]},
 87: {0: [1], 1: [0]},
 88: {0: [2], 1: [1], 2: [0]}}

def IsAttachedToNitrogen(atom):
    bond = atom.bonds[0]
    nbr = bond.getNbr(atom)
    return nbr.element == 7

def IsSpecialCase(atom):
    elem = atom.element
    if elem == 7:
        # Any exo-cyclic double bond from a N
        # e.g. pyridine N-oxide as the double bond form
        if atom.getExplicitDegree() + atom.implh == 3 and atom.charge == 0:
            return True
    elif elem == 16:
        # e.g. Cs1(=O)ccccn1
        if atom.getExplicitDegree() + atom.implh == 4 and atom.charge == 0:
            return True
    return False

def NeedsDblBond(atom):
    # Does it have an aromatic bond?
    if all(not bond.arom for bond in atom.bonds):
        return False

    # Does it already have an explicit double bond?
    for bond in atom.bonds:
        if bond.arom: continue
        nbr = bond.getNbr(atom)
        bo = bond.order
        if bo < 2: continue
        if bo == 2:
            return True if IsSpecialCase(atom) else False
        else: # bo > 2
            return False

    # Is it one of the cases where we know that it only has single bonds?
    chg = atom.charge
    deg = atom.getExplicitDegree() + atom.implh
    elem = atom.element
    if elem == 6:
        if deg == 3 and (chg == 1 or chg == -1):
            return False
    elif elem in [5, 7, 15, 33, 51, 83]:
        if chg == 0: # e.g. a pyrrole-type nitrogen
            if deg == 3 or deg > 4:
                return False
        elif chg == -1:
            if deg == 2:
                return False
        elif chg == 1:
            if deg > 3:
                return False
    elif elem in [8, 16, 34, 52]:
        if chg == 0:
            if deg == 2 or deg == 4 or deg > 5:
                return False
        elif chg == -1 or chg == 1:
            if deg == 3 or deg == 5 or deg > 6:
                return False

    return True # It needs a double bond
