common_valencies = {1: {0: [1], 1: [0]},
 2: {0: [0]},
 3: {0: [1], 1: [0]},
 4: {0: [2], 1: [1], 2: [0]},
 5: {-2: [3], -1: [4], 0: [3], 1: [2], 2: [1]},
 6: {-2: [2], -1: [3], 0: [4], 1: [3], 2: [2]},
 7: {-2: [1], -1: [2], 0: [3, 5], 1: [4], 2: [3]},
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
    nbr = bond.beg if bond.beg != atom else bond.end
    return nbr.element == 7

def AdjustForAromaticity(atom):
    if atom.element == 6:
        if atom.charge == 0:
            if atom.getExplicitValence() + atom.implh == 3:
                return 1
    elif atom.element == 7:
        if atom.charge == 0:
            if atom.getExplicitValence() + atom.implh == 2:
                return 1
        elif atom.charge == 1:
            if atom.getExplicitValence() + atom.implh == 3:
                return 1
    elif atom.element == 8 or atom.element == 34: # O/Se
        if atom.charge == 1:
            if atom.getExplicitValence() + atom.implh == 2:
                return 1
    elif atom.element == 16:
        if atom.charge == 1:
            if atom.getExplicitValence() + atom.implh == 2:
                return 1
    return 0

def HasCommonValence(atom):
    data = common_valencies.get(atom.element, None) 
    # How to handle elements not in the list?
    if data is None:
        return True
        return False # Alternatively, you may wish to return True
    allowed = data.get(atom.charge, None)
    if allowed is None:
        return False # unusual charge state
    explval = atom.getExplicitValence()
    if atom.arom:
        explval += AdjustForAromaticity(atom)
    totalbonds = explval + atom.implh
    if totalbonds not in allowed:
        if not(atom.element==8 and atom.charge==0 and explval==1 and atom.implh==0 and IsAttachedToNitrogen(atom)): # TEMPO-like
            return False # unusual valence (and not TEMPO-like)
    return True

