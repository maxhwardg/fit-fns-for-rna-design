from random import randrange
import vienna

RNA_ALPHA = "AUGC"


def valid_pair(a, b):
    if a > b:
        a, b = b, a
    if a == 'A':
        return b == 'U'
    elif a == 'C':
        return b == 'G'
    elif a == 'G':
        return b == 'C' or b == 'U'
    return False


def matching_to_db(match):
    strn = []
    for i in range(len(match)):
        if match[i] < i:
            strn.append(')')
        elif match[i] > i:
            strn.append('(')
        else:
            strn.append('.')
    return ''.join(strn)


def db_to_matching(db):
    stk = []
    match = [i for i in range(len(db))]
    for i in range(len(db)):
        if db[i] == '(':
            stk.append(i)
        elif db[i] == ')':
            j = stk.pop()
            match[j] = i
            match[i] = j
    return match


def random_primary(sz):
    bases = []
    for _ in range(sz):
        bases.append(RNA_ALPHA[randrange(0, len(RNA_ALPHA))])
    return ''.join(bases)


def encode_db(db, padding):
    enc = []
    for b in db:
        idx = 0
        if b == '(':
            idx = 1
        elif b == ')':
            idx = 2
        else:
            assert(b == '.')
        hot = [0]*4
        hot[idx] = 1
        enc.append(hot)
    for _ in range(padding-len(db)):
        hot = [0]*4
        hot[3] = 1
        enc.append(hot)
    return enc


def decode_db(encoding):
    char_list = []
    for lst in encoding:
        ch = '.'
        if lst[1] == 1:
            ch = '('
        elif lst[2] == 1:
            ch = ')'
        elif lst[3] == 1:
            break
        else:
            assert(lst[0] == 1)
        char_list.append(ch)
    return ''.join(char_list)


def encode_primary(pri, padding):
    enc = []
    for b in pri:
        idx = 4
        if b in RNA_ALPHA:
            idx = RNA_ALPHA.index(b)
        hot = [0]*6
        hot[idx] = 1
        enc.append(hot)
    for _ in range(padding-len(pri)):
        hot = [0]*6
        hot[5] = 1
        enc.append(hot)
    return enc


def decode_primary(encoding):
    char_list = []
    for lst in encoding:
        lst = list(lst)
        idx = lst.index(max(lst))
        if idx == 5:
            break
        elif idx == 4:
            char_list.append('X')
        else:
            char_list.append(RNA_ALPHA[lst.index(max(lst))])
    return ''.join(char_list)


def encode_seq_struc(prim, match, padding):
    one_hot = []
    for i in range(len(prim)):
        bid = 4
        if prim[i] != 'X':
            bid = RNA_ALPHA.index(prim[i])
        dbid = 0
        if match[i] < i:
            dbid = 1
        elif match[i] > i:
            dbid = 2
        hot = [0]*16
        hot[bid*3+dbid] = 1
        one_hot.append(hot)
    while len(one_hot) < padding:
        encoding = [0]*16
        encoding[15] = 1
        one_hot.append(encoding)
    return one_hot

def encode_seq_struc_probs(prim, match, padding):
    one_hot = []
    ctx = vienna.ViennaContext(prim)
    bppt = ctx.make_bppt()
    for i in range(len(prim)):
        bid = 4
        if prim[i] != 'X':
            bid = RNA_ALPHA.index(prim[i])
        dbid = 0
        if match[i] < i:
            dbid = 1
        elif match[i] > i:
            dbid = 2
        hot = [0]*17
        hot[bid*3+dbid] = 1
        hot[16] = 1-bppt[i,i]
        one_hot.append(hot)
    while len(one_hot) < padding:
        encoding = [0]*17
        encoding[15] = 1
        one_hot.append(encoding)
    return one_hot


def decode_seq_struc(encoding,):
    pri = []
    db = []
    for oh in encoding:
        oh = list(oh)
        idx = oh.index(1)
        if idx == 15:
            break
        bid = idx//3
        dbid = idx % 3
        if bid == 4:
            pri.append('X')
        else:
            pri.append(RNA_ALPHA[bid])
        db.append(".)("[dbid])
    return ''.join(pri), ''.join(db)

def valid_pairs(base):
    if base == 'A':
        return ['U']
    elif base == 'U':
        return ['A', 'G']
    elif base == 'G':
        return ['C', 'U']
    elif base == 'C':
        return ['G']
    else:
        return []
