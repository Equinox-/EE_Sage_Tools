def parallel(*args):
    q = 0
    for z in args:
        q += 1/z
    return 1/q

def series(*args):
    q = 0
    for z in args:
        q += z
    return q