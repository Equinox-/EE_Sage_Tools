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

# https://en.wikipedia.org/wiki/File:Wye-delta-2.svg
def deltaWye1(ra, rb, rc):
    return (rb*rc)/(ra+rb+rc)
def deltaWye2(ra, rb, rc):
    return (ra*rc)/(ra+rb+rc)
def deltaWye3(ra, rb, rc):
    return (ra*rb)/(ra+rb+rc)

def wyeDeltaA(r1,r2,r3):
    return (r1*r2+r2*r3+r3*r1)/r1
def wyeDeltaB(r1,r2,r3):
    return (r1*r2+r2*r3+r3*r1)/r2
def wyeDeltaC(r1,r2,r3):
    return (r1*r2+r2*r3+r3*r1)/r3


