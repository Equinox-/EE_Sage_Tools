class MealyMachine:
    def __init__(self):
        self.nodes=[]
        self.transitions=[]

    def node(self, nm):
        if nm in self.nodes:
            raise Exception("Bad node")
        self.nodes.append(nm)

    def add(self, fro, to, input, output):
        if not(fro in self.nodes):
            raise Exception("Bad from")
        if not(to in self.nodes):
            raise Exception("Bad to")
        self.transitions.append([fro, to, input, output])

    def nodemap(self):
        nm = []
        inputs=int(ceil(log(len(self.nodes)+1) / log(2)))
        for j in xrange(0, len(self.nodes)):
            nm.append([bin(j)[2:].zfill(inputs), self.nodes[j]])
        return nm

    def booltable(self):
        table = []
        # Fill table
        inputs=int(ceil(log(len(self.nodes)+1) / log(2)))
        ins=0
        outs=0
        if len(self.transitions) > 0:
            if isinstance(self.transitions[0][2], list):
                ins = len(self.transitions[0][2])
            else:
                ins = 1
            if isinstance(self.transitions[0][3], list):
                outs = len(self.transitions[0][3])
            else:
                outs = 1
        inputs += ins

        for i in xrange(0, 2**inputs):
            table.append([bin(i)[2:].zfill(inputs), "x"])

        for tr in self.transitions:
            src=self.nodes.index(tr[0])
            dest=self.nodes.index(tr[1])
            input=tr[2]
            output=tr[3]
            inval=(src << ins) | input
            outval=(dest << outs) | output
            table[inval] = [bin(inval)[2:].zfill(inputs), bin(outval)[2:].zfill(inputs-ins+outs)]
        return table


KARNAUGH_MAP_DECODE=["00", "01", "11", "10"]

def karnaughMap(btable):
    if (len(btable[0][0]) > 4 or len(btable[0][0]) <= 1):
        raise Exception()

    cbits=ceil(len(btable[0][0]) / 2)
    rbits=len(btable[0][0]) - cbits
    cols=2**cbits
    rows=2**rbits
    otable = [[""]]
    # Gen ROW 1
    for j in xrange(0, cols):
        otable[0].append(KARNAUGH_MAP_DECODE[j][2-cbits:2])
    for i in xrange(0, rows):
        kip=[KARNAUGH_MAP_DECODE[i][2-rbits:2]]
        for j in xrange(0, cols):
            kip.append("x")
        otable.append(kip)

    for item in btable:
        addr=item[0]
        val=item[1]
        col=addr[0:cbits]
        row=addr[cbits:]
        colV=int(col, 2)
        rowV=int(row, 2)
        if colV >= 2:
            colV ^^= 1
        if rowV >= 2:
            rowV ^^= 1
        otable[rowV+1][colV+1]=val
    return otable

def unpackKMapLayers(kmap):
    mlayer = 0
    for r in kmap:
        for c in r:
            mlayer = max(mlayer, len(c))

    layers=[None]*mlayer
    for l in xrange(0, mlayer):
        layers[l] = deepcopy(kmap)
        for i in xrange(1, len(layers[l])):
            r=layers[l][i]
            for j in xrange(1, len(r)):
                idx=min(l, len(r[j])-1)
                r[j] = r[j][idx:idx+1]
    return layers