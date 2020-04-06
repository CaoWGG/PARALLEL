#!/usr/bin/python

# Code to generate fractal slicings of X-Y plane
# Slicing represented by tree

# File format for tree
# Each node is line in file of the form:
# Type NodeNumber Width Height leftX upperY Degree Child[1] ... Child[Degree]
# Type is 'H' (horizontal) 'V' (vertical), 'L' (leaf)
# Ordered in file by descending node ID (such that Children occur before node)
# (and therefore root node is the final node)

import datetime
import getopt
import string
import sys

import rutil

def usage(name):
    print("%s: [-h][-g] [-X] [-e EXP] [-d DEXP]  [-W WIDTH] [-H HEIGHT] [-n N] [-s SEED] [-u X:Y] [-o OUT]" % name)
    print("  -h        Print this message")
    print("  -g        Write grid diagram, rather than tree file")
    print("  -X        Generate a Hilbert Tree")
    print("  -e EXP    Set area exponent in computing weight of node for splitting")
    print("  -d DEXP   Set distance exponent in computing weight of node for splitting")
    print("  -W WIDTH  Width of grid")
    print("  -H HEIGHT Height of grid")
    print("  -n N      Minimum number of leaf nodes")
    print("  -s SEED   Specify random seed")
    print("  -u X:Y    Generate uniform grid of X x Y regions")
    print("  -o OUT    Output file")

def trim(s):
    while len(s) > 0 and s[-1] in '\n\r':
        s = s[:-1]
    return s

def showComments(lineList, outfile = sys.stdout):
    for line in lineList:
        outfile.write("# %s\n" % line)

def errorMessage(line):
    if len(line) > 0 and line[-1] not in '\r\n':
        line += '\n'
    sys.stderr.write("Error:" + line)

def infoMessage(line):
    if len(line) > 0 and line[-1] not in '\r\n':
        line += '\n'
    sys.stderr.write(line)

# Check whether string is a comment
def isComment(s):
    # Strip off leading whitespace
    while len(s) > 0 and s[0] in string.whitespace:
        s = s[1:]
    return len(s) == 0 or s[0] == '#'

class FractalTree:
    # Instance variables
    areaExponent = 0.0
    distanceExponent = 0.0
    root = None
    rng = None
    seed = 1
    nodeCount = 0
    leafCount = 0
    uniformGrid = None

    def __init__(self, areaExponent = None, distanceExponent = None):
        self.clear()
        if areaExponent is not None:
            self.areaExponent = areaExponent
        if distanceExponent is not None:
            self.distanceExponent = distanceExponent
        self.root = FractalNode(self, width = 0, height = 0)
            
    def clear(self):
        self.rng = None
        self.seed = 1
        self.nodeCount = 0
        self.leafCount = 0
        self.uniformGrid = None
        self.leafCount = 0

    def finish(self):
        self.nodeCount  = self.root.nodeCount()
        self.leafCount = self.root.leafCount()

    def generateBasic(self, width, height):
        self.clear()
        self.root = FractalNode(self, width = width, height = height)
        self.finish()
        return self

    def generateTree(self, width, height, targetCount, seed = None, hilbert = False):
        self.clear()
        if seed is not None:
            self.seed = seed
        self.rng = rutil.RNG([seed])
        self.root = HilbertNode(self, width, height) if hilbert else FractalNode(self, width, height)
        leafSet = [self.root]
        while len(leafSet) > 0 and len(leafSet) < targetCount:
            weights = [n.weight(self.areaExponent, self.distanceExponent) for n in leafSet]
            idx = self.rng.weightedIndex(weights)
            node = leafSet[idx]
            leafSet = leafSet[:idx] + leafSet[idx+1:]
            if node.branch():
                leafSet += node.children
        if len(leafSet) == 0:
            print("Failed to generate tree")
            self.clear()
        self.finish()
        return self

    def generateUniformTree(self, width, height, xcount, ycount):
        self.clear()
        if width % xcount != 0:
            errorMessage("Width %d not divisible by count %d" % (width, xcount))
            return self
        if height % ycount != 0:
            errorMessage("Height %d not divisible by count %d" % (height, ycount))
            return self
        self.root = FractalNode(self, width, height)
        self.root.uniformBranch(xcount, ycount)
        self.uniformGrid = (xcount, ycount)
        self.finish()
        return self

    def loadTree(self, infileName = None):
        self.clear()
        lineNumber = 0
        indexList = {}
        if infileName is None:
            infile = sys.stdin
        else:
            try:
                infile = open(infileName, 'r')
            except Exception as ex:
                errorMessage("Couldn't open file '%s' (%s)" % (infileName, str(ex)))
                return False
        ok = True
        for line in infile:
            lineNumber += 1
            if isComment(line):
                continue
            line = trim(line)
            fields = line.split()
            if len(fields) < 6:
                errorMessage("Couldn't parse fractal tree node from line #%d '%s'" % (lineNumber, line))
                ok = False
                break
            try:
                ifields = [int(s) for s in fields[1:]]
            except:
                errorMessage("Couldn't parse fractal tree node from line #%d '%s'" % (lineNumber, line))
                ok = False
                break
            type = fields[0]
            (id, width, height, leftX, upperY) = ifields[:5]
            n = FractalNode(tree = self, width = width, height = height, leftX = leftX, upperY = upperY, id = id, type = type)
            indexList[id] = n
            for cid in ifields[5:]:
                if cid in indexList:
                    n.addChild(indexList[cid])
                else:
                    errorMessage("Invalid child ID '%d' in line #%d" % (cid, lineNumber))
                    ok = False
                    break
            if not ok:
                break
        if infile != sys.stdin:
            infile.close()
        if not ok:
            self.clear()
            return False
        if 1 in indexList:
            self.root = indexList[1]
            self.finish()
            return True
        else:
            errorMessage("No root node found")
            self.clear()
            return False

    def width(self):
        return self.root.width

    def height(self):
        return self.root.height

    def newId(self):
        self.nodeCount += 1
        return self.nodeCount

    def buildIndex(self):
        return self.root.buildIndex(dict = {})

    def leafList(self):
        return self.root.leafList()

    def headerList(self):
        hlist = []
        if self.width() == 0 or self.height() == 0:
            hlist.append("Empty Tree")
        hlist.append("Nodes: %d" % self.nodeCount)
        hlist.append("Leaves: %d" % self.leafCount)
        hlist.append("Width: %d" % self.width())
        hlist.append("Height: %d" % self.height())
        if self.uniformGrid is None:
            hlist.append("Fractal grid generated with initial seed: %d" % self.seed)
        else:
            hlist.append("Uniform %d X %d grid" % (self.uniformGrid[0], self.uniformGrid[1]))
        return hlist

    def showHeader(self, infoList, outfile):
        ilist = infoList + self.headerList()
        ilist.append("")
        ilist.append("Type Id Width Height leftX upperY Children")
        showComments(ilist, outfile)
        
    def showTree(self, infoList = [], outfile = sys.stdout):
        self.showHeader(infoList, outfile)
        index = self.buildIndex()
        nlist = [index[self.nodeCount-i] for i in range(self.nodeCount)]
        for node in nlist:
            outfile.write(str(node) + '\n')

    def showGrid(self, outfile = sys.stdout):
        g = GridDrawing(self.width(), self.height())
        nodes = self.leafList()
        id = 0
        for n in nodes:
            label = str(id)
            g.drawBox(n.width, n.height, n.leftX, n.upperY, label)
            id += 1
        g.show(outfile)

class FractalNode:
    # Class variables
    branchFactors = [2, 3]
    # Instance variables
    tree = None # Containing tree
    id = 0
    type = 'L'
    width = 0
    height = 0
    leftX = 0
    upperY = 0
    children = []
    parent = None

    def __init__(self, tree, width = 1, height = 1, leftX = 0, upperY = 0, id = None, type = 'L'):
        self.tree = tree
        self.id = tree.newId() if id is None else id
        self.type = type
        self.width = width
        self.height = height
        self.leftX = leftX
        self.upperY = upperY
        self.children = []
        self.parent = None

    def addChild(self, child):
        self.children.append(child)
        child.parent = self

    # Compute weight of node based on area, for biasing selection of node to split
    def weight(self, areaExponent, distanceExponent):
        area = float(self.height * self.width)
        avalue = area ** areaExponent
        dist = (self.leftX + self.width/2.0) ** 2 + (self.upperY + self.height/2.0) ** 2
        dvalue = dist ** distanceExponent
        return avalue * dvalue

    # Attempt to create children of this node
    def branch(self):
        # See what type this should become
        if self.parent is None:
            tset = ['H', 'V']
        else:
            tset = ['H'] if self.parent.type == 'V' else ['V']
        type = self.tree.rng.randElement(tset)
        dim = self.width if type == 'V' else self.height
        splits = [d for d in FractalNode.branchFactors if dim % d == 0]
        if len(splits) == 0:
            return False # Can't split this node
        self.type = type
        degree = self.tree.rng.randElement(splits)
        if type == 'V':
            width = self.width//degree
            height = self.height
            deltaX = width
            deltaY = 0
        else:
            width = self.width
            height = self.height/degree
            deltaX = 0
            deltaY = height
        self.type == type
        for i in range(degree):
            leftX = self.leftX + i * deltaX
            upperY = self.upperY + i * deltaY
            child = FractalNode(self.tree, width, height, leftX, upperY)
            self.addChild(child)
        return True
        
    def uniformHorizontal(self, xcount):
        self.type = 'V'
        degree = xcount
        width = self.width//degree
        height = self.height
        deltaX = width
        deltaY = 0
        for i in range(degree):
            leftX = self.leftX + i * deltaX
            upperY = self.upperY + i * deltaY
            child = FractalNode(self.tree, width, height, leftX, upperY)
            self.addChild(child)

    def uniformVertical(self, ycount):
        if ycount == 1:
            self.type = 'L'
            return
        self.type = 'H'
        degree = ycount
        width = self.width
        height = self.height//degree
        deltaX = 0
        deltaY = height
        for i in range(degree):
            leftX = self.leftX + i * deltaX
            upperY = self.upperY + i * deltaY
            child = FractalNode(self.tree, width, height, leftX, upperY)
            self.addChild(child)

    def uniformBranch(self, xcount, ycount):
        if ycount == 1:
            self.uniformHorizontal(xcount)
        else:
            self.uniformVertical(ycount)
            for child in self.children:
                child.uniformHorizontal(xcount)

    def leafCount(self):
        if self.type == 'L':
            return 1
        count = 0
        for child in self.children:
            count += child.leafCount()
        return count

    def nodeCount(self):
        count = 1
        for child in self.children:
            count += child.nodeCount()
        return count

    # Build dictionary of all nodes, indexed by id
    def buildIndex(self, dict = {}):
        dict[self.id] = self
        for child in self.children:
            child.buildIndex(dict)
        return dict
        
    # Build list of all leaves 
    def leafList(self):
        if self.type == 'L':
            return [self]
        else:
            ls = []
            for c in self.children:
                ls += c.leafList()
        return ls

    # Generate string representation of this node
    def __str__(self):
        vlist = [self.type, self.id, self.width, self.height, self.leftX, self.upperY]
        for child in self.children:
            vlist.append(child.id)
        slist = [str(v) for v in vlist]
        return " ".join(slist)
        
class HT:
    (typeA, typeB, typeC, typeD,
     typeA1, typeA2, typeB1, typeB2,
     typeC1, typeC2, typeD1, typeD2) = list(range(12))
    typeNames = ['A', 'B', 'C', 'D', 'A1', 'A2', 'B1', 'B2', 'C1', 'C2', 'D1', 'D2']
    splitRight, splitLeft, splitUp, splitDown = list(range(4))

    def typeName(self, ht):
        if ht is None:
            return "None"
        return self.typeNames[ht]

class HilbertNode(FractalNode):

    # Definition of how each recursive partition is performed
    # Each entry of form (mode, child[0], child[1])
    splitRules = {
        HT.typeA:  (HT.splitRight, HT.typeA1, HT.typeA2),
        HT.typeB:  (HT.splitDown,  HT.typeB1, HT.typeB2),
        HT.typeC:  (HT.splitLeft,  HT.typeC1, HT.typeC2),
        HT.typeD:  (HT.splitUp,    HT.typeD1, HT.typeD2),

        HT.typeA1: (HT.splitUp,    HT.typeD,  HT.typeA),
        HT.typeA2: (HT.splitDown,  HT.typeA,  HT.typeB),
        HT.typeB1: (HT.splitLeft,  HT.typeC,  HT.typeB),
        HT.typeB2: (HT.splitRight, HT.typeB,  HT.typeA),
        HT.typeC1: (HT.splitDown,  HT.typeB,  HT.typeC),
        HT.typeC2: (HT.splitUp,    HT.typeC,  HT.typeD),
        HT.typeD1: (HT.splitRight, HT.typeA,  HT.typeD),
        HT.typeD2: (HT.splitLeft,  HT.typeD,  HT.typeC),
        }

    hilbertType = None

    def __init__(self, tree, width = 1, height = 1, leftX = 0, upperY = 0, id = None, type = 'L', hilbertType = None):
        FractalNode.__init__(self, tree, width=width, height=height, leftX=leftX, upperY=upperY, id=id, type = 'L')
        self.hilbertType = hilbertType

    # Compute the largest power of 2 that divides x
    def maxPower2(self, x):
        if x == 0:
            return 0
        val = 0
        while x % 2 == 0:
            val += 1
            x = x//2
        return val

    # What Hilbert type should be used at the root?
    def possibleRootTypes(self):
        wp2 = self.maxPower2(self.width)
        hp2 = self.maxPower2(self.height)
        if wp2 > hp2:
            return (HT.typeB1, HT.typeB2, HT.typeD1, HT.typeD2)
        elif wp2 < hp2:
            return (HT.typeA1, HT.typeA2, HT.typeC1, HT.typeC2)
        else:
            return (HT.typeA, HT.typeB, HT.typeC, HT.typeD)
        
    def branch(self):
        # Fresh root.  See what is possible
        if self.hilbertType is None and self.parent is None:
            tset = self.possibleRootTypes()
            if len(tset) == 0:
                return False
            self.hilbertType = self.tree.rng.randElement(tset)
        (st, ht0, ht1) = HilbertNode.splitRules[self.hilbertType]
        if st in [HT.splitRight, HT.splitLeft] and self.width % 2 != 0:
            return False
        if st in [HT.splitUp, HT.splitDown] and self.height % 2 != 0:
            return False
        if st == HT.splitLeft:
            self.type = 'V'
            w = self.width//2
            x = self.leftX + w
            y = self.upperY
            h = self.height
            c0 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht0)
            x = self.leftX
            c1 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht1)
            self.addChild(c0)
            self.addChild(c1)
        elif st == HT.splitRight:
            self.type = 'V'
            w = self.width//2
            x = self.leftX
            y = self.upperY
            h = self.height
            c0 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht0)
            x = self.leftX + w
            c1 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht1)
            self.addChild(c0)
            self.addChild(c1)
        elif st == HT.splitDown:
            self.type = 'H'
            h = self.height//2
            y = self.upperY
            x = self.leftX
            w = self.width
            c0 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht0)
            y = self.upperY + h
            c1 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht1)
            self.addChild(c0)
            self.addChild(c1)
        elif st == HT.splitUp:
            self.type = 'H'
            h = self.height//2
            y = self.upperY + h
            x = self.leftX
            w = self.width
            c0 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht0)
            y = self.upperY
            c1 = HilbertNode(self.tree, w, h, x, y, None, 'L', ht1)
            self.addChild(c0)
            self.addChild(c1)
        else:
            print("Unknown split rule %d" % st)
            return False
        return True

    def show(self):
        prefix = FractalNode.__str__(self)
        hts = HT().typeName(self.hilbertType)
        return  prefix + " (%s)" % hts



# Show partionings of Rectangle using ASCII art
class GridDrawing:
    # Class variables
    showDX = 2
    showDY = 1
    # Instance variables
    charList = []

    def __init__(self, width, height):
        self.charList = []
        swidth = width * self.showDX + 1
        sheight = height * self.showDY + 1
        for r in range(sheight):
            self.charList.append([' '] * swidth)
        
    # In the following, x and y are scaled coordinates
    def makeVertical(self, x, y):
        oldc = self.charList[y][x]
        newc = '|'
        if oldc in "-+":
            newc = '+'
        self.charList[y][x] = newc

    def makeHorizontal(self, x, y):
        oldc = self.charList[y][x]
        newc = '-'
        if oldc in "|+":
            newc = '+'
        self.charList[y][x] = newc

    def makeCross(self, x, y):
        self.charList[y][x] = '+'

    def horizontalLine(self, x, y, length):
        for l in range(length):
            self.makeHorizontal(x+l, y)

    def verticalLine(self, x, y, length):
        for l in range(length):
            self.makeVertical(x,y+l)

    # Here, the dimensions are unscaled
    def drawBox(self, width, height, leftX, upperY, label = None):
        scaleLeftX = leftX * self.showDX
        scaleUpperY   = upperY * self.showDY
        scaleRightX = (leftX+width) * self.showDX
        scaleLowerY = (upperY+height) * self.showDY
        scaleWidth = width * self.showDX + 1
        scaleHeight = height * self.showDY + 1
        self.horizontalLine(scaleLeftX, scaleUpperY, scaleWidth)
        self.horizontalLine(scaleLeftX, scaleLowerY, scaleWidth)
        self.verticalLine(scaleLeftX, scaleUpperY, scaleHeight)
        self.verticalLine(scaleRightX, scaleUpperY, scaleHeight)
        if label is not None and len(label) < self.showDX * width and self.showDY * height > 1:
            x = self.showDX * leftX+1
            clist = self.charList[upperY+1]
            for c in label:
                clist[x] = c
                x += 1

    def show(self, outfile = sys.stdout):
        for chars in self.charList:
            line = ''.join(chars)
            outfile.write(line + '\n')
    

def run(name, args):
    showGrid = False
    width = 36
    height = 36
    targetCount = 20
    seed = 418
    outfile = sys.stdout
    uniformGrid = None
    areaExponent = None
    hilbert = False
    optlist, args = getopt.getopt(args, "hgXe:d:W:H:n:s:u:o:")
    for (opt, val) in optlist:
        if opt == '-h':
            usage(name)
            return
        if opt == '-g':
            showGrid = True
        if opt == '-X':
            hilbert = True
        if opt == '-e':
            areaExponent = float(val)
        if opt == '-d':
            distanceExponent = float(val)
        elif opt == '-W':
            width = int(val)
        elif opt == '-H':
            height = int(val)
        elif opt == '-n':
            targetCount = int(val)
        elif opt == '-s':
            seed = int(val)
        elif opt == '-u':
            fields = val.split(':')
            if len(fields) >= 1:
                xcount = int(fields[0])
            if len(fields) >= 2:
                ycount = int(fields[1])
            else:
                ycount = xcount
            uniformGrid = (xcount, ycount)
        elif opt == '-o':
            try:
                outfile = open(val, 'w')
            except Exception as ex:
                errorMessage("Couldn't open output file '%s' (%s)" % (val, str(ex)))
                return

    tree = FractalTree(areaExponent = areaExponent, distanceExponent = distanceExponent)
    if uniformGrid is None:
        tree.generateTree(width, height, targetCount, seed, hilbert = hilbert)
    else:
        tree.generateUniformTree(width, height, uniformGrid[0], uniformGrid[1])
    if showGrid:
        tree.showGrid(outfile = outfile)
    else:
        infoList = []
        tgen = datetime.datetime.now()
        infoList.append("Generated %s" % tgen.ctime())
        tree.showTree(infoList = infoList, outfile = outfile)
    if outfile != sys.stdout:
        outfile.close()

if __name__ == "__main__":
    run(sys.argv[0], sys.argv[1:])
