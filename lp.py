#!/usr/bin/python
# Version 0.05alpha
# June 23, 2026
# must run: python -m pip install sympy
# No AI used here


#MIT License

#Copyright (c) 2026 InstrumentPurple(Github)

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

DEFAULT_ERROR_MARGIN=0.00000000001

from math import sqrt
import sympy
import re
import copy



def withinError(subj: float, prist:float, delta:float) -> float:
    dist = abs(subj - prist)
    if dist >= 0.0 and dist <= delta:
        return prist
    return subj

class Vec:
    def __init__(self, dim: int,undefined:bool=False ):
        self.dimen = dim
        self.vals: list[float] = []
        self.undef = undefined
        for i in range(0,dim):
            self.vals.append(0.0)

    def print(self):
        print()
        if self.undef:
            print("undefined")
        for val in self.vals:
            print(val)

    def dotp(self, rhs: "Vec") -> float:
        total = 0.0
        for (val1,val2) in zip(self.vals, rhs.vals):
            total += val1*val2
        return total

    def length(self) -> float:
        total = 0.0
        for val in self.vals:
            total += val**2.0

        return sqrt(total)

    def magnitude(self) -> float:
        return self.length()

    def norm(self) -> "Vec":
        if len(self.vals) < 1:
            u=Vec(0)
            u.undef=True
            return u

        l = self.length()

        if withinError(l,0.0, DEFAULT_ERROR_MARGIN) == 0.0:
            z=Vec(self.dimen) # zero vector
            return z

        nw = Vec(self.dimen)
        for i in range(0,self.dimen):
            nw.vals[i] = self.vals[i] / l

        return nw

    def isPerpen(self, rhs: "Vec") -> bool:
        return withinError(self.dotp(rhs),0.0,DEFAULT_ERROR_MARGIN) == 0.0

    def isUnitVec(self):
        diff = 1.0 - self.length()
        return (diff >= 0.0) and (diff < DEFAULT_ERROR_MARGIN)

    def subt(self, rhs: "Vec") -> "Vec":
        result = Vec(len(rhs.vals))
        for i in range(0,len(result.vals)):
            result.vals[i] = self.vals[i] - rhs.vals[i]

        return result

    def cross(self, rhs: "Vec"):
        result = Vec(3)
        if len(self.vals) < 3 or len(rhs.vals) < 3:
            result.undef=True
            return result

        result.vals[0] = self.vals[1] * rhs.vals[2] - self.vals[2]*rhs.vals[1]

        result.vals[1] = self.vals[2] * rhs.vals[0] - self.vals[0]*rhs.vals[2]

        result.vals[2] = self.vals[0] * rhs.vals[1] - self.vals[1] * rhs.vals[0]

        return result

    def toRowMatrix(self) -> "Matrix":
        m=Matrix(1,len(self.vals),0)
        m.values[0]=self.vals[0:]
        return m

    def toColMatrix(self) -> "Matrix":
        m=Matrix(len(self.vals),1,0)
        for c in range(0,len(self.vals)):
            m.values[c][0] = self.vals[c]
        return m

    def sympyValueOfX(self) -> dict["str",float]:
        x = "x"
        i = 0

        tmpstr=""
        cur = {}
        while i < len(self.vals):
            tmpstr = x + str(i)
            cur[tmpstr] = self.vals[i]
            i += 1

        return cur

    def add(self, rhs: "Vec") -> "Vec":
        tmp=Vec(rhs.dimen)
        if self.dimen == rhs.dimen:
            for i in range(0,self.dimen):
                tmp.vals[i] = self.vals[i] + rhs.vals[i]
        return tmp

    def scale(self, scalar: float) -> "Vec":
        scaled=copy.deepcopy(self)
        for i in range(0, self.dimen):
            scaled.vals[i] *= scalar
        return scaled

    def eq(self, rhs: "Vec") -> bool:
        evalu = True

        for (l,r) in zip(self.vals,rhs.vals):
            if withinError(l,r,DEFAULT_ERROR_MARGIN) != r:
                evalu = False
                break

        return evalu

    def resetDimension(self):
        self.dimen = len(self.vals)

class Matrix:
    def __init__(self, r: int, c: int, defaultInit: bool=True, isUndefined: bool=False):
        self.values: list[list[float]] = []
        row: list[float] = []
        self.rows=r
        self.cols=c
        self.undef=isUndefined
        for i in range(0,c):
            row.append(0.0)

        for i in range(0,r):
            self.values.append(row[0:len(row)]) # deep copy

        if defaultInit:
            if c == r: # is a square
                for i in range(0,c):
                    self.values[i][i] = 1.0


    def transpose(self) -> "Matrix":
        done = Matrix(self.rows,self.cols)
        for r in range(0,self.rows):
            tmp=self.toRowVecAt(r)
            done.setPlaceVecCol(r,tmp)
        return done

    def print(self):
        if self.undef:
            print("undefined")

        maximum = 0
        for r in range(0,self.rows):
            for c in range(0,self.cols):
                if len(str(self.values[r][c])) > maximum:
                    maximum = len(str(self.values[r][c]))

        for r in range(0,self.rows):
            for c in range(0,self.cols):
                leng = len(str(self.values[r][c]))
                delta = maximum - leng
                print(" " * delta, end="")
                print(str(self.values[r][c]) + " ", end="")
            print()

    def mul(self, rhs: "Matrix") -> "Matrix":
        dest=Matrix(self.rows, rhs.cols,0)

        curTotal=0.0
        if self.cols == rhs.rows:
            #for each entry in the dest
            for r in range(0, dest.rows):
                for c in range(0, dest.cols):
                    #calc dot product
                    curTotal = 0.0
                    for lhsc in range(0, self.cols):
                        curTotal += self.values[r][lhsc] * rhs.values[lhsc][c]
                    dest.values[r][c] = curTotal
        else:
            dest=Matrix(0,0,0,True)

        return dest

    def setPlaceVecCol(self,ident: int, col: Vec):
        for r in range(0, self.rows):
            self.values[r][ident] = col.vals[r]
        return

    def setPlaceVecRow(self,ident: int, col: Vec):
        for r in range(0, self.cols):
            self.values[ident][r] = col.vals[r]
        return

    def placeVecCol(self,ident: int, col: Vec) -> "Matrix":
        newmat = copy.deepcopy(self)
        for r in range(0, newmat.rows): #changed to -1
            #print("r",r)
            newmat.values[r][ident] = col.vals[r]
        return newmat

    def placeVecRow(self, ident: int, row: Vec) -> "Matrix":
        newmat = copy.deepcopy(self)
        for c in range(0,newmat.cols):
            newmat.values[ident][c] = row.vals[c]
        return newmat

    def toColVec(self) -> Vec:
        cur = Vec(self.rows)
        for i in range(0,self.rows):
            cur.vals[i] = self.values[i][0]
        return cur

    def toRowVec(self) -> Vec:
        v=Vec(len(self.values[0]))
        v.vals = self.values[0][0:]
        return v

    def toColVecAt(self, at: int) -> Vec:
        cur = Vec(self.rows)
        for i in range(0,self.rows):
            cur.vals[i] = self.values[i][at]
        return cur

    def toRowVecAt(self, at: int) -> Vec:
        v=Vec(len(self.values[at]))
        v.vals = self.values[at][0:]
        return v

    def load(self, path: str):
        fh = open(path, "r")
        rawData = fh.read()
        fh.close()
        rawData = re.sub(r' +', ' ', rawData)
        lines = rawData.split('\n')

        cleanLines = []
        for line in lines:
            line = re.sub("\r?", "", line)
            cleanLines.append(line)

        almost = []
        for cleanLine in cleanLines:
            if cleanLine == "":
                continue
            almost.append(cleanLine.split(' '))

        finished = []
        for row in almost:
            workingRow = []

            for elm in row:
                if elm == "":
                    continue
                workingRow.append(float(elm))
            finished.append(workingRow)

        self.values = finished
        self.rows = len(finished)
        self.cols = len(finished[0])

    def save(self, path: str):
        with open(path, "w") as f:
            for i in range(0,self.rows):
                for j in range(0, self.cols):
                    if j + 1 == len(self.values[j]):
                        f.write(str(self.values[i][j]))
                    else:
                        f.write(str(self.values[i][j]) + " ")
                f.write("\n")

    def scaleRow(self, rowNum:int, by:float):
        if not(rowNum < len(self.values)):
            self.undef = True
            return

        i = 0
        while i < len(self.values[rowNum]):
            self.values[rowNum][i] *= by
            i += 1

    def scaleSubtractRows(self,scalar:float, srcRow:int, destRow:int):
        for spot in range(0,self.cols):
           self.values[destRow][spot] -= (scalar*self.values[srcRow][spot])

    def scaleAddRows(self,scalar:float, srcRow:int, destRow:int):
        for spot in range(0,self.cols):
           self.values[destRow][spot] += (scalar*self.values[srcRow][spot])

    def rowReduce(self) -> "Matrix":
        workingMat = copy.deepcopy(self)

        # ensure non zeros in pivot entries
        found = False
        for pivot in range(0,min([workingMat.rows, workingMat.cols])):
            if withinError(workingMat.values[pivot][pivot],0.0,DEFAULT_ERROR_MARGIN) == 0.0:
                #try to find non-zero in rows below the pivot
                for i in range(pivot + 1, workingMat.rows):
                    if withinError(workingMat.values[i][pivot],0.0,DEFAULT_ERROR_MARGIN) != 0.0:
                        workingMat.scaleAddRows(1.0, i, pivot)
                        found = True
                        break
                if not found:
                    #try to find non-zero above pivots
                    j = pivot -1
                    while j >= 0:
                        if withinError(workingMat.values[j][pivot],0.0,DEFAULT_ERROR_MARGIN) != 0.0:
                            workingMat.scaleAddRows(1.0, j, pivot)
                            found = True
                            break
                        j -= 1
                    if not found: # entire column is zero
                        workingMat.undef=True
                        return workingMat
                found=False

        for pivot in range(0,min([workingMat.rows, workingMat.cols])):
            pivotVal = workingMat.values[pivot][pivot]
            if withinError(pivotVal,0.0,DEFAULT_ERROR_MARGIN) == 0.0:
                continue
            workingMat.scaleRow(pivot, 1.0/pivotVal)
            for spot in range(pivot+1,min([workingMat.rows, workingMat.cols])):
                workingMat.scaleSubtractRows(workingMat.values[spot][pivot],pivot, spot)

        for pivot in range(1,min([workingMat.rows, workingMat.cols])):
            i = pivot -1
            while i >= 0:
                workingMat.scaleSubtractRows(workingMat.values[i][pivot],pivot,i)
                i -= 1

        return workingMat

    def niceRowReduce(self) -> "Matrix":
        dest = copy.deepcopy(self)
        dest = dest.rowReduce()
        return dest.closeInts()

    def findx(self) -> Vec:
        got = self.rowReduce()
        #if got.rows+1 != got.cols:
        #    ve=Vec(got.rows) # if you don't check undef you might have an invalid array access if set to 0
        #    ve.undef=True
        #    return ve
        #got.print()
        return got.toColVecAt(got.cols-1)

    def rank(self) -> int:
        workingMat = self.rowReduce()
        count = 0
        for i in range(0,min([workingMat.rows, workingMat.cols])):
            if withinError(workingMat.values[i][i],1.0,DEFAULT_ERROR_MARGIN) == 1.0:
                count += 1
            else:
                break
        return count

    def scale(self, scalar: float) -> "Matrix":
        dest = copy.deepcopy(self)
        for i in range(0, dest.rows):
            dest.scaleRow(i, scalar)
        return dest

    def eq(self, rhs: "Matrix") -> bool:
        if self.cols !=  rhs.cols or self.rows != rhs.rows:
            return False
        for rowl,rowr in zip(self.values, rhs.values):
            if len(rowl) != len(rowr):
                return False
            for (l,r) in zip(rowl,rowr):
                if withinError(l,r,DEFAULT_ERROR_MARGIN) != r:
                    return False
        return True

    def isSymetric(self):
        t = self.transpose()
        return self.eq(t)

    def add(self, rhs: "Matrix") -> "Matrix":
        dest = Matrix(rhs.rows, rhs.cols)
        if rhs.cols != self.cols or rhs.rows != self.rows:
            return rhs

        i,j=0,0
        while i < self.cols:
            while j < self.rows:
                dest.values[j][i] = self.values[j][i] + rhs.values[j][i]
                j += 1
            i += 1
            j=0

        return dest

    def inverse(self) -> "Matrix":
        dest = Matrix(self.rows, self.cols*2)

        i,j=0,0

        while i < self.rows:
            while j < self.cols:
                dest.values[i][j] = self.values[i][j]
                j += 1
            j=0
            i += 1


        i=0
        j=self.cols

        while i < self.rows and j < (self.cols*2):
            dest.values[i][j]=1.0

            i += 1
            j += 1

        done = dest.rowReduce()
        dest2 = Matrix(self.rows, self.cols)

        i = self.cols
        j = 0
        while j < done.rows:
            while i < done.cols:
                dest2.values[j][i-self.cols] = done.values[j][i]
                i += 1
            i = 0
            j += 1

        return dest2

    # shave ints
    def closeInts(self) -> "Matrix":
        dest=Matrix(self.rows,self.cols)

        for i in range(0,self.rows):
            for j in range(0,self.cols):
                intg = int(self.values[i][j])
                delta = self.values[i][j] - float(intg)
                if delta < DEFAULT_ERROR_MARGIN and delta > -DEFAULT_ERROR_MARGIN:
                    dest.values[i][j] = float(intg)
                else:
                    dest.values[i][j] = self.values[i][j]

        return dest

    def trace(self) -> float:
        total = 0.0
        for val in range(0,min(self.rows, self.cols)):
            total += self.values[val][val]
        return total

def pointBetween(a: Vec, b: Vec, lamb: float) -> Vec:
    dest=Vec(a.dimen)
    if lamb >= 0.0 and lamb <= 1.0:
        aprime = a.scale(lamb)
        bprime = b.scale(1.0 - lamb)
        dest = aprime.add(bprime)
        return dest
    dest.undef = True
    return dest


#I don't know how usefull it is but it's something
class LinearProgram:
    def __init__(self, objective: sympy.Mul, restrictions: list[sympy.LessThan]):
        self.obj = objective
        self.rest = restrictions
        self.addedSlack = False


    def feasible(self, xVals: dict[str,float]) -> bool:
        #enforce non-negativity
        stillNonNeg = True
        for (s,f) in xVals.items():
            stillNonNeg = stillNonNeg and (f >= 0.0)

        if not stillNonNeg:
            return False

        #check user inputed ineqs
        stillGood = True
        for rst in self.rest:
            stillGood = stillGood and rst.subs(xVals)

        return stillGood

    def maxSolution(self, allVert: list[dict[str,float]]) -> tuple[float, dict[str,float]]:
        solved: list[tuple[float, dict[str,float]]] = []

        for vert in allVert:
            solved.append( (self.obj.subs(vert),vert) )

        maximum=solved[0]
        for s in solved:
            if s[0] >= maximum[0]:
                maximum=s

        return maximum

    def addSlackVars(self): # really should be called rewrite to cannonical form
        pre = "s"
        count = 0

        equalities: list[sympy.Eq] = []

        for ineq in self.rest:
            equalities.append(sympy.Eq(ineq.lhs + sympy.Symbol(pre + str(count)), ineq.rhs))
            count += 1

        self.rest = equalities
        self.addedSlack = True

    # NOT COMPLETED
    def simplexMin(self) -> Vec:
        pass
        #if not self.addedSlack:
        #    self.addSlackVars()
        # turn into a tablue represented by a Matrix object
        #mat=Matrix(len(self.rest)+1,len(self.obj))

        #d = self.obj.as_coefficients_dict()
        #cur = 0
        #for term, coeff in d.items():
        #    mat[len(mat.values)-1][cur] = coeff
        #    cur += 1





# sets all free variables to zero when len(points) > rows
def AffineComb(*points: Vec) -> Vec:
    destMat = Matrix(len(points[0].vals),len(points)-1)

    track=0
    if len(points)+1 == len(points[0].vals):
            for point in points:
                destMat.setPlaceVecCol(track, point)
            solus = destMat.findx()
            if sumsToOne(solus.vals):
                return solus
            else:
                tmp=Vec(len(points[0].vals))
                tmp.undef=True
                return tmp


    # place translated vectors V_2 through V_N into a matrix
    track=1
    while track < len(points):
        destMat.setPlaceVecCol(track-1,points[track].subt(points[0]))
        track += 1

    solus = destMat.findx()
    # fill in the last coeffient for V_1
    c_1=0.0
    for c_k in solus.vals:
        c_1 -= c_k
    c_1 += 1

    #pour the solutions and c_1 into the same vector
    destVec = Vec(len(solus.vals)+1)
    destVec.vals[0] = c_1

    track=1
    while track < len(destVec.vals):
        destVec.vals[track]= solus.vals[track-1]
        track += 1

    return destVec

def sumsToOne(coeffs: list[float]) -> bool:
    total = 0.0
    for ce in coeffs:
        total += ce

    return withinError(total,1.0,DEFAULT_ERROR_MARGIN) == 1.0

# for Matrix(n,n+1)
def checkAffine(aug: Matrix) -> bool:
    coeffs = aug.findx()
    return sumsToOne(coeffs.vals)

def convertColVecs(*vecs) -> Matrix:
    dest=Matrix(len(vecs[0].vals),len(vecs))
    for i in range(0,dest.cols):
        dest.setPlaceVecCol(i, vecs[i])
    return dest

def convertRowVecs(*vecs) -> Matrix:
    dest=Matrix(len(vecs),len(vecs[0].vals))
    for i in range(0,dest.rows):
        dest.setPlaceVecRow(i, vecs[i])
    return dest

def convertMatrixToColVecs(source: Matrix) -> list[Vec]:
    dest = []

    for i in range(0, source.cols):
        dest.append(source.toColVecAt(i))

    return dest


if __name__=="__main__":
    print("Begining tests / examples... ")

    I = Matrix(2,2)
    I.values=[[1,1],[0,2]]
    P = Matrix(2,1)
    P.values=[[3],[8]]

    I.mul(P).print()

    print()
    g=Matrix(2,2)
    g.values = [[1.5,1.3],[1.6, 0.4]]

    g.print()
    print()
    g.transpose().print()

    print()
    print("rank of g: ", end="")
    print(g.rank())



    m=Vec(2)
    m.vals=[5.0,0.0]

    m_ = m.toColMatrix()
    Ax=g.mul(m_) # m_ is rhs because it is x in Ax


    n=Vec(2)
    n.vals=[1.3,1.0]
    n.print()
    n.scale(2.0).print()
    n.add(m).print()


    #A=Matrix(0,0)
    #load automatically adjusts the dimensions to the matrix in the text file
    #A.load('./testmat2.txt')
    #A.print()
    #print()

    #ve=A.findx()
    #ve.print()

    #ve2=A.rowReduce()
    #ve2.print()

    #Q=Matrix(0,0)
    #Q.load("./testrank.txt")
    #print("rank of Q: ", end="")
    #print(Q.rank())

    #H=Matrix(0,0)
    #H.load("./testmat3.txt")
    #H.rowReduce().print()


    x0 = sympy.Symbol('x0')
    x1 = sympy.Symbol('x1')

    #objective function set up
    y = x0*8.5 - x1/3.0

    # restrictions
    k = (5.9*x0 + 2.2*x1 <= 260)
    k2 = (7.8*x0 + 3.1*x1 <= 100)

    lp = LinearProgram(y,[k,k2])

    # a and b are possible solutions to the objective function
    a = Vec(2)
    a.vals=[1.0,2.0]

    b=Vec(2)
    b.vals=[3.0,8.0]

    print(lp.feasible(a.sympyValueOfX()))
    print(lp.feasible(b.sympyValueOfX()))
    print("maximum tested solution")
    print(lp.maxSolution([a.sympyValueOfX(),b.sympyValueOfX()]))
    lp.addSlackVars()

    u=Vec(2)
    u.vals=[1,2]
    v=Vec(2)
    v.vals=[-2,2]
    w=Vec(2)
    w.vals=[0,4]
    p=Vec(2)
    p.vals=[3,7]
    q=Vec(2)
    q.vals=[5,3]
    aff=AffineComb(u,v,w,p,q)
    print(aff.vals)

    V = [Vec(3),Vec(3),Vec(3),Vec(3)]
    V[0].vals = [-3.0,1.0,1.0]
    V[1].vals = [0.0,4.0,-2.0]
    V[2].vals = [4.0,-2.0,6.0]
    V[3].vals = [17.0,1.0,5.0]
    aff=AffineComb(*V)
    print(aff.vals)

    print("check affine")
    g=Matrix(3,4)
    g.values=[[1,1,4,1],[2,6,4,1],[1,1,1,1]]
    print(checkAffine(g))

    term = Matrix(2,2)
    term2 = Matrix(2,2)
    term2.values = [[1,1],[1,1]]
    s = term.add(term2)
    s.print()

    print()

    A=Matrix(2,2)
    A.values=[[1,2],[7,9]]
    invr = A.inverse()
    A.mul(invr).closeInts().print() # should see identity

    print()
    U=Matrix(2,2)
    U.values=[[1.2,41.000000000001],[-500.0000000000034,5]]
    U.closeInts().print()

    print()

    U.niceRowReduce().print() # should see identity again since U is invertible

    print()

    L=pointBetween(V[0],V[1],0.5)
    print("exactly halfway between V[0] and V[1]")
    L.print()

    print()

    cv = convertColVecs(*V)
    cv.print()

    print()

    rv = convertRowVecs(*V)
    rv.print()

    # defaultInit controls how the matrix is initialized if you want all 0s
    # You must have False as your third argument. it's diagonal will be 1s otherwise
    z = Matrix(2,2,False)
    z.print()

    o = Matrix(2,2,True)
    o.print()

    # it only produces identity when it is square
    ni = Matrix(2,3)
    ni.print() # non-square

    print("Is identity symetric?")
    symet = Matrix(20,20,True)
    print(symet.isSymetric())

    print("trace of identity of size 20")
    print(symet.trace())

    print()

    got = convertMatrixToColVecs(o)
    got[0].print()
    got[1].print()
    #o.save("./test.txt")
    #B = Matrix(0,0)
    #B.load("./test.txt")
    #B.scale(20).print()
