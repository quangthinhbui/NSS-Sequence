import numpy as np
##############################
# Introduction
##############################
def about():
    print("=====================================================")
    print("Python Package for Sequence of Neutrosophic Soft Sets")
    print("Authors: Quang-Thinh Bui, My-Phuong Ngo, and Bay Vo")
    print("Version: Command Line v.1.0")
    print("Paper: The Sequence of Neutrosophic Soft Sets and Its")
    print("       Decision-Making Problem in Medical Diagnosis")
    print("=====================================================")
##############################
# Read from keyboard or file
##############################
def readfromkey():
    while True:
        try:
            m = input(" >Number of sets, objects, parameters in the sequence with structure <s,o,p>: ").split(",")
        except ValueError:
            continue
        if len(m) == 3:
            m = [int(i) for i in m]
            break
    t = []
    for i in range(m[0]):
        t1 = []
        print("  ******")
        print("  )SET {}".format(i+1))
        print("  ******")
        for j in range(m[1]):
            t2 = []
            print("    *Object [{}]:".format(j+1))
            for k in range(m[2]):
                while True:
                    try:
                        w = input("     -Parameter {} with structure <T,I,F>: ".format(k+1)).split(",")
                    except ValueError:
                        continue
                    if len(w) == 3:
                        w = [int(i) for i in w]
                        t2.append(w)
                        break
            t1.append(t2)
        t.append(t1)
    return m, np.array(t)
#########
def readfromfile(m):
    with open(m) as fm:
        content = fm.readlines()
        n = [x.strip() for x in content]
    s = [int(k) for k in n[0].split()]
    f = []
    for i in range(s[0]):
        fi = []
        for j in range(i*s[1]+1,(i+1)*s[1]+1):
            t = [float(k) for k in n[j].split()]
            fi.append(np.resize(t,(s[2],3)).tolist())
        f.append(fi)
    return s, np.array(f)
##############################
# Subsequence or supsequence?
##############################
def issubnssseq(m,a,b):
    count = 0
    for i in range(m[0]):
        for j in range(m[1]):
            for k in range(m[2]):
                if (a[i][j][k][0] <= b[i][j][k][0]) and (a[i][j][k][1] <= b[i][j][k][1]) and (a[i][j][k][2] <= b[i][j][k][2]):
                    count += 1
    if count == m[0]*m[1]*m[2]:
        return 1
    else:
        return 0
##############################
# Operations on NSS-sequences
##############################
def create(m):
    temp = []
    for i in range(m[0]):
        t1 = []
        for j in range(m[1]):
            t2 = []
            for k in range(m[2]):
                t2.append(np.array([0,0,0]))
            t1.append(t2)
        temp.append(t1)
    return np.array(temp)
#########
def complement(m,a):
    temp = create(m) + a
    for i in range(m[0]):
        for j in range(m[1]):
            for k in range(m[2]):
                temp[i][j][k][0] = a[i][j][k][2]
                temp[i][j][k][1] = 1 - a[i][j][k][1]
                temp[i][j][k][2] = a[i][j][k][0]
    return temp
#########
def intersection(m,a,b):
    temp = create(m) + a
    for i in range(m[0]):
        for j in range(m[1]):
            for k in range(m[2]):
                temp[i][j][k][0] = min(a[i][j][k][0],b[i][j][k][0])
                temp[i][j][k][1] = min(a[i][j][k][1],b[i][j][k][1])
                temp[i][j][k][2] = max(a[i][j][k][2],b[i][j][k][2])
    return temp
#########
def union(m,a,b):
    temp = create(m) + a
    for i in range(m[0]):
        for j in range(m[1]):
            for k in range(m[2]):
                temp[i][j][k][0] = max(a[i][j][k][0],b[i][j][k][0])
                temp[i][j][k][1] = max(a[i][j][k][1],b[i][j][k][1])
                temp[i][j][k][2] = min(a[i][j][k][2],b[i][j][k][2])
    return temp
#########
def difference(m,a,b):
    return intersection(m,a,complement(m,b))
##############################
# Distance on NSS-Sequences
##############################
# Manhattan distance
def manhattan_distance(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([np.abs(i) for i in step])
            d1 += d2
        d += d1
    return d
#########
# Euclidean distance
def euclid_distance(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([i**2 for i in step])
            d1 += d2
        d += np.sqrt(d1)
    return d
#########
# Lp-norm distance
def lpnorm_distance(m,p,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([i**p for i in step])
            d1 += d2
        d += d1**(1/p)
    return d
#########
# Chebyshev distance
def chebyshev_distance(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = []
        for j in range(m[1]):
            d2 = []
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2.append(max([np.abs(i) for i in step]))
            d1.append(max(d2))
        d += max(d1)
    return d
#########
# Normalized Manhattan distance
def manhattan_dist(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([np.abs(i) for i in step])
            d1 += d2
        d += d1
    return d/(m[0]*m[1]*m[2])
#########
# Normalized Euclid distance
def euclid_dist(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([i**2 for i in step])
            d1 += d2
        d += np.sqrt(d1)
    return d/(m[0]*np.sqrt(m[1]*m[2]))
#########
# Normalized Lp-norm distance
def lpnorm_dist(m,p,a,b):
    d = 0
    for i in range(m[0]):
        d1 = 0
        for j in range(m[1]):
            d2 = 0
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2 += np.mean([i**p for i in step])
            d1 += d2
        d += d1**(1/p)
    return d/(m[0]*(m[1]*m[2])**(1/p))
#########
# Normalized Chebyshev distance
def chebyshev_dist(m,a,b):
    d = 0
    for i in range(m[0]):
        d1 = []
        for j in range(m[1]):
            d2 = []
            for l in range(m[2]):
                step = np.subtract(a[i][j][l],b[i][j][l])
                d2.append(max([np.abs(i) for i in step]))
            d1.append(max(d2))
        d += max(d1)
    return d/m[0]
##############################
# Maji Algorithm on NSS-Sequences
##############################
def transform(s,a):
    l0 = []
    for i in range(s[2]):
        l1 = []
        for j in range(s[1]):
            l2 = np.array([0,0,0])
            for k in range(s[0]):
                l2 = l2 + a[k][j][i]
            n = l2/s[0]
            l1.append(n)
        l0.append(l1)
    t = []
    for i in range(s[1]):
        p = []
        for j in range(s[2]):
            p.append(l0[j][i])
        t.append(p)
    return np.around(np.array(t),2)
#########
def matrix(s,a):
    l0 = []
    for i in range(s[2]):
        l1 = []
        for j in range(s[1]):
            m0 = 0
            m1 = 0
            m2 = 0
            for k in range(s[1]):
                if a[j][i][0] >= a[k][i][0]:
                    m0 +=1
                if a[j][i][1] >= a[k][i][1]:
                    m1 +=1
                if a[j][i][2] >= a[k][i][2]:
                    m2 +=1
            m = m0 + m1 - m2 - 1
            l1.append(m)
        l0.append(l1)
    return np.array(l0).transpose()
#########
def score(s,a):
    return np.sum(matrix(s,a),axis=1)
########
def maji(s,a):
    return np.argmax(np.sum(matrix(s,a),axis=1)) + 1