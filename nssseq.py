import numpy as np
#########
def welcome():
    print("=====================================================")
    print("Python Package for Sequence of Neutrosophic Soft Sets")
    print("Authors: Quang-Thinh BUI, My-Phuong NGO, and Bay Vo")
    print("Version: Command Line v.1.0")
    print("Paper: The Sequence of Neutrosophic Soft Sets and Its")
    print("       Decision-Making Problem in Medical Diagnosis")
    print("=====================================================")
#########
def size():
    while True:
        try:
            m = input(" >Number of sets, objects, parameters in the sequence with structure <s,o,p>: ").split(",")
        except ValueError:
            continue
        if len(m) == 3:
            m = [int(i) for i in m]
            break
    return m
#########
def matrix(m):
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
    return np.array(t)
#########
def file(m,s):
    with open(m) as fm:
        content = fm.readlines()
        n = [x.strip() for x in content]
    f = []
    for i in range(s[0]):
        fi = []
        for j in range(i*s[1],(i+1)*s[1]):
            t = [float(k) for k in n[j].split()]
            fi.append(np.resize(t,(s[2],3)).tolist())
        f.append(fi)
    return np.array(f)
#########
def transform(a,s):
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
def compmat(a,s):
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
def score(a,s):
    return np.sum(compmat(a,s),axis=1)
########
def maji(a,s):
    return np.argmax(np.sum(compmat(a,s),axis=1)) + 1