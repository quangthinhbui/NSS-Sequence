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
def transform(a,s,w):
    l0 = []
    for i in range(s[2]):
        l1 = []
        for j in range(s[1]):
            l2 = np.array([0,0,0])
            for k in range(s[0]):
                m = a[k][j][i] * w[k]
                l2 = l2 + m
            n = l2/(s[0]*max(w))
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