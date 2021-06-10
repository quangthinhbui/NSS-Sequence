import numpy as np
##############################
# Introduction
##############################
def about():
    print("=====================================================")
    print("Python Package for Sequence of Neutrosophic Soft Sets")
    print("Authors: Quang-Thinh Bui, My-Phuong Ngo, Vaclav Snasel,")
    print("         Witold Pedrycz, and Bay Vo")
    print("Version: Command Line v.1.0")
    print("Paper: The Sequence of Neutrosophic Soft Sets and A")
    print("       Decision-Making Problem in Medical Diagnosis")
    print("=====================================================")
##############################
# Presentation Matrix
##############################
def transform(s,a,w):
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
##############################
# Comparison Matrix
##############################
def matrix(s,a):
    l0 = []
    for i in range(s[1]):
        l1 = []
        for j in range(s[2]):
            m0 = 0
            m1 = 0
            m2 = 0
            for k in range(s[1]):
                if a[i][j][0] >= a[k][j][0]:
                    m0 +=1
                if a[i][j][1] >= a[k][j][1]:
                    m1 +=1
                if a[i][j][2] >= a[k][j][2]:
                    m2 +=1
                m = m0 + m1 - m2 - 1
            l1.append(m)
        l0.append(l1)
    return np.array(l0)
##############################
# Object Scores
##############################
def score(s,a):
    return np.sum(matrix(s,a),axis=1)
##############################
# Maji Algorithm on NS-Set
##############################
def maji(s,a):
    return np.argmax(np.sum(matrix(s,a),axis=1)) + 1