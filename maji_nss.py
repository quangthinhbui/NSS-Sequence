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
# Comparison Matrix
##############################
def matrix(s,a):
    l0 = []
    for i in range(s[0]):
        l1 = []
        for j in range(s[1]):
            m0 = 0
            m1 = 0
            m2 = 0
            for k in range(s[0]):
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