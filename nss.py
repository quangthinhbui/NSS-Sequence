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
def compmat(a,s):
    l0 = []
    for i in range(s[1]):
        l1 = []
        for j in range(s[0]):
            m0 = 0
            m1 = 0
            m2 = 0
            for k in range(s[0]):
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
#########
def maji(a,s):
    return np.argmax(np.sum(compmat(a,s),axis=1)) + 1