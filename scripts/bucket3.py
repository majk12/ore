# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:48:01 2020

@author: mirey
"""
import collections
import math
import pandas as pd
import numpy as np
import scipy
from scipy.special import perm, comb # for nPk
from sympy.functions.combinatorial.numbers import stirling

def trailing_zero(n):
    x = int(bin(n)[2:])
    count = 0
    while ((x & 1 ) == 0):
        x = x >> 1
        count += 1
    return count

def report_number(K, I):
    r = [[0]* (I + 1) for _ in range(K+1)]
    r[1][1] = 1
    for k in range(2, K+1):
        b = trailing_zero(k-1)
        for i in range(1, I+1):
            r[k][i] = r[k-1][i] + r[k-1][i-1] - r[k-1- int(math.pow(2,b))][i-1]
    return r

def c_bell(i, n):
    c = 0
    for u in range(0, i+1):
        c = c + pow(-1, (i-u)) * scipy.special.comb(i, u, exact = True) * pow(u, n)
    return c

def uni_ore_bayes(n, k):
    outputs = 0
    for c in range(n):
        if k < c+1:
            break
        outputs = outputs + (stirling(n, c+1) * math.factorial(c+1))
    inputs = pow(k, n)
    return outputs / inputs

def uni_ore_clww_bayes(n, k, M):
    outputs = 0
    for i in range(n):
        if k < i+1:
            break
        outputs = outputs + c_bell(i+1, n) * M[k][i+1]
    inputs = pow(k, n)
    return outputs / inputs

def g4(i, b, s):
    f = math.floor(i/b)
    m = i % b
    left = pow( (scipy.special.comb(s, f+1, exact = True)) , m)
    right = pow( (scipy.special.comb(s, f, exact = True)) , (b-m))
    return left * right


def big_cat4(n, k, b):
    prior = pow(k, n) # k, -n
    v = 0
    for i in range(1,min(n, k)+1):
        v = v + c_bell(i, n) * g4(i, b, k/b)
    if v <= 0 or prior <= 0:
        r = v / prior
    else:
        v_log = math.log2(v)
        prior_log = math.log2(prior)
        r_log = v_log - prior_log
        r = pow(2, r_log)
    return r

def gen_df(nl, i4l, i7l, i10l, c4l, c7l, c10l, n):
    df = pd.DataFrame(
        {'n':nl,
         'I4':i4l,
         'I7':i7l,
         'I10':i10l,
         'C4':c4l,
         'C7':c7l,
         'C10':c10l
        }
    )
    name = str(n)
    filename = 'back'+'%s.csv' % name
    df.to_csv(filename, index = None, header=True)
    return df

def df_3line(n_range, k, M): # hardcoded for b = 16, 128, 1024
    n = 0
    n_list = []
    I4_list = []
    I7_list = []
    I10_list = []
    C4_list = []
    C7_list = []
    C10_list = []
    for x in range(n_range):
        n = n + 1
        if n % 20 == 0: print("n : ", n)
        if n ==1024 or n == 2000 or n == 3000: gen_df(n_list, I4_list, I7_list, I10_list, C4_list, C7_list, C10_list, n)
        n_list.append(n)
        i4 = big_cat4(n, k, 16)
        I4_list.append(i4)
        i7 = big_cat4(n, k, 128)
        I7_list.append(i7)
        i10 = big_cat4(n, k, 1024)
        I10_list.append(i10)
        c4 = uni_ore_clww_bayes(n, 16, M)  # we plug b in for k for clww bayes
        C4_list.append(c4)
        c7 = uni_ore_clww_bayes(n, 128, M)
        C7_list.append(c7)
        c10 = uni_ore_clww_bayes(n, 1024, M)
        C10_list.append(c10)
    df = gen_df(n_list, I4_list, I7_list, I10_list, C4_list, C7_list, C10_list, n)
    return df

def main():
    M = report_number(1024, 1024)
    df = df_3line(4000, 1024, M)

if __name__== "__main__":
    main()
