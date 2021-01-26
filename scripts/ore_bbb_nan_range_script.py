# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:48:01 2020

@author: mirey
"""
import sys
import numpy as np
import pandas as pd
from math import floor, log2
from itertools import combinations # for ki_combo

def npPosG(n1, n2):
    xorval = np.bitwise_xor(n1, n2)
    pos = floor(log2(xorval))+1
    return pos

def unique_ki(k,i):
    ki_combo = combinations(range(k),i)
    big_pos_list = []
    for s in ki_combo: #for each set in ki_combo
        pos_list = []
        for n1, n2 in zip(s, s[1:]):
            n1n2_pos = npPosG(n1, n2)
            pos_list.append(n1n2_pos)
        if pos_list not in big_pos_list:
            big_pos_list.append(pos_list)
    return len(big_pos_list)

def nan_triangle_range(k_start, k_end, i_limit):
    df = pd.DataFrame()

    for k in range(k_start, k_end):
        for i in range(i_limit):
            if i > k:
                continue
            else:
                u = unique_ki(k+1,i+1)
                df.at[k+1, i+1] = u
                print("at ", k+1, ", ", i+1, ", u: ", u)
    return df

def main():
    ks = input("k start?")
    ke = input("k end?")
    il = input("i lim?")
    try:
        k_start = int(ks)
        k_end = int(ke)
        i_limit = int(il)
    except ValueError:
        print("cannot cast to int")

    df = nan_triangle_range(k_start, k_end, i_limit)
    df.to_csv(r'range.csv', index = None, header=True)

if __name__== "__main__":
    main()
