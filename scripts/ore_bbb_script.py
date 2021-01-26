# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:48:01 2020

@author: mirey
"""
import numpy as np
import pandas as pd
from math import floor, log2
from itertools import combinations # for ki_combo

def npPos(n1, n2, w):
    xorval = np.bitwise_xor(n1, n2)
    pos = w - floor(log2(xorval))
    return pos

def unique_ki(k,i):
    ki_combo = combinations(range(k),i)
    big_pos_list = []
    w = (k-1).bit_length()
    for s in ki_combo: #for each set in ki_combo
        pos_list = []
        for n1, n2 in zip(s, s[1:]):
            n1n2_pos = npPos(n1, n2, w)
            pos_list.append(n1n2_pos)
        if pos_list not in big_pos_list:
            big_pos_list.append(pos_list)
    return len(big_pos_list)

def bit_pattern_fork(k, i_limit):
    zero_data = np.zeros(shape=(1,i_limit)) # 1 row for 1 k 
    df = pd.DataFrame(zero_data)
    df.rename_axis('k')

    for i in range(i_limit):
        if i > k:
            continue
        else:
            u = unique_ki(k,i+1)
            print("at ", i+1, " u = ", u)
            df.at[0, i+1] = u
    return df

def main():
    k26df = bit_pattern_fork(26, 26)
    k26df.to_csv(r'data\k26df.csv', index = None, header=True)

if __name__== "__main__":
    main()