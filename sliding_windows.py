#!env python3
# encoding=UTF8
# 该脚本用于编写 sliding windows 生成器
import multiprocessing

def sliding_window(sequence,window_size,start=0,step_size=0):
    if step_size==0:
        step_size=window_size
    yield sequence[start:start+window_size]
    start+=step_size