# encoding=UTF-8
# 这个脚本用于测试并行度

from multiprocessing.dummy import Pool as Pl
from multiprocessing import Array
import time

global list1
list1 = []


def fileinputiterator(filehandle):
    multinum = []
    conter = 0
    for elment in filehandle:
        elment.strip()
        # conter+=1
        # multinum.append(elment)
        # if conter==10:
        #     yield multinum
        #     conter=0
        #     multinum=[]
        yield elment


def printnumbers(element):
    global list1
    list1.append(element)
    pass


if __name__ == "__main__":
    t1 = time.time()
    with Pl(1) as pool:
        with open("testdata.txt") as ipf:
            # file=fileinputiterator(ipf)
            pool.map(printnumbers, ipf)
            # for element in ipf:
            #     printnumbers(element)
    t2 = time.time()
    print(t2 - t1)
    print(len(list1))
