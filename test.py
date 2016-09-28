# encoding=UTF-8
# 这个脚本用于测试并行度

from multiprocessing.dummy import Pool as Pl
from multiprocessing import Array
import time
import numpy
import random
global list1
# size=[]
import plotslidingwindow
# for x in range(30000):
#     size.append([x,random.randint(1,100)])
size=[[1,5],[2,5],[3,6],[4,7],[5,7],[6,8],[7,8],[8,9],[9,9],[10,9],[11,7],[12,7],[13,6],[14,6],[17,6],[18,3],[19,7]]

plotslidingwindow.plotfig(size)
# list1 = [[1,2],[3,4],[4,6]]
# # list2=list1[:,0]
# b=numpy.array(list1)
# print(b[:,0])
# tmp = []
# if len(tmp) == 0:
#     print('haha')
#
#
# # x=tmp.split('\t')
#
# def fileinputiterator(filehandle):
#     multinum = []
#     conter = 0
#     flag = True
#     while True:
#         # for elment in filehandle:
#         elment = filehandle.readline()
#         if not elment:
#             raise StopIteration()
#         elment = elment.strip()
#         x = filehandle.tell()
#         print(conter)
#         conter += 1
#         # multinum.append(elment)
#         # if conter==10:
#         #     yield multinum
#         #     conter=0
#         #     multinum=[]
#         yield elment
#         if conter >= 20000:
#             raise StopIteration()
#
#
# def printnumbers(element):
#     global list1
#     list1.append(element)
#     print(element)
#     pass
#
#
# if __name__ == "__main__":
#     t1 = time.time()
#     with Pl(1) as pool:
#         with open("testdata.txt") as ipf:
#             file = fileinputiterator(ipf)
#             # x=pool.map(printnumbers, file)
#             for element in file:
#                 x = ipf.tell()
#                 # printnumbers(x)
#     t2 = time.time()
#     print(t2 - t1)
#     print(len(list1))
