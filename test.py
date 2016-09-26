#encoding=UTF-8
#这个脚本用于测试并行度

from multiprocessing.dummy import Pool as Pl

pool=Pl(2)


def fileinputiterator(filehandle):
    for elment in filehandle:
        yield elment




if __name__=="__main__":
    with open()