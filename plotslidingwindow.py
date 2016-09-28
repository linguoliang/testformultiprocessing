#codig=UTF-8
import seaborn as sb
from matplotlib import pylab
from matplotlib import pyplot as plt
import matplotlib
import numpy
import optparse
import sys
import time
import GTF_decoding
# from sliding_windows import sliding_window

_author__ = 'Guoliang Lin'
Softwarename = 'plotslidingwindow'
version = '1.0.0'
bugfixs = ''
__date__ = '2016-09-28'


seqdepth=5
def printinformations():
    print("%s software version is %s in %s" % (Softwarename, version, __date__))
    print(bugfixs)
    print('Starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))


def programends():
    print('Ends at :' + time.strftime('%Y-%m-%d %H:%M:%S'))

#
# class Loaction:
#     """store loaction for fpkm intems"""
#
#     def __init__(self, string):
#         """init values"""
#         assert isinstance(string, str)
#         self.schaffold = string.split(':')[0]
#         self.start = int(string.split(':')[1].split('-')[0])
#         self.end = int(string.split(':')[1].split('-')[1])


def _parse_args():
    """Parse the command line for options."""
    usage = 'usage: %prog -i FILE.bcf -g FILE.gff3 -o OUTPREFIX'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',
                      '--input', dest='input', type='string',
                      help='input depth file ')
    #    parser.add_option('-f','--fpkm',dest='fpkm_file',type='string',help='input fpkm file')
    #    parser.add_option('-v','--variation', dest='variation', type='string', help='input variation information file')
    parser.add_option('-g', '--gtf', dest='gtf', help='gtf file')
    parser.add_option('-o', '--output', dest='output', type='string', help='output prefix')
    parser.add_option('-w','--windows-size',dest='window_size',type='int',default=3000000,help='window size')
    options, args = parser.parse_args()
    # positional arguments are ignored
    return options


def slidingwindow(filehindle,schaffold,window_size):
    filepointer = 0
    last=[]
    data=[]
    end=window_size
    while True:
        if len(end)!=0:
            data.append(last)
            last=[]
        elment=filehindle.readline()
        if not elment:
            yield data
            raise StopIteration()
        # filepointer=filehindle.tell()
        elment=elment.strip()
        tmp=elment.split('\t')
        if tmp[0]!=schaffold:
            yield data
            data=[]
            filehindle.seek(filepointer)
            # raise StopIteration()
        else:
            data.append(tmp[1:3])
        filepointer=filehindle.tell()
        # end+=window_size
        if int(data[-1][0])>window_size:
            last=data.pop()
            yield data
            data=[]
            end+=window_size


def plotfig(windowdata):
    data=numpy.array(windowdata)
    fig=plt.bar(data[:,0],data[:,1],linewidth=0,width=1,color='r')
    fig[0].set_color("r")
    plt.show()



def mian():
    options=_parse_args()
    GTF_decoding.decodegff(options.gtf)
    with open(options.input) as filehindle:
        for key in GTF_decoding.genomeDict.keys():
            itrator=slidingwindow(filehindle,key,options.window_size)
            for windowdata in itrator:
                plotfig(windowdata)



if __name__ == '__main__':
    mian()