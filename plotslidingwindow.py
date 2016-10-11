# codig=UTF-8
import seaborn as sb
from matplotlib import pylab
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import numpy
import optparse
import sys
import time
import GTF_decoding
import random

# from sliding_windows import sliding_window

_author__ = 'Guoliang Lin'
Softwarename = 'plotslidingwindow'
version = '1.0.0'
bugfixs = ''
__date__ = '2016-09-28'

seqdepth = 5


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
    parser.add_option('-w', '--windows-size', dest='window_size', type='int', default=3000000, help='window size')
    options, args = parser.parse_args()
    # positional arguments are ignored
    return options


def slidingwindow(filehindle, window_size):  # generate sliding_window data

    # init all parameters
    start = 0
    end = window_size
    elment = filehindle.readline()
    if not elment:
        raise StopIteration()
    elment = elment.strip()
    tmp = elment.split('\t')
    scaffold = tmp[0]
    data = [[int(tmp[1]), int(tmp[2])]]

    while True:
        elment = filehindle.readline()
        if not elment:  # if file end ,raise stop
            if len(data) != 0:
                yield scaffold, data, start, end
            print("reach ends")
            raise StopIteration()
        elment = elment.strip()
        tmp = elment.split('\t')

        if tmp[0] != scaffold:  # if scaffold is not the same
            if len(data) != 0:
                yield scaffold, data, start, end
            data = [[int(tmp[1]), int(tmp[2])]]
            scaffold = tmp[0]
        elif int(tmp[1]) > end:
            if len(data) != 0:
                yield scaffold, data, start, end
            start = end
            end += window_size
            data = [[int(tmp[1]), int(tmp[2])]]
        else:
            data.append([int(tmp[1]), int(tmp[2])])


def setcolorpattern(scaffld, pos):
    '''
    define the color of region,intergeneic is blue 'b', intron is yellow 'y' exon is red 'r'

    :param scaffld: str
    :param pos: int
    :return: str
    '''
    lens = len(GTF_decoding.genomeDict[scaffld])
    idx = GTF_decoding.B_Search(lens, pos, scaffld)
    while True:
        if GTF_decoding.genomeDict[scaffld][idx].start > pos:
            return 'b'
        elif GTF_decoding.genomeDict[scaffld][idx].start < pos and GTF_decoding.genomeDict[scaffld][idx].end > pos:
            if GTF_decoding.genomeDict[scaffld][idx].IsinCommonIntrons(pos):
                return 'y'
            else:
                return 'r'
        elif GTF_decoding.genomeDict[scaffld][idx].end < pos and GTF_decoding.genomeDict[scaffld][idx + 1].start > pos:
            return 'b'
        else:
            idx += 1

def setxtick(start,end):
    #设置x标签,使得两者范围落在[1,10)之间
    counter=0
    divisor=0
    rangei=end-start
    if rangei<=3:
        return None
    while not (1<=rangei<10):
        rangei=rangei/10
        counter+=1
    if counter<3:
        unit=''
        divisor=10**0
    elif 3<=counter<6:
        unit='K'
        divisor=10**3
    elif 6<=counter<9:
        unit='M'
        divisor=10**6
    elif 9<=counter<12:
        unit='G'
        divisor=10**9
    elif 12<=counter<15:
        unit='T'
        divisor=10**12
    else:
        unit='P'
        divisor=10**15
    return divisor,unit,

def plotfig(windowdata, scaffold, start, end,pdfhindle):
    data = numpy.array(windowdata)
    fig = plt.bar(data[:, 0], data[:, 1], linewidth=0, width=1, color='r')
    plt.title("scaffold "+scaffold)
    for x in range(data.shape[0]):
        color = setcolorpattern(scaffold, data[x][0])
        fig[x].set_color(color)
    xticks=setxtick(start,end)
    plt.ylim(0,20)
    plt.ylabel('Depth')
    plt.xlabel('Position from '+start+' to '+end)
    m=range(start,end+(end-start)//10,(end-start)//10)
    q=map(lambda x:str(x/xticks[0])+xticks[1],m)
    plt.xticks(m,q)
    plt.show()
    # pdfhindle.savefig(fig)


def mian():
    pdf=PdfPages("results.pdf")
    options = _parse_args()
    GTF_decoding.decodegff(options.gtf)
    with open(options.input) as filehindle:
        itrator = slidingwindow(filehindle, options.window_size)
        count = random.randint(20,100)
        print(count)
        for scaffold, windowdata, start, end in itrator:
            if scaffold in GTF_decoding.genomeDict:
                # plotfig(windowdata, scaffold, start, end)
                count -= 1
                if count == 0:
                    plotfig(windowdata, scaffold, start, end,pdf)
                    break


if __name__ == '__main__':
    mian()
