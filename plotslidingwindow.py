# codig=UTF-8
import seaborn as sb
from matplotlib import pylab
from matplotlib import pyplot as plt
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


def plotfig(windowdata, scaffold, start, end):
    data = numpy.array(windowdata)
    fig = plt.bar(data[:, 0], data[:, 1], linewidth=0, width=1, color='r')
    for x in range(data.shape[0]):
        color = setcolorpattern(scaffold, data[x][0])
        fig[x].set_color(color)
    plt.show()


def mian():
    options = _parse_args()
    GTF_decoding.decodegff(options.gtf)
    with open(options.input) as filehindle:
        itrator = slidingwindow(filehindle, options.window_size)
        count = random.randint(50, 100)
        print(count)
        for scaffold, windowdata, start, end in itrator:
            if scaffold in GTF_decoding.genomeDict:
                # plotfig(windowdata, scaffold, start, end)
                count -= 1
                if count == 0:
                    plotfig(windowdata, scaffold, start, end)
                    break


if __name__ == '__main__':
    mian()
