#codig=UTF-8
import seaborn as sb
from matplotlib import pylab as pl
from matplotlib import pyplot as plt
import numpy
import optparse
import sys
import time
import GTF_decoding

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
    parser.add_option('-o', '--output', dest='output', type='string', help='input variation information file')
    parser.add_option('-d','--depth',dest='depth',type='int',default=5,help='depth')
    options, args = parser.parse_args()
    # positional arguments are ignored
    return options









def mian():
    options=_parse_args()



if __name__ == '__main__':
    mian()