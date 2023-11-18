#!/usr/bin/python

import sys
import gzip

cutoff = float(sys.argv[1])
nperbatch = int(sys.argv[2])
cm = sys.argv[3]
runs = sys.argv[4]
print "Pairwise global IBD cut-off: %f" % cutoff
print "Number of lines processed in one batch: %d" % nperbatch
print "Length: %s cm" % cm
print "Replicate runs: %s" % runs

fds = []
buf = []
curline = []
for chr in range(1, 23): # 1-22
    path = cm + 'cm/' + runs + '/all_' + cm + '.max.global.chr' + str(chr) + '.gz'
    fd = gzip.open(path, 'rb')
    fds.append(fd)
    for i in range(nperbatch):
        line = fd.readline()
        line.strip()
        items = line.split('\t')
        buf.append([int(items[1]), int(items[2]), float(items[3])])
    curline.append([int(items[1]), int(items[2])])
curline_sorted = sorted(curline, key=lambda x: (x[0], x[1]))
buf.sort(key=lambda x: (x[0], x[1]))
outpath = cm + 'cm/' + runs + '/all_' + cm + '.max.global.cut' + str(cutoff) + '.gz'
outfile = gzip.open(outpath, 'wb')
data = []
while True:
    if buf[0][0] > curline_sorted[0][0]:
        break
    if buf[0][0] == curline_sorted[0][0] and buf[0][1] >= curline_sorted[0][1]:
        break
    tmpdata = buf.pop(0)
    if not data:
        data = tmpdata
    elif tmpdata[0] == data[0] and tmpdata[1] == data[1]:
        data[2] = data[2] + tmpdata[2]
    else:
        if data[0] == data[1] or data[2] > cutoff:
            outfile.write(str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\n')
        data = tmpdata
while curline_sorted:
    for chr in range(22):
        if curline[chr] != curline_sorted[0]:
            continue
        for i in range(nperbatch):
            line = fds[chr].readline()
            if not line:
                curline[chr] = [0, 0]
                break
            line.strip()
            items = line.split('\t')
            buf.append([int(items[1]), int(items[2]), float(items[3])])
        if curline[chr][0] != 0:
            curline[chr] = [int(items[1]), int(items[2])]
    curline_sorted = sorted(curline, key=lambda x: (x[0], x[1]))
    while curline_sorted and curline_sorted[0][0] == 0:
        curline_sorted.pop(0)
    buf.sort(key=lambda x: (x[0], x[1]))
    while True:
        if curline_sorted and buf[0][0] > curline_sorted[0][0]:
            break
        if curline_sorted and buf[0][0] == curline_sorted[0][0] and buf[0][1] >= curline_sorted[0][1]:
            break
        if not buf:
            if data[0] == data[1] or data[2] > cutoff:
                outfile.write(str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\n')
            data = []
            break
        tmpdata = buf.pop(0)
        if not data:
            data = tmpdata
        elif tmpdata[0] == data[0] and tmpdata[1] == data[1]:
            data[2] = data[2] + tmpdata[2]
        else:
            if data[0] == data[1] or data[2] > cutoff:
                outfile.write(str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\n')
            data = tmpdata
for fd in fds:
    fd.close()
outfile.close()


