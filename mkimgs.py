#!/usr/bin/python
# syntax: ./mkimgs.py dir obsid elowlist ehighlist

import sys
import os
import string

try:
    fobs = open(sys.argv[1]+'/'+sys.argv[2]+'/'+sys.argv[2]+'.dat', 'r')
except IOError:
    obsids=[sys.argv[2]]
    dir=sys.argv[1]+"/"+obsids[0]+"/event_cl/"
else:
    obsids = fobs.readlines()
    for i in range(len(obsids)):
        obsids[i]=obsids[i].rstrip()
    fobs.close()
    dir=sys.argv[1]+'/'+sys.argv[2]+'/'

slow=sys.argv[3].split(',')
shigh=sys.argv[4].split(',')
elow=[float(i) for i in slow]
ehigh=[float(i) for i in shigh]
clow=[str(int((i-1.6)/0.04+1)) for i in elow]
chigh=[str(int((i-1.6)/0.04)) for i in ehigh]

for det in ["A","B"]:
    xsel=open(dir+"xsel.xco","w")
    xsel.write("session1\n")
    iobs=0
    for obsid in obsids:
        edir=sys.argv[1]+'/'+obsid+'/'
        xsel.write("read event "+edir+"/event_cl/nu"+obsid+det+"01_cl.evt\n")
        if iobs == 0:
            xsel.write("./\n")
            xsel.write("yes\n")
        iobs=iobs+1

    for i in range(len(slow)):
        try:
            fobs=open(dir+"im"+det+slow[i]+"to"+shigh[i]+"keV.fits",'r')
        except IOError:
            blah=1
        else:
            os.system("rm -f -r "+dir+"im"+det+slow[i]+"to"+ \
                  shigh[i]+"keV.fits")
        xsel.write('filter pha_cutoff '+clow[i]+' '+chigh[i]+'\n')
        xsel.write("extract image\n")
        xsel.write("save image\n")
        xsel.write(dir+"im"+det+slow[i]+"to"+shigh[i]+"keV.fits\n")
        xsel.write("clear\n")
        xsel.write("pha_cutoff\n")
    xsel.write('set xyname det1x det1y\n')

    for i in range(len(slow)):
        try:
            fobs=open(dir+"im"+det+slow[i]+"to"+shigh[i]+"keVdet.fits",'r')
        except IOError:
            blah=1
        else:
            os.system("rm -f -r "+dir+"im"+det+slow[i]+"to"+ \
                  shigh[i]+"keVdet.fits")
        xsel.write('filter pha_cutoff '+clow[i]+' '+chigh[i]+'\n')
        xsel.write('filter column "RAWX=1:30 RAWY=1:30"\n')
        xsel.write("extract image\n")
        xsel.write("save image\n")
        xsel.write(dir+"im"+det+slow[i]+"to"+shigh[i]+"keVdet.fits\n")
    xsel.write("exit\n")
    xsel.write("no\n")
    xsel.close()
    os.system("xselect @"+dir+"xsel.xco")
    os.system("rm -r -f "+dir+"xsel.xco")


