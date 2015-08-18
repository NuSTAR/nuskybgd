#!/usr/bin/python
# syntax: ./getspecnoarf.py dir obsid regcore outdir det

import sys
import os
import string

dir=sys.argv[1]
obsid=sys.argv[2]
regcore=sys.argv[3]
outdir=sys.argv[4]
det=sys.argv[5]
dir=dir+"/"+obsid+"/"
rmfarf=" runmkarf=no runmkrmf=yes "

os.system("nuproducts "+ \
          "infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
          " srcregionfile="+dir+"event_cl/"+regcore+".reg "+ \
          " indir="+dir+"event_cl "+ \
          " outdir="+dir+"event_cl/"+outdir+rmfarf+ \
          " bkgextract=no "+ \
          " lcfile=NONE "+ \
          " instrument=FPM"+det+" steminputs=nu"+obsid+ \
          " stemout="+regcore+" boxsize=20 "+ \
          " clobber=yes")

os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.rmf "+ \
          dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig.rmf")
os.system("cmprmf "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig.rmf "+ \
          dir+"event_cl/"+outdir+"/"+regcore+"_sr.rmf 1e-6")
os.system("rm -f "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_g30.pha")
os.system("grppha "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.pha "+ \
          dir+"event_cl/"+outdir+"/"+regcore+ \
         "_sr_g30.pha 'chkey RESPFILE "+dir+"event_cl/"+outdir+"/"+ \
          regcore+"_sr.rmf & group min 30 & exit'")
