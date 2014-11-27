#saxs_multiavg.py---------------------------------------------------------------------------
#
#   SAXS multiple file targeted average program
#   find best linear combination of N SAXS dats sets to Target data set
#	if buffer file(s) present (num_buf>0) will subtract buffer first
#	data files assumed to have a title line unless
#		 -notitle is used or name ends in .out which is assumed to be Gnom p(r) file
#	if output_file2 is specified, it is simple x y 2 column space separated  (useful for Gnom)
#
# to run:
#
#  python saxs_multiavg.py command file  output_file [output_file2] [-x=min:max] [-merge] [-notitle] [fit=min:max]
#
# switches:
#    -x=min:max     will limit x_axis info range example -x=.05:.25  will select q values between .05 and .25
#    -fit=xmin:xmax   will only using selected range for fitting
#    -merge	      will include target(s) in the final averages, otherwise targets used for scaling only
#    -notitle	      must use if input file does not have a title line
#
# command file:
#	num_buffer_files num_data_files 
#	directory
#	[buf target file (or avg buffer if num_buf=1)]
#	[buffer files1-N]
#	data target file
#	data files 1-N
#
#  NOTE: the target files are included in the file counts (num_buffer_files & num_data_files)
#
#   v1.00       7/21/06     DAA
#   v1.10       12/28/12    DAA  changed command file, added merge, -notitle, Gnom file handling
#   v1.20	1/11/13	    DAA  added 2nd output file
#   v1.25       1/18/13     DAA make -x= work and add -fit= switches
#
#..................................................................................
import sys
import numpy as np

def parse_line(files,switches):
    nargs = len(sys.argv)
#    print nargs,sys.argv
    nfiles = 0
    nswitches = 0
    for j in range(1,nargs):
#      print 'argv',sys.argv[j]
      if sys.argv[j].find('-') < 0:
          if nfiles== 0:
              files[0] = sys.argv[j]
          else:
            files.append(sys.argv[j])
#          print 'file:', j,nfiles,files[nfiles]
          nfiles+=1
      else:
          sw = sys.argv[j]
          neq = sw.find('=')
#          print 'neq,sw',neq,sw
          switch = ['',0,0.,0.,0.]
          if neq > 0:
            switch[0] = sw[:neq]
            nums = sw[neq+1:]
            n = 0
            nc = nums.find(':')
            while nc >= 0:
              switch[n+2] = float(nums[:nc])
              n+=1
              switch[1] = int(n)
              nums = nums[nc+1:]
              nc = nums.find(':')
            switch[1] = int(n+1)
            switch[n+2] = float(nums)
          else:
              switch[0] = sw
          if nswitches == 0:            
              switches[0] = switch
          else:
              switches.append(switch)
#          print 'switch:', j,nswitches,switches[nswitches]
          nswitches+=1
    return nfiles,nswitches
#end parseline


#interpolate data onto same q spacing as target
def interp(xy_t,xy_s,xmin_fit,xmax_fit,ind_o):
    done = False
    x1 = xy_s[0][0]
    start = 0
    low_fit = -1
    len_s = len(xy_s)-1
#    print 'lent,s',len(xy_t),len_s
    while xy_t[start][0] < x1:
        start+=1
#    print 'starting',x1,start,xy_t[start][0]
    indx = 0
    for j in range(start,len(xy_t)):
		x2 = xy_t[j][0]
		if x2<=xmax_fit: hi_fit=j
		if (low_fit<0) & (x2>=xmin_fit): low_fit = j
		while xy_s[indx][0] <= x2 and indx < len_s:
			indx+=1
#			print 'indx',indx,len_s,x2,xy_s[indx][0]
	    	if indx==len_s: 
				done = True
				break
		if done: break
		dx = xy_s[indx][0] - xy_s[indx-1][0]
		dy = xy_s[indx][1] - xy_s[indx-1][1]
#		if j<20: print 'dx,dy',dx,dy,len(xy_t[j]),ind_o
		if ind_o >= len(xy_t[j]):
	  		xy_t[j].append(xy_s[indx-1][1] + dy*(x2-xy_s[indx-1][0])/dx)
		else:
	  		xy_t[j][ind_o] =  xy_s[indx-1][1] + dy*(x2-xy_s[indx-1][0])/dx
#		if j<20: print 'after',xy_t[j]
    lims=[start,j,low_fit,hi_fit]
    return lims
#end interp()

#
# scale two curves together (ind_2 is target)
def buf_scale(data,low,hi,low_fit,hi_fit,ind_1,ind_2):
	sum_1 = 0.0
	sum_2 = 0.0
	sum_11 = 0.0
	sum_12 = 0.0
	rf1 = 0
	rf2 = 0
	for j in range(low_fit,hi_fit):
		sum_1 += data[j][ind_1]
		sum_2 += data[j][ind_2]
		sum_11 += data[j][ind_1]**2
		sum_12 += data[j][ind_1]*data[j][ind_2]
		rf1 += abs(data[j][ind_2] - data[j][ind_1])
	num = float(hi - low)
	d = num*sum_11 - sum_1**2
	b = (sum_11*sum_2 - sum_1*sum_12)/d
	a = (num*sum_12 - sum_1*sum_2)/d

	for j in range(low_fit,hi_fit):
		data[j][ind_1] = data[j][ind_1]*a + b
		rf2 += abs(data[j][ind_2] - data[j][ind_1])
	return a,b,rf1/sum_2,rf2/sum_2
#end buf_scale()


# scale two data-buffer curves together
def data_scale(data,low,hi,low_fit,hi_fit,ind_dat,ind_buf,ind_targ):
#	if (ind_buf<0):
#		a,b,rf1,rf2 =buf_scale(data,low,hi,low_fit,hi_fit,ind_dat,ind_targ)
#		return a,0,b,rf1,rf2
	sum_y = 0.0
	sum_d = 0.0
	sum_b = 0.0
	sum_d2 = 0.0
	sum_b2 = 0.0
	sum_yd = 0.0
	sum_yb = 0.0
	sum_db = 0.0
	rf1 = 0.0
	rf2 = 0.0
	a = np.zeros((3,3))
	vec = np.zeros((3,1))
	for j in range(low_fit,hi_fit):
		buf = 0.0
#		if j<10:print j,data[j][ind_targ],data[j][ind_dat]
		if ind_buf>0: buf = data[j][ind_buf]
		sum_y += data[j][ind_targ]
		sum_d += data[j][ind_dat]
		sum_b += buf
		sum_d2 += data[j][ind_dat]**2
		sum_b2 += buf**2
		sum_yd += data[j][ind_targ]*data[j][ind_dat]
		sum_yb += data[j][ind_targ]*buf
		sum_db += data[j][ind_dat]*buf
		rf1 += abs(data[j][ind_targ] - data[j][ind_dat])
	num = float(hi - low)
#	print 'ab',sum_y,sum_d,sum_b,sum_d2
#	print 'bc',sum_b2,sum_yd,sum_yb,sum_db
#	a = np.array(((sum_d2,sum_db,sum_d),(sum_db,sum_b2,sum_b),(sum_d,sum_b,num)))
	a[0][0] = sum_d2
	a[0][1] = sum_db
	a[0][2] = sum_d
	a[1][0] = sum_db
	a[1][1] = sum_b2
	a[1][2] = sum_b
	a[2][0] = sum_d
	a[2][1] = sum_b
	a[2][2] = num
	vec[0][0] = sum_yd
	vec[1][0] = sum_yb
	vec[2][0] = sum_y
	a_inv = np.linalg.pinv(a,1.e-6)
#	print 'matrix a',a
#	print 'inv',a_inv
#	print 'mul',np.matrixmultiply(a_inv,a)
#	vec = np.reshape(np.array(((sum_yd,sum_yb,sum_y))),(3,1))
#	print 'vec',vec
	abc = np.dot(a_inv,vec)
#	print 'abc',abc

	for j in range(low_fit,hi_fit):
		buf = 0.0
		if ind_buf>0: buf = data[j][ind_buf]
		data[j][ind_dat] = data[j][ind_dat]*abc[0] + buf*abc[1] + abc[2]
		rf2 += abs(data[j][ind_targ] - data[j][ind_dat])
	return float(abc[0]),float(abc[1]),float(abc[2]),float(rf1/sum_y),float(rf2/sum_y)
#end data_scale()

#avg set of curves, append on end of list
def data_avg(data,low,hi,ind_1,ncurves,ind_out,ind_targ,merge_target):

#	print 'avg: ',ind_1,ncurves,ind_out,ind_targ,merge_target,len(data[low])
	if ind_out >= len(data[low]):
		appnd = True
	else:
		appnd = False
	scl_t = 0.0
	if merge_target: scl_t = 1.0 
	scl = 1./(ncurves + scl_t)
#	print 'avg: ',scl,scl_t,len(data)
	for j in range(low,hi):
		avg = scl_t*data[j][ind_targ]
		for k in range(ind_1,ind_1+ncurves):
	  		avg += data[j][k]
		avg *= scl
		if appnd:
			data[j].append(avg)
		else:
			data[j][ind_out] = avg
	return
# end data_avg()


#afunc------------------------------------------------------------------
#
#  function being minimized by simplex maximizer
#
def afunc(var,data):
	ndata = len(data)
	error = 0
#	for j in range(7):	  print 'vars',j,var[j]
	for j in range(ndata):
	  ind = j + 7
	  o1 = data[j][3] - var[5]*data[j][9]
	  o2 = var[0]*(data[j][4] - var[6]*data[j][9])
	  o3 = var[1]*(data[j][5] - var[5]*data[j][9])
	  o4 = var[2]*(data[j][6] - var[6]*data[j][9])
	  o5 = var[3]*(data[j][7] - var[5]*data[j][9])
	  o6 = var[4]*(data[j][8] - var[6]*data[j][9])
	  obs = (o1+o2+o3+o4+o5+o6)/6.0
	  error += (obs - o1)**2
	  error += (obs - o2)**2
	  error += (obs - o3)**2
	  error += (obs - o4)**2
	  error += (obs - o5)**2
	  error += (obs - o6)**2
	return -error
# end(afunc)

#main here--------------------------------------------------------------------

print '\n\nSAXS multiple file averaging program v1.25\n'


#parse command line for files (should be 4) and optional switches
files=['']
switches=[['',0,0.,0.,0.]]
nfiles,nswitches = parse_line(files,switches)
#print '#files,switches',nfiles,nswitches
#print files

# look for switches
xmin = 0.0
xmax = 1.e10
xmin_fit = 0.0
xmax_fit = 1.e10
merge_target = False
low = 0
hi = 100000
low_fit = 0
hi_fit = 100000
notitle = False
if nswitches > 0:
     # print switches
     for j in range(nswitches):
       if switches[j][0] == '-x':
           if switches[j][1] >= 1: xmin = switches[j][2]
           if switches[j][1] >= 2: xmax = switches[j][3]
       if switches[j][0] == '-fit':
           if switches[j][1] >= 1: xmin_fit = switches[j][2]
           if switches[j][1] >= 2: xmax_fit = switches[j][3]
       if switches[j][0].find('-merg')>=0: merge_target = True
       if switches[j][0].find('-notit')>=0: notitle = True
if nfiles<2:
   print '\n**** Insufficient number of files found:',nfiles,files
   sys.exit('************ERROR: incorrect number of files ***********')

print 'Command file:   ',files[0]
print 'Output file: ',files[1]
f_c = open(files[0],'r')
f_out = open(files[1],'w')


nl = 0
nf = 0
conc =[0,0,0]
ifile = ['']
for line in f_c:
	ls = line.split()
	if nl == 0:
		nfil_b = long(ls[0])
		nfil_d = long(ls[1])
		print '# buffer,data files',nfil_b,nfil_d
	elif nl == 1:
		dir = line
	else:
		nam = dir[:len(dir)-1]+'/'+ls[0]
		print 'file #',nf+1,':',nam
		if nf == 0:
		  ifile = [nam]
#		  print 'nf',nf,ifile
		else:
		  ifile.append(nam)
#		  print 'nf',nf,ifile
		nf += 1
	nl +=1

nread1= int(-1)
if notitle: nread1 = int(0)
#print 'nam',0,ifile[0]
f_in = open(ifile[0],'r')
#if Gnom .out file then skip useless info
if ifile[0].find('.out') >0:
	for line in f_in:
		if line.find('P(R)') > 0:
			print 'target data file title',line
			nread1 = 0
			break
for line in f_in:
	if len(line) > 1:
		if nread1 <0: print 'target data file title:',line
		ls = line.split()
		if nread1==0:
			if float(ls[0]) >= xmin: data = [[float(ls[0]),float(ls[1])]]
			else: nread1 -=1
		elif nread1 >0:
			if (float(ls[0])>=xmin)&(float(ls[0])<= xmax):data.append([float(ls[0]),float(ls[1])])
			else: nread1 -=1
		nread1 +=1
nread1 = len(data)	
f_in.close()
#for k in range(5): print 'dat0',data[k]

for j in range(1,nf):
#	print 'nam',j,ifile[j]
	f_in = open(ifile[j],'r')
	#if Gnom .out file then skip useless info
	if ifile[0].find('.out') >0:
		for line in f_in:
			if line.find('P(R)') > 0:
				print 'title file #',j,line
				nread2 = 0
				break
	nread2 = int(-1)
	if notitle: nread2 = int(0)
	for line in f_in:
		if len(line) > 1:
			if nread2 <0: print 'title file #',j,line
			ls = line.split()
			if nread2==0:
				if float(ls[0]) >= xmin: data2 = [[float(ls[0]),float(ls[1])]]
				else: nread2 -=1
			elif nread2>0:
				if (float(ls[0])>=xmin)&(float(ls[0])<= xmax):data2.append([float(ls[0]),float(ls[1])])
				else: nread2 -=1
			nread2 +=1
	nread2 = len(data2)	
	f_in.close()
	
	lim = interp(data,data2,xmin_fit,xmax_fit,j+1)
	low = max(low,lim[0])
	hi = min(hi,lim[1])
	low_fit = max(low_fit,lim[2])
	hi_fit = min(hi_fit,lim[3]) 
	print 'file,lim',j+1,lim,low_fit,hi_fit,low,hi
st = 20*['']
st2 = 2*['']
sbuf = ['Buf1','Buf2','Buf3','Buf4','Buf5','Buf6','Buf7','Buf8','Buf9','Buf10','Buf11','Buf12']
sdat = ['Data1','Data2','Data3','Data4','Data5','Data6','Data7','Data8','Data9','Data10','Data11','Data12']
st[0] = 'X'
print '\nX axis selection range: ',xmin,xmax
if merge_target: print 'Will merge target(s) into final output'
print ' '
#
# scale buffers together
#
ncols = 1
#print 'ncols b',ncols,len(data[low])
if nfil_b > 1:
	ind_buf_targ =  1
	ind_buf_st = ind_buf_targ + 1
	ind_buf_avg = ind_buf_st + nfil_b + nfil_d -1
	ncols +=2
	print 'ind buffer targ,st,avg ',ind_buf_targ,ind_buf_st,ind_buf_avg
	for j in range(nfil_b-1):
		print 'buffer: scl,off,Rf_in,Rf_out ',j+1,buf_scale(data,low,hi,low_fit,hi_fit,ind_buf_st+j,ind_buf_targ)
		st[j+ind_buf_st] = sbuf[j]
		ncols +=1
	data_avg(data,low,hi,ind_buf_st,nfil_b-1,ind_buf_avg,ind_buf_targ,merge_target)
	st[ind_buf_targ] = 'Buf_target'
	st[ind_buf_avg] = 'Buf_avg'
elif nfil_b==1:
	ind_buf_avg = ind_buf_targ
	ncols +=1
else:
	ind_buf_avg = -1
#print 'ncols b',ncols,len(data[low]),len(st)
#		
# scale data - buffer
if nfil_d > 0:
	
	if nfil_b > 0: 
		ind_dat_targ = ind_buf_st + nfil_b -1
		ind_dat_st = ind_dat_targ + 1
		ind_dat_avg = ind_buf_avg + 1
	else:
		ind_dat_targ = 1
		ind_dat_st = ind_dat_targ + 1
		ind_dat_avg = ind_dat_st + nfil_d - 1
	print '\nind data targ,st,avg ',ind_dat_targ,ind_dat_st,ind_dat_avg
	ncols +=2
	for j in range(nfil_d-1):
		print 'data: scl_d,scl_b,off,Rf_in,Rf_out ',j+1,data_scale(data,low,hi,low_fit,hi_fit,j+ind_dat_st,ind_buf_avg,ind_dat_targ)
		st[j+ind_dat_st] = sdat[j]
		ncols +=1
	data_avg(data,low,hi,ind_dat_st,nfil_d-1,ind_dat_avg,ind_dat_targ,merge_target)
	st[ind_dat_targ] = 'Data_target'
	st[ind_dat_avg] = 'Data_avg'		
 #    print '\nRfactors (a,b,merged for  q range):',round(r[3],7),round(r[4],7),
  #  round(r[5],7)
  #  print 'Rfactors (a,b,merged for all data):',round(r[0],7),round(r[1],7),round(r[2],7)
#print 'ncols d',ncols,len(data[low]),len(st)
sep = chr(9)
f_out.write(sep.join(st)+'\n')

for j in range(low,hi):
	for k in range(ncols):
		st[k] = str(round(data[j][k],5))
	f_out.write(sep.join(st)+'\n')
f_out.close()

sep = ' '
if (nfiles >2) & (nfil_d>0):
	print '\nWriting second output file'
	f_out2 = open(files[2],'w')
	for j in range(low,hi):
		st2[0] = str(round(data[j][0],5))
		st2[1] = str(round(data[j][ind_dat_avg],5))
		f_out2.write(sep.join(st2)+'\n')
	f_out2.close()
print 'Done'
