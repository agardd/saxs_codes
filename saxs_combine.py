#saxs_combine.py---------------------------------------------------------------------------
#
#   SAXS  multiple file merging  program
#   find best linear combination of N SAXS dats sets to Target data set
#	output file is merged sum of all data files (only includes target file if -merge is used) 
#	data files assumed to have a title line unless
#		 -notitle is used or name ends in .out which is assumed to be Gnom p(r) file
#
# to run:
#
#  python saxs_combine.py command_file  output_file [-x=xmin:xmax] [-merge] [-notitle] [-pad=] [-fit=]
#
#switches:
#	-x=xmin:xmax	sets min,max selection for input data
#	-fit=xmin:xmax	sets min,max selection for fitting data
#	-pad=xpad	will pad data to xmax for output purposes, fits only done over real data
#	-merge		will include target in final combined output
#	-notitle	needed if there is no title on input data
#
#command_file
#	directory
#	dat target_file
#	data files1-N
#
#   v1.00       7/4/07     DAA
#   v1.20	12/21/12   DAA   modified to handle Gnom.out files
#   v1.25	12/27/12   DAA   aded -merge -offset (not yet working) -notitle switches
#   v1.30	1/16/13	   DAA   added -pad -fit switches, bug fixes
#..................................................................................
import sys
from numpy import *
from numpy.linalg import inv
from numpy.linalg import pinv

def parse_line(files,switches):
    nargs = len(sys.argv)
    print nargs,sys.argv
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
#		print 'interp,j,start,len,indx,x2',j,start,len(xy_t),indx,x2
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
    print 'interp lims',lims
    return lims
#end interp()

# scale multiple data curves together and merge
def data_merge(data,low,hi,low_fit,hi_fit,nmerge,ind_st,ind_targ,ind_out,merge_target,calc_offset):

	use = nmerge*[1]
	offset = 0.0
	rf = 0.0
	sum = 0.0
	if calc_offset:
	  nuse,scl,offset = scales_offset(data,low,hi,low_fit,hi_fit,nmerge,use,ind_st,ind_targ)
	else:
	  nuse,scl = scales(data,low,hi,low_fit,hi_fit,nmerge,use,ind_st,ind_targ)
	if merge_target:
	  scl_t = 1./(nuse+1)
	  scl_c = nuse*scl_t
	else:
	  scl_c = 1.0
	  scl_t = 0.0
	for i in range(low_fit,hi_fit):
	    calc = offset
	    for k in range(nmerge):
	    	calc += data[i][k+ind_st]*scl[k]
	    rf += abs(data[i][ind_targ] - calc)
	    sum += data[i][ind_targ]
	rf /= sum
#	print 'scls',ok,scl,float(rf)
#
	for j in range(low,hi):
	  if ind_out >= len(data[j]):
	  	data[j].append(offset*scl_c + data[i][ind_targ]*scl_t)
	  else:
	      data[j][ind_out] = offset*scl_c + data[i][ind_targ]*scl_t
	  for k in range(nmerge):
		data[j][ind_out] += data[j][k+ind_st]*scl[k]*scl_c
	return float(rf)
#end data_merge()


#
#  calc best scales, if any scale is < 0 then set to zero and rescale remainder
def scales(data,low,hi,low_fit,hi_fit,nmerge,use,ind_st,ind_targ):
	inds = nmerge*[-1]
	scls = nmerge*[0.0]
	scli = nmerge*[0.0]
	off = nmerge*[0.0]
	for j in range (nmerge):
	  scli[j],off[j],rf1,rf2 = pair_scale(data,low,hi,low_fit,hi_fit,j+ind_st,ind_targ,False)
	  print 'individual scales: set,scl,offset: %d %.6f, %.6f, rf: %.4f' % (j+1,scli[j],off[j],rf2)
	for loop in range (nmerge):
		nuse = 0
		for j in range (nmerge):
	  		if use[j] >0:
	    			inds[nuse] = j
	    			nuse += 1
		a = zeros((nuse,nuse))
		vec = zeros((nuse,1))
		print 'lo,hi,nmerge,ind_st.ind_target,nuse ',low,hi,nmerge,ind_st,ind_targ,nuse 
		for i in range(low_fit,hi_fit):
	  		for j in range(nuse):
	      			vec[j][0] += data[i][inds[j]+ind_st]*data[i][ind_targ]
	      			for k in range(nuse):
	        			a[j][k] += data[i][inds[j]+ind_st]*data[i][inds[k]+ind_st]
#		print 'matrix a',a
		a_inv = pinv(a,1.e-6)
#		print 'pinv',a_inv
#		print 'mul',dot(a_inv,a)
#		print 'vec',vec
		scl = dot(a_inv,vec)
		ok = 1
		for j in range(nuse):
			if scl[j] < 0:	
				use[j] = 0
				ok = -1
			else: scls[inds[j]] = scl[j]
	  	if ok >0: break
	nuse = 0
	st = 0.0
	print 'nmerge',nmerge,len(scl)
	for j in range (len(scl)):
	  st += scl[j]
	  if use[j] >0: nuse += 1
	if nuse<nmerge: print 'sets used',use
	print ' '
	print 'Scale, scale_total, percentages:'
	for j in range (nmerge):
	  print '%d,  %.4f,  %.4f, %.4f'% (j+1,scls[j],scli[j]*scls[j],scls[j]/st)
	return nuse,scls
# end scales()

#
#  calc best scales, if any scale is < 0 then set to zero and rescale remainder
def scales_offset(data,low,hi,nmerge,use,ind_st,ind_targ):
	inds = nmerge*[-1]
	scls = nmerge*[0.0]
	scli = nmerge*[0.0]
	off = nmerge*[0.0]
	offset = 0.0
	for j in range (nmerge):
	  scli[j],off[j],rf1,rf2 = pair_scale(data,low,hi,low_fit,hi_fit,j+ind_st,ind_targ,True)
	  print 'individual scales: set,scl,offset: %d %.4f, %.4f, rf: %.4f' % (j+1,scli[j],off[j],rf2)
	for loop in range (nmerge):
		nuse = 0
		for j in range (nmerge):
	  		if use[j] >0:
	    			inds[nuse] = j
	    			nuse += 1
		a = zeros((nuse+1,nuse+1))
		vec = zeros((nuse+1,1))
		print 'lo,hi,nmerge,ind_st.ind_target,nuse ',low,hi,nmerge,ind_st,ind_targ,nuse 
		for i in range(low,hi):
			vec[nuse][0] += data[i][ind_targ]
	  		for j in range(nuse):
	      			vec[j][0] += data[i][inds[j]+ind_st]*data[i][ind_targ]
	      			a[nuse][j] += data[i][inds[j]]
	      			for k in range(nuse):
	        			a[j][k] += data[i][inds[j]+ind_st]*data[i][inds[k]+ind_st]
		a[nuse][nuse] = hi - low
#		print 'matrix a',a
		a_inv = pinv(a,1.e-6)
#		print 'pinv',a_inv
#		print 'mul',dot(a_inv,a)
#		print 'vec',vec
		scl = dot(a_inv,vec)
		ok = 1
		for j in range(nuse):
			if scl[j] < 0:	
				use[j] = 0
				ok = -1
			else: scls[inds[j]] = scl[j]
	  	if ok >0: break
	print 'offset,scale factors:',scl[nuse]
	print scls
	nuse = 0
	for j in range (nmerge):
	  if use[j] >0: nuse += 1
	if nuse<nmerge: print 'sets used',use
	return nuse,scls,scl[nuse]
# end scales_offset()

# scale two curves together (ind_2 is target)
def pair_scale(data,low,hi,low_fit,hi_fit,ind_1,ind_2,calc_offset):
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
	if calc_offset:
		d = num*sum_11 - sum_1**2
		b = (sum_11*sum_2 - sum_1*sum_12)/d
		a = (num*sum_12 - sum_1*sum_2)/d
	else:
		a = sum_12/sum_11
		b = 0.0
	for j in range(low,hi):
		data[j][ind_1] = data[j][ind_1]*a + b
		rf2 += abs(data[j][ind_2] - data[j][ind_1])
	return a,b,rf1/sum_2,rf2/sum_2
#end pair_scale()


#main here--------------------------------------------------------------------

print '\n\n SAXS linear combination program v1.30\n'


#parse command line for files (should be 4) and optional switches
files=['']
switches=[['',0,0.,0.,0.]]
nfiles,nswitches = parse_line(files,switches)
print '#files,switches',nfiles,nswitches
print files

# look for switches
xmin = 0.0
xmax = 1.e10
xmin_fit = 0.0
xmax_fit = 1.e10
merge_target = False
calc_offset = False
notitle = False
pad = False
pad_hi = 0
low = 0
low_fit = 0
hi = 100000
hi_fit = hi
nf = 2
if nswitches > 0:
     print switches
     for j in range(nswitches):
       if switches[j][0] == '-x':
           if switches[j][1] >= 1: xmin = switches[j][2]
           if switches[j][1] >= 2: xmax = switches[j][3]
       if switches[j][0] == '-fit':
           if switches[j][1] >= 1: xmin_fit = switches[j][2]
           if switches[j][1] >= 2: xmax_fit = switches[j][3]
       if switches[j][0].find('-mer')>=0: merge_target = True
###       if switches[j][0].find('-off')>=0: calc_offset = True
       if switches[j][0].find('-notit')>=0: notitle = True
       if switches[j][0].find('-pad')>=0: 
           pad = True
           if switches[j][1] >= 1: pad_hi = switches[j][2]
if nfiles!=nf:
   print '\n**** Incorrect number of files, should be:',nf, 'found:',nfiles,files
   sys.exit('************ERROR: incorrect number of files ***********')

print '\nCommand file:   ',files[0]
print 'Output file: ',files[1]
f_c = open(files[0],'r')
f_out = open(files[1],'w')
nl = 0
nfils = 0
ifile = ['']
for line in f_c:
	ls = line.split()
	if nl == 0:
		dir = line
	else:
		nam = dir[:len(dir)-1]+'/'+ls[0]
		print 'nf,nam',nfils+1,nam
		if nfils == 0:
		  ifile = [nam]
		else:
		  ifile.append(nam)
		nfils += 1
	nl +=1

nread1= int(-1)
if notitle: nread1 = int(0)
print 'target file',ifile[0]
f_in = open(ifile[0],'r')
#if Gnom .out file then skip useless info
if ifile[0].find('.out') >0:
	for line in f_in:
		if line.find('P(R)') > 0:
			print 'target title',line
			nread1 = 0
			break
# read in data target file
for line in f_in:
	if len(line) > 1:
		if nread1 <0: print 'target title:',line
		ls = line.split()
		if nread1==0:
			if float(ls[0]) >= xmin: data = [[float(ls[0]),float(ls[1])]]
			else: nread1 -=1
		elif nread1 >0:
			if line.find('Reciprocal') > 0: break
			if (float(ls[0])>=xmin)&(float(ls[0])<= xmax):data.append([float(ls[0]),float(ls[1])])
			else: nread1 -=1
		nread1 +=1
nread1 = len(data)	
f_in.close()
xmax_fit = min(xmax_fit,max(pad_hi,data[nread1-1][0])) #****
if pad:
	x = data[nread1-1][0]
	dx = x - data[nread1-2][0]
	while (x + dx) <= pad_hi:
	  	x += dx
		data.append([float(x),float(0.0)])
		nread1 +=1
print 'target xmax for fit, last: ',xmax_fit,data[nread1-1][0]	
#now read in data files
for j in range(1,nfils):
	print '\nnam',j+1,ifile[j]
	f_in = open(ifile[j],'r')
	nread2 = int(-1)
	if notitle: nread2 = int(0)
	#if Gnom .out file then skip useless info
	if ifile[j].find('.out') >0:
		for line in f_in:
			if line.find('P(R)') > 0:
				print 'data title',line
				nread2 = 0
				break
	for line in f_in:
		if len(line) > 1:
			if nread2 <0: print 'data title:',line
			ls = line.split()
			if nread2==0:
				if float(ls[0]) >= xmin: data2 = [[float(ls[0]),float(ls[1])]]
				else: nread2 -=1
			elif nread2 >0:
				if line.find('Reciprocal') > 0: break
				if (float(ls[0])>=xmin)&(float(ls[0])<= xmax):data2.append([float(ls[0]),float(ls[1])])
				else: nread2 -=1
			nread2 +=1
	nread2 = len(data2)	
	f_in.close()
	xmax_fit = min(xmax_fit,max(pad_hi,data2[nread2-1][0]))  #****
	
	if pad:
		x = data2[nread2-1][0]
		while (x + dx) <= pad_hi:
	  		x += dx
			data2.append([float(x),float(0.0)])
			nread2 +=1
	print 'file, xmax for fit, last: ',j+1,xmax_fit,data2[nread2-1][0]
	lim = interp(data,data2,xmin_fit,xmax_fit,j+1)
	low = max(low,lim[0])
	hi = min(hi,lim[1]) #hi for output
	low_fit = max(low_fit,lim[2])
	hi_fit = min(hi_fit,lim[3]) 
	print 'file,lims',j+1,lim,low_fit,hi_fit,low,hi
	

print' '
print 'X axis selection range: ',xmin,xmax
print 'X axis fit range: ',xmin_fit,xmax_fit,low_fit,hi_fit
if pad: print 'data padded to ',pad_hi
if calc_offset: print 'Offset being calculated'
if merge_target: print 'Will merge target into final output'
print ' '
nmerge = nfils - 1
st = 20*['']
sdat = ['Data1','Data2','Data3','Data4','Data5','Data6','Data7','Data8','Data9','Data10','Data11','Data12']
st[0] = 'X'
for j in range(nmerge):
	st[j+ 1] = sdat[j]
st[nmerge+1] = 'Data_target'
st[nmerge+2] = 'Data_merge'
#		
# merge data
#
ind_st = 2
ind_targ = 1
ind_out = ind_st + nmerge
print 'inds',ind_st,ind_targ,ind_out
r = data_merge(data,low,hi,low_fit,hi_fit,nmerge,ind_st,ind_targ,ind_out,merge_target,calc_offset)
print 'Merge Rfactor: %.4f'% r
print '\nDone\n'
#
sep = chr(9)
f_out.write(sep.join(st)+'\n')
for j in range(low,hi):
	st[0] = str(round(data[j][0],5))
	for k in range(nmerge):
		st[k+1] = str(round(data[j][k+2],5))
	st[nmerge+1] = str(round(data[j][1],5))
	st[nmerge+2] = str(round(data[j][ind_out],5))
	f_out.write(sep.join(st)+'\n')
f_out.close()
