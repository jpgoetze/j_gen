#!/usr/bin/env python2
#encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#This will compute the spectral overlap between an absorption spectrum given in a file (absorption, acceptor) and a fluorescence spectrum from a different file (donor). Both files must by in X Y format, X being the wavelength in nm and Y the molar extinction coeff in cm-1/M (for the acceptor) or an arbitrary unit (for the donor). One value pair per line. Take care to only include the region of your spectrum containing your states of interest.
#needs inputfile_acc, inputfile_donor. Output will be printed to prompt.

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$" #during a rain storm

def read_xyz_file(input):
	content=[]
	with open(input) as inp:
		for line in inp:
			vec=[float(field) for field in line.split()]
			content.append(vec)
	return content

def fast_eq(sp):
        llimit=350.00
        ulimit=550.00
        eq_sp=[]
        curr_nm=float(llimit)
        #curr_sp=float(sp[0][0])
        count=0
        while curr_nm<=ulimit:
            eq_vec=[]
            eq_vec.append(curr_nm)
            if curr_nm<float(sp[0][0]) or curr_nm>float(sp[len(sp)-1][0]):
                eq_vec.append(0.0)
            else:
                if float(sp[count][0])==curr_nm:
                    eq_vec.append(sp[count][1])
                else:
                    while not (float(sp[count][0])<curr_nm and float(sp[count+1][0])>curr_nm):
                        count+=1
                        if count+1==len(sp):
                            if float(sp[count][0])==curr_nm:
                                eq_vec.append(sp[count][1])
                                break
                            print("Grid equal error for "+str(curr_nm)+".")
                            exit(1)
                        elif float(sp[count][0])==curr_nm:
                            eq_vec.append(sp[count][1])
                            break
                    if float(sp[count][0])!=curr_nm:
                        dx=float(sp[count+1][0])-float(sp[count][0])
                        dy=float(sp[count+1][1])-float(sp[count][1])
                        dy/=dx
                        new_y=float(sp[count][1])+dy*(curr_nm-float(sp[count][0]))
                        eq_vec.append(new_y)
            eq_sp.append(eq_vec)
            curr_nm+=1.0
        return eq_sp

def compute_spec_comp(sysargs):
	args=[]
	for element in sysargs:
		args.append(str(element))
	del args[0]
	numbers=[]
	for i in range(0,len(args),2):
		numbers.append(int(args[i]))
	spectrafiles=[]
	for i in range(1,len(args),2):
		spectrafiles.append(str(args[i]))
	spectra=[]
	int_read=read_xyz_file("photon_int_terrestrial.txt")
	#int_read=[]
	for i in range(350,551):
		#int_read.append([int(i),1.0])
		spectra.append([int(i)])
	int_eq=fast_eq(int_read)
	intensities=[]
	for element in int_eq:
		intensities.append([element[1]])
	for element in spectrafiles:
		data=read_xyz_file(element)
		spectrum=fast_eq(data)
		for i in range(0,201):
			spectra[i].append(float(spectrum[i][1]))
	extfrac=[]
	import math
	for i in range(0,len(spectra)):
		evec=[]
		evec.append(spectra[i][0])
		curr_ext_sum=0.0
		for j in range(1,len(spectra[i])):
			curr_ext_sum+=float(spectra[i][j])
		for j in range(1,len(spectra[i])):
			if float(curr_ext_sum)==0.0:
				evec.append(0.0)
				continue
			evec.append(float(spectra[i][j])/float(curr_ext_sum))
		extfrac.append(evec)
	integral_vector=[]
	for j in range(1,len(extfrac[i])):
		overlap=0.0
		last_vals=[-1.0,-1.0]
		for i in range(0,len(extfrac)):	
			if last_vals[0]>0.0:
				left_prod=float(last_vals[1])
				right_prod=float(extfrac[i][j])
				avg_prod=(left_prod+right_prod)/2.0
				width=float(abs(extfrac[i][0])-float(last_vals[0]))
				overlap+=float(avg_prod)*width*float(intensities[i][0])
			last_vals[0]=float(extfrac[i][0])
			last_vals[1]=float(extfrac[i][j])
		integral_vector.append(overlap)
	overlap=0.0
	last_vals=[-1.0,-1.0]
	for i in range(0,len(intensities)):
		if last_vals[0]>0.0:
			left_prod=float(last_vals[1])
			right_prod=float(intensities[i][0])
			avg_prod=(left_prod+right_prod)/2.0
			width=float(abs(extfrac[i][0])-float(last_vals[0]))
			overlap+=float(avg_prod)*width*float(intensities[i][0])
		last_vals[0]=float(extfrac[i][0])
		last_vals[1]=float(intensities[i][0])
	for element in integral_vector:
		print('{:8.5E}\n'.format(float(element)/float(overlap)))

if __name__ == '__main__':
	import sys
	if len(sys.argv)<3:
		print("Must have at least 2 arguments. Exiting\n")
		exit(1)
	compute_spec_comp(sys.argv)
