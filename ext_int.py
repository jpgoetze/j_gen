#!/usr/bin/env python2
#encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#This will this will compute the spectral integral ( int(\epsilon(\nu)d\nu) ) for a given spectrum. Spectrum in wavenumbers.
#needs inputfile. Output will be printed to prompt.

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$" #during a rain storm

def read_xyz_file(input):
	content=[]
	with open(input) as inp:
		for line in inp:
                        vec=[float(field) for field in line.split()]
                        content.append(vec)
	return content

def compute_extinction_integral(infile1):
	data=read_xyz_file(infile1)
	data.sort()
	overlap=0.0
	last_vals=[-1.0,-1.0]
	for element in data:
		if last_vals[0]>0.0:
			curr_lambda=(float(element[0])+float(last_vals[0]))/2.0
			left_prod=float(last_vals[1])
			right_prod=float(element[1])
			avg_prod=(left_prod+right_prod)/2.0
			width=float(abs(float(element[0])-float(last_vals[0])))
			overlap+=float(avg_prod)*width
		last_vals[0]=element[0]
		last_vals[1]=element[1]
	print('{:8.5E}'.format(float(overlap)))

if __name__ == '__main__':
        import sys
        if len(sys.argv)!=2:
                print("Must have 1 argument. Exiting\n")
                exit(1)
        compute_extinction_integral(sys.argv[1])
