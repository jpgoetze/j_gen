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

def compute_overlap_full(infile1,infile2):
        data=read_xyz_file(infile1)
        data2=read_xyz_file(infile2)
        import imp
        eq_grid = imp.load_source("operations", str("../blue_light/integrate_photons.py"))
        eq_data=eq_grid.create_equalized_grid(data,data2)
        overlap=0.0
        flu_norm=0.0
        last_vals=[-1.0,-1.0,-1.0]
        for element in eq_data:
                if last_vals[0]>0.0:
                        curr_lambda=(float(element[0])+float(last_vals[0]))/2.0
                        left_prod=float(last_vals[1])*float(last_vals[2])
                        right_prod=float(element[1])*float(element[2])
                        avg_prod=(left_prod+right_prod)/2.0
                        width=float(abs(float(element[0])-float(last_vals[0])))
                        overlap+=float(avg_prod)*float(pow(float(curr_lambda),4.0))*width
                        avg_flu=(float(last_vals[2])+float(element[2]))/2.0
                        flu_norm+=avg_flu*width
                last_vals[0]=element[0]
                last_vals[1]=element[1]
                last_vals[2]=element[2]
        print('{:8.5E}'.format(float(overlap/flu_norm)))
        with open("temp","a") as out:
                out.write(('{:.5E}'.format(float(overlap/flu_norm))))
                out.write("\n")

if __name__ == '__main__':
        import sys
        if len(sys.argv)!=3:
                print("Must have 2 arguments. Exiting\n")
                exit(1)
        compute_overlap_full(sys.argv[1],sys.argv[2])
