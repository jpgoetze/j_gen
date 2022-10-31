#!/usr/bin/env python2
#encoding: ISO-8859-15

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#This will combine a photonic irradiation spectrum with a molar extinction spectrum to produce absorbed photons per wavelength. Requires a wavelength to normalize, with set this specific wavelength to absorbance of 1/T of the incoming light.
#needs input(photons s-1 m-2), input2(molar extinction spectrum), norm_wavelength, T, output file.

__author__ = "jangoetze"
__date__ = "$15-May-2018 17:02:17$" #during a rain storm

def read_xyz_file(input):
	content=[]
	with open(input) as inp:
		for line in inp:
                        vec=[float(field) for field in line.split()]
                        content.append(vec)
	return content

def create_equalized_grid(data,data2):
        full_set_x=[]
        data2_count=0
        for i in range(0,len(data)):
                if data2_count<len(data2):
                        while float(data2[data2_count][0])<float(data[i][0]):
                                full_set_x.append(float(data2[data2_count][0]))
                                data2_count+=1
                                if data2_count==len(data2):
                                        break
                full_set_x.append(float(data[i][0]))
                if data2_count<len(data2):
                        if float(data[i][0])==float(data2[data2_count][0]):
                                data2_count+=1
        while data2_count<len(data2):
                full_set_x.append(float(data2[data2_count][0]))
                data2_count+=1
        interval=-1.0
        for i in range(1,len(full_set_x)):
                curr_x=full_set_x[i]
                last_x=full_set_x[i-1]
                if interval<0.0:
                        interval=(full_set_x[i]-full_set_x[i-1])
                else:
                        new_interval=(full_set_x[i]-full_set_x[i-1])
                        if new_interval<interval:
                                interval=new_interval
        if interval<0.0:
                print("Interval < 0.0 Exiting")
                exit(1)
        curr_x=full_set_x[0]
        eq_data=[]
        while curr_x<full_set_x[len(full_set_x)-1]:
                curr_y=0.0
                curr_z=0.0
                for i in range(0,len(data)):
                        if float(data[i][0])==curr_x:
                                curr_y=float(data[i][1])
                                break
                        if float(data[i][0])>curr_x:
                                if i==0:
                                        break
                                else:
                                        prev_x=float(data[i-1][0])
                                        prev_y=float(data[i-1][1])
                                        next_x=float(data[i][0])
                                        next_y=float(data[i][1])
                                        delta_prev_x=curr_x-prev_x
                                        delta_x=next_x-prev_x
                                        delta_y=next_y-prev_y
                                        dydx=delta_y/delta_x
                                        curr_y=(dydx*delta_prev_x)+prev_y
                                        break
                for i in range(0,len(data2)):
                        if float(data2[i][0])==curr_x:
                                curr_z=float(data2[i][1])
                                break
                        if float(data2[i][0])>curr_x:
                                if i==0:
                                        break
                                else:
                                        prev_x=float(data2[i-1][0])
                                        prev_y=float(data2[i-1][1])
                                        next_x=float(data2[i][0])
                                        next_y=float(data2[i][1])
                                        delta_prev_x=curr_x-prev_x
                                        delta_x=next_x-prev_x
                                        delta_y=next_y-prev_y
                                        dydx=delta_y/delta_x
                                        curr_z=(dydx*delta_prev_x)+prev_y
                                        break
                eq_data_vec=[]
                eq_data_vec.append(curr_x)
                eq_data_vec.append(curr_y)
                eq_data_vec.append(curr_z)
                eq_data.append(eq_data_vec)
                curr_x+=interval
        return eq_data

def compute_overlap2(infile,infile2,norm,trans,outfile):
        import numpy as np
        T=0.0
        if trans=="e":
                T+=np.exp(1.0)
        else:
                T+=float(trans)
        data=read_xyz_file(infile)
        data2=read_xyz_file(infile2)
        eq_data=create_equalized_grid(data,data2)
        e_photons_wavelength=float(norm)
        norm_ext=0.0
        for element in eq_data:
                if float(element[0])==e_photons_wavelength:
                        norm_ext=element[2]
                        break
        if norm_ext==0.0:
                print("The normalization wavelength was not found in your spectrum that has the highest resolution in the overlapping region. Please restate the value with a properly existing number. Exiting.")
                exit(1)
        final_data=[]
        for element in eq_data:
                final_vec=[]
                final_vec.append(element[0])
                final_vec.append(float(element[1])*np.exp(-1.0*float(element[2])*np.log(1.0/T)/float(norm_ext)))
                final_data.append(final_vec)
        with open(outfile,"w") as out:
                for element in final_data:
                        out.write(str(float(element[0]))+" "+str(float(element[1]))+"\n")

if __name__ == '__main__':
        import sys
        if len(sys.argv)!=6:
            print("Must have 5 arguments. Exiting\n")
            exit(1)
        compute_overlap2(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])


