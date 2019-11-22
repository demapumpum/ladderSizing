# Script for extracting features of the peaks in all fsa files (data extracted for Hsc40 loci only)
# To extract data for other loci, change line 57 to 'data = record.annotations['abif_raw'][*]'	* - loci of interest, see lines 20-24
# Features:
# 			area_of_peaks	: Area under the curve in the region
# 			number_of_peaks	: Number of peaks in the region
# 			length_in_bp	: Corresponding length (in bp) of the detected peak
# 			height_of_peaks	: Height of the detected peak

from Bio import SeqIO
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
import pandas as pd
import math
from ladder_fit import convert_to_bp, convert_to_index, find_lower, find_upper
from findpeaks import findpeaks as fp


file = pd.read_csv('Channel1_Data_Normalized.csv')	# 'Channel1_Data_Normalized.csv'
a = 'DATA1'		# Hsc 40
b = 'DATA2'		# Hsc 24
c = 'DATA3'		# empty
d = 'DATA4'		# Hsc 20
e = 'DATA105'
dye = [500, 490, 450, 400, 350, 340, 300, 250, 200, 160, 150, 139, 100, 75, 50, 35]

area_of_peaks = []
label = []
file_name = []
# number_of_peaks = []
length_in_bp = []
height_of_peaks = []
noise_count = 0

for i in range(len(file)):
	filename = file.iat[i, 0]
	alelle1 = file.iat[i, 1]
	alelle2 = file.iat[i, 2]

	if filename == 'ROM_12_41_Hos' or filename == 'ROM_12_42_Hos' or filename == 'SAM_12_34_Hos' or filename == 'SAM_12_36_Hos':
		continue

	if filename != 'POP':
		if filename[0] == 'A':
			filename = filename.replace('A_', '')
		if filename[4] == 'C':
			filename = filename.replace('_CON', '')
		elif filename[4] == 'T':
			filename = filename.replace('_TIG', '')
		abif_file = 'A_' + filename + '.fsa'
		record = SeqIO.read(abif_file, 'abi')
		print('\nopening fsa file ({}/{}): {}'.format(i,len(file),abif_file))


		index_min = find_lower(record.annotations['abif_raw'][e], dye)
		index_max = find_upper(record.annotations['abif_raw'][e], dye)
		data = record.annotations['abif_raw'][a]

		height = []
		index_of_peaks = []
		for x in range(index_min, index_max):
			# print("#",end='')
			converted_bp = convert_to_bp(x, record.annotations['abif_raw'][e], dye)
			in_range = converted_bp >= alelle1 - 1.5 and converted_bp <= alelle1 + 1.5
			if in_range:
				index_of_peaks.append(x)
				height.append(data[x])
			elif converted_bp > alelle1+1.5:
				break
		alelle1_index = index_of_peaks[height.index(max(height))]


		if alelle1 != alelle2:

			height = []
			index_of_peaks = []
			for x in range(index_min, index_max):
			# print("#",end='')
				converted_bp = convert_to_bp(x, record.annotations['abif_raw'][e], dye)
				in_range = converted_bp >= alelle2 - 1.5 and converted_bp <= alelle2 + 1.5
				if in_range:
					index_of_peaks.append(x)
					height.append(data[x])
				elif converted_bp > alelle2+1.5:
					break

			alelle2_index = index_of_peaks[height.index(max(height))]

			if alelle2_index > alelle1_index:
				if alelle1_index+80>=alelle2_index:
					two_peaks = True
				else:
					two_peaks = False
			else:
				if alelle1_index-80<=alelle2_index:
					two_peaks = True
				else:
					two_peaks = False

		all_peaks = fp.findpeaks(data, spacing=25, limit=25)

		for i in range(len(all_peaks)-1):
			lower = all_peaks[i] - 40
			upper = all_peaks[i] + 40
			if lower < index_min or upper > index_max:
				continue

			if noise_count >= 1268:
				if alelle1 != alelle2:
					if all_peaks[i] == alelle1_index or all_peaks[i] == alelle2_index:
						file_name.append(filename)
						area_of_peaks.append(np.trapz(data[lower:upper]))
						# number_of_peaks.append(len(fp.findpeaks(data[lower:upper], spacing=5,limit=15)))
						length_in_bp.append(round(convert_to_bp(all_peaks[i], record.annotations['abif_raw'][e], dye)))
						height_of_peaks.append(data[all_peaks[i]])
					else:
						continue
				elif alelle1 == alelle2:
					if all_peaks[i] == alelle1_index:
						file_name.append(filename)
						area_of_peaks.append(np.trapz(data[lower:upper]))
						# number_of_peaks.append(len(fp.findpeaks(data[lower:upper], spacing=5,limit=15)))
						length_in_bp.append(round(convert_to_bp(all_peaks[i], record.annotations['abif_raw'][e], dye)))
						height_of_peaks.append(data[all_peaks[i]])
					else:
						continue
			else:
				file_name.append(filename)
				area_of_peaks.append(np.trapz(data[lower:upper]))
				# number_of_peaks.append(len(fp.findpeaks(data[lower:upper], spacing=5,limit=15)))
				length_in_bp.append(round(convert_to_bp(all_peaks[i], record.annotations['abif_raw'][e], dye)))
				height_of_peaks.append(data[all_peaks[i]])


			if alelle1 != alelle2:
				if two_peaks == True:
					if all_peaks[i] == alelle1_index or all_peaks[i] == alelle2_index:
						label.append(1)
					else:
						label.append(0)
						noise_count += 1
				elif two_peaks == False:
					if all_peaks[i] == alelle1_index or all_peaks[i] == alelle2_index:
						label.append(1)
					else:
						label.append(0)
						noise_count += 1
			elif alelle1==alelle2:
				# print(alelle1_index)
				# print(all_peaks[i])
				if all_peaks[i] == alelle1_index:
					label.append(1)
				else:
					label.append(0)
					noise_count += 1

	else:
		continue


dataset = {'Filename': file_name, 'Label': label, 'Area': area_of_peaks, 'Length': length_in_bp, 'Height of peaks': height_of_peaks}

df = pd.DataFrame(dataset, columns=['Filename', 'Label', 'Area', 'Length', 'Height of peaks'])
df.to_csv('training.csv', index=False)	# 'Hsc40_improved.csv'
