# Main file for AuSOME-ML
# Input: Name of fsa file
# Output: Suggested peaks with corresponding length (in bp)
# Workflow:
# 			1. Load model created from Hsc40_model.py
# 			2. Read the input fsa file (access channel 1; the data for Hsc 40 loci, and channel 105; data for the size standard/ladder)
# 			3. Detect all peaks in the file within the range of the size standard/ladder (35 bp - 500 bp)
# 			4. Create a region (-40 +40) around each peak from which the 3 features will be extrcted
# 			5. Feed each of the regions in the model
# 			6. Mark the peaks for which the model returned a value of 1 (true peak) and display the corresponding length (in bp)

from Bio import SeqIO
from ladder_fit import convert_to_bp, convert_to_index, find_lower, find_upper
from findpeaks import findpeaks as fp
import numpy as np
import pandas as pd
import pickle
from matplotlib import pyplot as plt

model_name = 'Hsc40_model_RandomForest_improved.sav'														# load model
model = pickle.load(open(model_name, 'rb'))

file = pd.read_csv('test.csv')
a = 'DATA1'
b = 'DATA2'
c = 'DATA3'
d = 'DATA4'
e = 'DATA105'
dye = [500, 490, 450, 400, 350, 340, 300, 250, 200, 160, 150, 139, 100, 75, 50, 35]

label = []
area_of_peaks = []
length = []
peaks = []
height = []
correct = 0
total = 0

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

		data = record.annotations['abif_raw'][a]
		ladder = record.annotations['abif_raw'][e]

		if alelle1 == alelle2:
			total += 1
		else:
			total += 2

		index_35 = find_lower(ladder, dye)																# the index of the 35 bp peak
		index_500 = find_upper(ladder, dye)																# the index of the 500 bp peak

		detected_peaks = fp.findpeaks(data, spacing=25, limit=25)

		for i in range(len(detected_peaks)-1):
			lower_end_of_window = detected_peaks[i] - 40												# Create region from which the features will be extracted
			upper_end_of_window = detected_peaks[i] + 40
			if lower_end_of_window < index_35 or upper_end_of_window > index_500:						# ignore peaks outside the range of the ladder
				continue

			peaks = detected_peaks[i]																	# store the index of the peak in an array
			area_of_peaks = np.trapz(data[lower_end_of_window:upper_end_of_window])						# get area under the curve in the region 'lower_end_of_window' to 'upper_end_of_region'
			length = round(convert_to_bp(detected_peaks[i], ladder, dye))
			height = data[detected_peaks[i]]

			label = model.predict([[area_of_peaks, length, height]])									# the height of the peak (in RFU)	note: remember 'data' is an array of length ~7000 where each value of the element corresponds to intensity (height)

			if label == 1:
				if length >= alelle1-2 and length <= alelle1+2 or length >= alelle2-2 and length <= alelle2+2:
					correct += 1
				# else:
				# 	correct -= 1


		# for i in range(len(peaks)-1):
		# 	label.append(model.predict([[area_of_peaks[i], length[i], height[i]]]))		# feed the model the features of each region


		# print(len(label))
		# print(len(length))
		# for x in range(len(label)-1):
		# 	if label[x] == 1:
		# 		if length[x] in [alelle1-1, alelle1+1, alelle2-1, alelle2+1]:
		# 			correct += 1


accuracy = correct/total
print(f'\nTotal no. of correct calls is:	{correct}')
print(f'Total no. of calls:		{total}')
print(f'\nThe accuracy is:		{accuracy}')

# for x in range(len(label)-1):
# 	if label[x] == 1:
# 		plt.plot(peaks[x], height[x], 'o', markersize=3, color='red', label=str(int(length[x])) + ' bp')			# Mark the peak for which the model classified as true peak and display the corresponding length (in bp)
# 		print(length[x])


# plt.plot(data, color='blue')
# plt.plot(ladder, color='black')
# plt.xlabel('Index', fontsize=16)
# plt.ylabel('RFU (Relative Fluorescent Units)', fontsize=16)
# plt.legend()
# plt.show()