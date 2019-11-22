# Necessary functions for extracting data in fsa files

# 	ladder_sizing		: Finds the corresponding indexes of each peak in the ladder/size standard ex.: 35 bp peak - 1632, 50 bp peak - 1845, 75 bp peak - 2031 ...
# 	find_lower			: Finds the index of the first peak in the ladder (or the 35-bp allele)
# 	find_upper			: Finds the index of the last peak in the ladder (or the 500-bp allele)
# 	convert_to_bp		: Takes index as input, returns the corresponding length in base pair as output
# 	convert_to_index	: Takes length (in bp) as input, returns the corresponding index as output

from findpeaks import findpeaks

# LIZ_500 is the data for the ladder: [500, 490, 450, 400, 350, 340, 300, 250, 200, 160, 150, 139, 100, 75, 50, 35]

def ladder_sizing(data, LIZ_500):
	data_105 = list(data)	# the data for the channel 105 of the fsa file containing data for the ladder/size standard. Array of length ~7000 with each element corresponding to value of intensity (in Relative Fluorescent Units) 
	indexes = findpeaks.findpeaks(data_105, spacing=35, limit=100)	# returns array containing indexes of peaks detected in data_105 (ladder)
	ind = []				# initialize array that will contain the indexes of the ladder
	j = len(indexes) - 1	# stop when all elements in 'indexes' have been checked
	k = 0					# stop when 'ind' contains 16 elements (16 peaks in the ladder)

	# read 'indexes' starting from the last element (start from last since data_105(ladder) contains a lot of noise in the first few indexes)
	while j >= 0:
		if k == 15:								# stop search for 16th element in 'ind', tune if needed for other size standards
			ind.append(indexes[j])
			break
		counter = 1
		next_index = j - counter
		del_ratio = (indexes[j] - indexes[next_index])/(LIZ_500[k] - LIZ_500[k+1])	# gets the ratio of the difference between the current index and the next index to the difference between the corresponding value of that index and the next index in LIZ_500 
		in_range = del_ratio > 8.5 and del_ratio < 13.5								# true if 'del_ratio' is within 8.5 to 13.5  note: arbitrary value that works for the 500 LIZ ladder
		if in_range:
			ind.append(indexes[j])													# append 'indexes[j]' in 'ind' if 'in_range' satisfied
			k += 1																	# add 1 to 'ind' element counter
			j -= 1																	# move to next element in 'indexes'
		else:
			while del_ratio < 8.5 or del_ratio > 13.5:								# find the next pair of indexes that will satisfy 'in_range'
				counter += 1
				next_index = j - counter
				if next_index < 0:
					j -= 1
					counter = 1
					next_index = j - counter
				del_ratio = (indexes[j] - indexes[next_index])/(LIZ_500[k] - LIZ_500[k+1])
			ind.append(indexes[j])
			j = next_index
			k += 1

	return ind

def find_lower(data, LIZ_500):
	ind = ladder_sizing(data, LIZ_500)			# get the indexes of the ladder

	if ind[0] == ind[len(ind)-1]:				# ladder_sizing sometimes returns 'ind' with the same first and last element  hint: has something to do with 8.5 - 13.5 threshold range or sometimes ladder data (data_105) too noisy e.g. multiple peaks
		ind.pop(len(ind)-1)

	return ind[len(ind)-1]						# return last element of 'ind' which corresponds to the index of the 35 bp peak

def find_upper(data, LIZ_500):
	ind = ladder_sizing(data, LIZ_500)

	if ind[0] == ind[len(ind)-1]:
		ind.pop(len(ind)-1)

	return ind[0]								# return first element of 'ind' which corresponds to the index of the 500 bp peak

def convert_to_bp(index, data, LIZ_500):
	ind = ladder_sizing(data, LIZ_500)

	for c in range(0, 16):
		if index > ind[c] or index == ind[c]:
			bp_pred = ((index - ind[c])/((ind[c-1] - ind[c])/(LIZ_500[c-1] - LIZ_500[c]))) + LIZ_500[c]	# contact igdemavivas@up.edu.ph or gaabrasaldo@up.edu.ph for explanation of this line
			break
		else:
			bp_pred = 0							# returns 0 if something wrong

	return bp_pred

def convert_to_index(bp, data, LIZ_500):
	ind = ladder_sizing(data, LIZ_500)

	for c in range(0, 16):
		if bp > LIZ_500[c]:
			index_pred = ((bp - LIZ_500[c])/((LIZ_500[c-1] - LIZ_500[c])/(ind[c-1] - ind[c]))) + ind[c]
			break

	return index_pred
