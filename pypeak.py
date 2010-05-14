#!/usr/bin/python

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#       
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#       
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.

import csv
import operator
from optparse import OptionParser

def main():
	# option parser
	usage = 'python pypeak.py <ip_bed_file> <control_bed_file>'
	parser = OptionParser(usage=usage)
	parser.add_option("-w", "--window_length",\
	type="int", dest="window_length", default="135",\
	help="scanning window size in bp for each strand: default=135")
	parser.add_option("-s", "--signal_threshold",\
	type="int", dest="signal_threshold", default="5",\
	help="signal threshold calling enriched regions: default=5")
	parser.add_option("-o", "--out_file",\
	type="string", dest="out_file", default="pypeaks_out.bed",\
	help="name of the bed file output is written to")
	parser.add_option("-q", "--quiet",
	action="store_false", dest="verbose", default=True,
	help="don't print status messages to stdout")
	parser.add_option("-d", "--dip_threshold",
	type="float", dest="dip_threshold", default="0.6",\
	help="dip threshold for splitting enriched regions: default=0.6")
	
	# read arguments and options
	(options, args) = parser.parse_args()
	if len(args) != 2:
		parser.error('Please provide an ip and control file!') 
	ip_bed = args[0]
	control_bed = args[1]
	
	# parse alignments
	print_status('Parsing IP alignments ...', options.verbose)
	ip_dict = build_coverage_dict(ip_bed)
	print_status('Parsing input alignments ...', options.verbose)
	control_dict = build_coverage_dict(control_bed)
	
	# subtract background
	print_status('Calculating fold enrichment ...', options.verbose)
	enrichment_dict = build_enrichment_dict(ip_dict, control_dict)
	
	# find peaks and store in a dictionary
	peak_dict = {}
	for chr in enrichment_dict.keys():
		print_status('Finding peaks on chromosome %s' % chr, options.verbose)
		plus_strand = enrichment_dict[chr]['+']
		minus_strand = enrichment_dict[chr]['-']
		peak_dict[chr] = peak_scanner(plus_strand, minus_strand, chr, options)
	
	# create outfile and write peaks to it
	f = open(options.out_file, 'w')
	f.close()
	for chr in peak_dict.keys():
		for i in peak_dict[chr]:
			f = open(options.out_file, 'a')
			f.write('%s\t%d\t%d\t%s\t%.2f\n' % i)
			f.close()

def find_peaks(peak_region, dip_threshold):
	# extract peaks from enriched regions
	s = [i[0] for i in peak_region]
	p = [i[1] for i in peak_region]
	ds = diff(s)
	l = len(ds)
	maxima = []
	# find all local maxima
	for i in range(l-1):
		# the sign of the slope changes from + we have a max or inflect
		if ds[i] > 0 and ds[i+1] <= 0:
			maxima.append((s[i+1], i+1))
	# if enriched region is too small for max finding just return first
	if maxima == []:
		tiny_peak = peak_region[0]
		return [tiny_peak]
	# sort maxima in descending order
	maxima.sort()
	maxima.reverse()
	# the first peak is defined as the highest maximum
	highest_max = maxima[0]
	first_peak = (highest_max[0], p[highest_max[1]])
	# position of first peak within enriched region
	max_local_pos = highest_max[1]
	for i in maxima[1:]:
		second_max_score = i[0]
		second_max_pos = i[1]
		# find the deepest dip between the two peaks
		if second_max_pos < max_local_pos:
			dip = min(s[second_max_pos:max_local_pos])
		else:
			dip = min(s[max_local_pos:second_max_pos])
		# return 2 peaks if dip is deep enough
		if dip/second_max_score < dip_threshold:
			second_peak = (i[0], p[i[1]])
			return [first_peak, second_peak]
	return [first_peak]

def peak_scanner(plus_strand, minus_strand, chr, options):
	# scan alignment, extract enriched regions and find peaks
	l = options.window_length
	t = options.signal_threshold
	peaklist = []
	peak_region = []
	count = 0
	# scan both strands with shifted double windows
	for i in range(0 + l, len(plus_strand) - l):
		plus_enr = sum(plus_strand[i-l:i])
		minus_enr = sum(minus_strand[i:i+l])
		peak_enr = (plus_enr + minus_enr) / (2.0 * l)
		# extract enriched regions and call peak finding function
		if peak_enr > t:
			peak_region.append((peak_enr,i))
		else:
			if len(peak_region) > 0:
				peaks = find_peaks(peak_region, options.dip_threshold)
				for i in peaks:
					count += 1
					peak = (chr, i[1] - l, i[1] + l,'%s_Peak_%d' % (chr, count), i[0])
					peaklist.append(peak)
				peak_region = []
	return peaklist
	
def build_enrichment_dict(ip_dict, c_dict):
	# create dict with subtracted enrichment values
	enr_dict = {}
	# adjust chromosome length in both dictionaries
	for chr in ip_dict.keys():
		enr_dict[chr] = {}
		for strand in ip_dict[chr].keys():
			ip_length = len(ip_dict[chr][strand])
			c_length = len(c_dict[chr][strand])
			if ip_length < c_length:
				ip_dict[chr][strand].extend([0]*(c_length - ip_length))
			else:
				c_dict[chr][strand].extend([0]*(ip_length - c_length))
	# subtract enrichment value from one from the other
	for chr in ip_dict.keys():
		enr_dict[chr]['+'] = map(operator.sub, ip_dict[chr]['+'], c_dict[chr]['+'])
		enr_dict[chr]['-'] = map(operator.sub, ip_dict[chr]['-'], c_dict[chr]['-'])
	return enr_dict
	
def build_coverage_dict(file_name):
	# parse a bed file and create a strand specific coverage dict
	cov_dict = {}
	for i in csv.reader(open(file_name), delimiter='\t'):
		chr = i[0]
		start = int(i[1])
		end = int(i[2])
		strand = i[5]
		# if chromosome is not in the dict, required keys
		if not chr in cov_dict:
			cov_dict[chr] = {}
			cov_dict[chr]['+'] = []
			cov_dict[chr]['-'] = []
		# dynamically adjust chromosome length if not long enough
		chr_len = len(cov_dict[chr][strand])
		if 	chr_len < end + 1:
			cov_dict[chr][strand].extend([0]*(end - chr_len + 1))
		# add alignment to the computed coverage
		for i in range(start, end + 1):
			cov_dict[chr][strand][i] += 1
	return cov_dict

def print_status(string, boolean):
	if boolean:
		print(string)

def diff(x):
	d = []
	for i in range(len(x)-1):
		d.append(x[i+1] - x[i])
	return d

if __name__ == '__main__':
	main()
