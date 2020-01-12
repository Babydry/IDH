#!/usr/bin/env python
# Copyright 2018 BARRIOT Roland
# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
from os.path import isfile
import numpy as np
from scipy.stats import binom, hypergeom, chi2_contingency

# SCRIPT PARAMETERS
# e.g. ./blastset.py --sets EcolA.biocyc.sets --query 'ALAS ARGS ASNS ASPS CYSS GLTX GLYQ GLYS HISS ILES'
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets (categories).')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage. chi2 and coverage are NOT YET IMPLEMENTED')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
param = parser.parse_args()

class ComparedSet(object):
	def __init__(self, id, name = '', common = 0, size = 0, pvalue = 1, elements = [], common_elements = []):
		self.id = id
		self.name = name
		self.common = common
		self.size = size
		self.pvalue = pvalue
		self.elements = elements
		self.common_elements = common_elements

# LOAD QUERY
text = param.query
query = set()
if isfile(text):
	with open(text) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if l!='':
				query |= set(l.split())
else: # parse string
	query |= set(text.split())

# LOAD REFERENCE SETS
def load_sets(filename):
	sets = {}
	ids = set()
	with open( filename ) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			words = l.split('\t')
			if len(words) > 2 and not words[0].startswith('#'):
				id = words.pop(0)
				name = words.pop(0)
				words = set(words)
				sets[ id ] = { 'name': name, 'elements': words}
				ids |= words
	return [ sets, len( ids ) ]
(sets, population_size) = load_sets(param.sets)

# EVALUATE SETS
results = []
query_size = len(query)
for id in sets:
	elements = sets[ id ][ 'elements' ]
	common_elements = elements.intersection( query )
	if param.measure=='binomial': # binom.cdf(>=success, attempts, proba)
		# p_success = 384/2064, 152 attempts, 61 success
		#~ pval = binom.pmf(61, 152, 384.0/2064)
		pval = binom.cdf( query_size - len(common_elements), query_size, 1 - float(len(elements))/population_size)
	elif param.measure=='coverage':
		#create coverage statistic
		pval = 1 - (len(common_elements)/len(elements))*(len(common_elements)/query_size)
	elif param.measure=='chi2':
		#create contingency table
		cont_tab=np.array([[len(common_elements),query_size-len(common_elements)],[len(elements)-len(common_elements),population_size-query_size-len(elements)+len(common_elements)]])
		pval = chi2_contingency(cont_tab)[1]
	elif param.measure=='hypergeometric': # hypergeom.sf(common-1, population, target, query) = 1-p( X <= x-1 ) = p( X >= x )
		pval = hypergeom.sf(len(common_elements)-1, population_size, len(elements), query_size)
	else:
		print('sorry, %s not (yet) implemented' % ( param.measure ))
		exit(1)
	r = ComparedSet( id, sets[id]['name'], len(common_elements), len(elements), pval, elements, common_elements)
	results.append( r )

# PRINT SIGNIFICANT RESULTS
results.sort(key=lambda an_item: an_item.pvalue)
i=1
for r in results:
	# FDR
	if param.adjust and r.pvalue > param.alpha * i / len(results): break
	# limited output
	if param.limit > 0 and i>param.limit: break
	# alpha threshold
	elif r.pvalue > param.alpha and param.measure!='coverage' : break
	# Dont show when c/q * c/t  = 0
	if param.measure == 'coverage' and r.pvalue == 1 : break
	# OUTPUT
	print("%s\t%s\t%s/%s\t%s\t%s" % ( r.id, r.pvalue, r.common, r.size, r.name, ', '.join(r.common_elements)))
	i+=1
