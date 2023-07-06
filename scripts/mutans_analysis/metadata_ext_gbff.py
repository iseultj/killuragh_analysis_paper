#!/usr/bin/python
# quick script to extract metadata from gbff files

# usage: python metadata_ext_gbff.py *gbff > ${out}.csv
import sys
from Bio import SeqIO
gb_files = sys.argv[1:]

print("ID,Strain,Isolation_Source,Host,Date,Country")
for gb_file in gb_files:
	file_id = gb_file.split("/")[0]


	# parse genbank file
	for gb_record in SeqIO.parse(gb_file, "genbank") :
    	# now do something with the record
        	for feat in gb_record.features:
                	if feat.type == 'source':
                        	source = gb_record.features[0]
# reset all variables
				country = 'not_reported'
				strain = 'not_reported'
				iso_source = 'not_reported'
				host = 'not_reported'
				date = 'not_reported'
                        	for qualifiers in source.qualifiers:
#			    	print(qualifiers)
                            		if qualifiers == 'country':
                                		country = str(source.qualifiers['country'][0])
                            			if country == '':
                                			country = 'none_reported'
				#	else:
				#		country = 'not_reported'
			    		if qualifiers == 'strain':
						strain = str(source.qualifiers['strain'][0])
				#	else:
				#		strain = 'not_reported'
			    		if qualifiers == 'isolation_source':
						iso_source = str(source.qualifiers['isolation_source'][0])
				#	else:
				#		iso_source = 'not_reported'
			    		if qualifiers == 'host':
						host = str(source.qualifiers['host'][0])
				#	else:
				#		host = 'not_reported'
			    		if qualifiers == 'collection_date':
						date = str(source.qualifiers['collection_date'][0])
				#	else:
				#		date = 'not_reported'
#			    if qualifiers == 'organism':
#				org = str(source.qualifiers['organism'][0])

        	
        	print(file_id +  ',' +  strain + ',' + iso_source + ',' + host + ',' + date +  ',' + country)

