import sys
import os
import gzip
import fileinput

args = sys.argv
# process vcf
for each in fileinput.input():

        each = each.strip()
        # split vcf to fields
        col = each.split()
        # print header
        if each.startswith('#'):
                print(each)
        else:
            # header cols: CHR-REF (unchanged)
            static_cols = '\t'.join(col[0:4])
            alt = col[4].split(',')
            # qual filter info and format are constant always but might as well save these also
            static2_cols = '\t'.join(col[5:9])
            #### Processing multiallelic sites
	    # do not allow triallelic sites
	    if len(alt) > 2:
	    	continue
            elif len(alt) > 1:
                if alt[0] == "*":
                    newalt = alt[1]
                    newgt = []
                    for gt in col[9:]:
                        if gt == "0":
                            newgt.append(gt)
                        elif gt == "1":
                            newgt.append(".")
                        elif gt == "2":
                            newgt.append("1")
                    
                    # now you've changed encoding of all genotypes, want to print this modified line to new VCF
                    newgt_print = '\t'.join(newgt)
                    print(static_cols + '\t' + newalt + '\t' + static2_cols + '\t' + newgt_print)
                elif alt[1] == "*":
                    newalt = alt[0]
                    newgt = []
                    for gt in col[9:]:
                        if gt == "2":
                            newgt.append(".")
                        else:
                            newgt.append(gt)
                    newgt_print = '\t'.join(newgt)
                    print(static_cols + '\t' + newalt + '\t' + static2_cols + '\t' + newgt_print)
                # don't allow other multiallelic sites????
                else:
                    continue

            else:
                print(each) # unchanged record
