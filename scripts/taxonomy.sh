#!/bin/bash

REFERENCE_FASTA=$1

bash downloadTaxonomy.sh

if [ ! -f ../Taxonomy/all.accession2taxid ]
        then 
			cat ../Taxonomy/*.accession2taxid > ../Taxonomy/all.accession2taxid
fi

echo "Getting accession numbers from reference FASTA file"
echo ""

#python checker.py $REFERENCE_FASTA ../Taxonomy/accessions.txt
grep ">" $REFERENCE_FASTA | awk '{print $1}' | awk -F ">" '{print $2}' > ../Taxonomy/accessions.txt

echo "Finding accession number TaxIDs"
echo ""

python accessions2taxids.py ../Taxonomy/ ../Taxonomy/accessions.txt ../Taxonomy/all.accession2taxid

echo "Removing unnecassary files"
echo ""

#rm ../Taxonomy/accessions.txt
rm ../Taxonomy/dead_nucl.accession2taxid
rm ../Taxonomy/dead_wgs.accession2taxid
rm ../Taxonomy/nucl_gb.accession2taxid
rm ../Taxonomy/nucl_wgs.accession2taxid
