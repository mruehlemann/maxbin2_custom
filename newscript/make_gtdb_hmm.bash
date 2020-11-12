GTDB_DIR="/work_ifs/ikmb_repository/shared/microbiome/GTDB/release95"

while read line; do
read group id <<< $(echo $line)
cat $GTDB_DIR/markers/$group/individual_hmms/$id
done < bac120_hmm.tsv > gtdb_bac120.hmm

while read line; do
read group id <<< $(echo $line)
cat $GTDB_DIR/markers/$group/individual_hmms/$id
done < ar122_hmm.tsv > gtdb_ar122.hmm
