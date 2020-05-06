allergome 1-100
group 3 has pep 66 multiple


rebuild pkgdown when internet

add annotation files to dropbox phip

push to dev version

debug 1:100


epitopefindr log

--

indexGroups line 135 section still allows merge of duplicated
should fix

also some peptides havae NA group, thus fails groupMSA

add progress into to log
add job complete to log
add newline after pregress bars / remove progress bars when they are done


fix indexgroups algorithm

--
epitoep 88 has 768 splits


make indexgorups faster

alinging groups are always the same. 1 from each?. not necessarily


groups.csv only assigning each peptide ot its last group.

also creates thousands of groups
make peptide.once == FALSE by default until ready

try crop out one peptide option

--
make epitopefindr method
write.files = TRUE
if false, skil MSA
skip fwrite in epfind
add return of the three spreadsheets 
epitopekey
epitopesummary
finalalignments
(also initial sequences and starting alignments?)
add groups.csv to output?
make groupMSA compatible with epitopekey
