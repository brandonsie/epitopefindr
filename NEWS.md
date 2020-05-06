# epitopefindr 1.1.27 (2020-05-06)
## Minor changes
- epitopeSummary output table sorting modification. group number coerced to numeric.

# epitopefindr 1.1.26
## Breaking changes
- output table changes (fixed bug with epitope_summary$id number of rows mismatch, outputTable() function now returns a list with epitope_key and epitope_summary tables as elements, epfind() writes csvs with "NA" to fill empty cells)  

## Minor changes
- added parameter `remove.intermediates = FALSE` to epfind(). If set to true, then the /intermediate_files/ folder is deleted from epitopefindr's output directory at the end of the run.

# epitopefindr 1.1.25
## Minor changes
- Added support for command line utility pdfunite as an alternative to pdftk.
- Removed some defunct code and comments.
- Beginnings of framework for parallelization.