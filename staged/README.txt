This is a placeholder directory with stub values for a minimal set of
data in order to get the app running (globals.Rdata, bag data in tsne/
and gene names in markers/).

Override the value of dropviz.prep.dir in options.R and either

- run prep-all.R to write data to prep.dir
- set prep.dir to a directory containing existing prep'd data
- run prep-local.sh to copy a subset from a remote prep'd data dir
