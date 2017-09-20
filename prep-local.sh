#!/bin/sh
#
# When developing locally, syncronize a subset of the atlas and staged data from a remote host
# for a subset of the data.

REMOTE=dropviz.dizz.org:/mnt/disks/dropseq
REMOTE_STAGED=${REMOTE}/staged
REMOTE_ATLAS=${REMOTE}/atlas_ica
LOCAL=/cygdrive/d/dropviz
LOCAL_STAGED=${LOCAL}/staged
LOCAL_ATLAS=${LOCAL}/atlas_ica

EXPERIMENTS="GRCm38.81.P60Hippocampus GRCm38.81.P60Striatum"

mkdir -p ${LOCAL_STAGED}/expr
mkdir -p ${LOCAL_STAGED}/metacells
mkdir -p ${LOCAL_STAGED}/tsne
mkdir -p ${LOCAL_ATLAS}

rsync -v ${REMOTE_STAGED}/expr/gene.map.RDS ${LOCAL_STAGED}/expr/
rsync -rav ${REMOTE_STAGED}/markers ${LOCAL_STAGED}/
rsync -v ${REMOTE_STAGED}/tsne/\*.Rdata ${LOCAL_STAGED}/tsne

for exp in ${EXPERIMENTS}; do
    mkdir -p ${LOCAL_STAGED}/expr/${exp}
    rsync -v ${REMOTE_STAGED}/expr/${exp}/genes.RDS ${LOCAL_STAGED}/expr/${exp}/
    rsync -rav ${REMOTE_STAGED}/expr/${exp}/gene ${LOCAL_STAGED}/expr/${exp}/
    rsync -rav ${REMOTE_STAGED}/metacells/${exp}\* ${LOCAL_STAGED}/metacells/
    rsync -rav ${REMOTE_STAGED}/tsne/${exp} ${LOCAL_STAGED}/tsne/

    mkdir -p ${LOCAL_ATLAS}/${exp}
    rsync -rav ${REMOTE_ATLAS}/${exp}/components ${LOCAL_ATLAS}/${exp}/
done

# create exp_sets.txt
cat > exp_sets.txt <<EOF
exp.label	exp.title	exp.dir
EOF

for exp in ${EXPERIMENTS}; do
    short=`echo $exp | sed -e 's/GRCm38.81.P60//'`
    echo "${exp}\t${short}\t${LOCAL_ATLAS}/F_${exp}\n"
done

cat >> options.R <<EOF
options(dropviz.prep.dir='${LOCAL_STAGED}')
options(dropviz.experiments='exp_sets.txt')
EOF

# generate a globals.RDS that only includes EXPERIMENTS
R -e "source('prep-global.R')"
