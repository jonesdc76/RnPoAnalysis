#!/bin/bash
IBD=0
RUNLIST=2018C_GoodRuns.txt
DATA_DIR=/p/lustre2/psptexp/user/jonesdc/data/Analyzed/
#DATA_DIR=/p/lustre2/psptexp/converted/analyzed/
RELEASE=Phys_2018C
ANARELEASE=Analyzed_2018C
FILE=AD1_BiPo_DT.root
FILE=AD1_RnPo.root
if [ $IBD -eq 1 ];then
    FILE=AD1_Wet_PhysL.root
fi
#SSHSOCKET=~/.ssh/jones291@borax.llnl.gov
#ssh -M -f -N -o ControlPath=$SSHSOCKET jones291@borax.llnl.gov
ssh -M -f -N borax
for i in $(cat $RUNLIST);do
    dir=$BIPO_OUTDIR/$ANARELEASE/$i
    if [ ! -d $dir ];then
	mkdir -p $dir
    fi
    echo "Copying $i"
    #scp borax:${DATA_DIR}/${ANARELEASE}/${i}/${FILE} ${dir}/AD1_BiPo_wDT.root
    scp -p borax:${DATA_DIR}/${ANARELEASE}/${i}/${FILE} ${dir}/${FILE}
    #scp oslic:${DATA_DIR}/${ANARELEASE}/${i}/${FILE} ${dir}/${FILE}
done
#ssh -S $SSHSOCKET -O exit jones291@borax.llnl.gov
