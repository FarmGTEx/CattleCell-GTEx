#!/bin/sh
###Set workspace
work_path=$1
type=$2
Tissues=$3
cd ${work_path}

###Cycle submit task
#for tissue in ${Tissues[*]}
#do
	for ttt in {1..29}
	do
		awk 'NR>1{print $0}' /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/chr${ttt}/eQTL.CD4Tcells_Model_training_chr${ttt}_covariances.txt >> /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/chr.tmp #Remove column names and merge files
	done
	sed -i '1i\GENE RSID1 RSID2 VALUE' /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/chr.tmp #Add column name
	gzip -c /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/chr.tmp > /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/eQTL.CD4Tcells_Model_training_covariances.txt.gz
	rm /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/TWAS/Blood/CD4Tcells/covariances/chr.tmp #Delete temporary files
#done
