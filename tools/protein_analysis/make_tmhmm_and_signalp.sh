#!/bin/sh
echo "This will create a tar-ball suitable to upload to the toolshed."

if [ -f "tools/protein_analysis/make_tmhmm_and_signalp.sh" ]
then
echo "Good, in the expected directory"
else
echo "ERROR. Run this from the galaxy root directory."
exit 1
fi

if [ -f "tmhmm_and_signalp.tar.gz" ]
then
echo "ERROR. File tmhmm_and_signalp.tar.gz already exists."
exit 1
fi

#Create tar file with core XML wrappers
if [ -f "tmhmm_and_signalp.tar" ]
then
rm tmhmm_and_signalp.tar
fi

#Create tar file (-cf then -rf to add to it)
tar -cf tmhmm_and_signalp.tar tools/protein_analysis/LICENSE
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/README
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/suite_config.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/seq_analysis_utils.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/signalp3.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/signalp3.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/tmhmm2.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/tmhmm2.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/promoter2.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/promoter2.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/psortb.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/psortb.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/wolf_psort.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/wolf_psort.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/rxlr_motifs.xml
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/rxlr_motifs.py
tar -rf tmhmm_and_signalp.tar tools/protein_analysis/whisson_et_al_rxlr_eer_cropped.hmm
tar -rf tmhmm_and_signalp.tar test-data/four_human_proteins.fasta
tar -rf tmhmm_and_signalp.tar test-data/four_human_proteins.signalp3.tabular
tar -rf tmhmm_and_signalp.tar test-data/four_human_proteins.tmhmm2.tabular
tar -rf tmhmm_and_signalp.tar test-data/four_human_proteins.wolf_psort.tabular
tar -rf tmhmm_and_signalp.tar test-data/empty.fasta
tar -rf tmhmm_and_signalp.tar test-data/empty_tmhmm2.tabular
tar -rf tmhmm_and_signalp.tar test-data/empty_signalp3.tabular
tar -rf tmhmm_and_signalp.tar test-data/empty_psortb_terse.tabular
tar -rf tmhmm_and_signalp.tar test-data/empty_wolf_psort.tabular
tar -rf tmhmm_and_signalp.tar test-data/empty_promoter2.tabular
tar -rf tmhmm_and_signalp.tar test-data/k12_ten_proteins.fasta
tar -rf tmhmm_and_signalp.tar test-data/k12_ten_proteins_psortb_p_terse.tabular
tar -rf tmhmm_and_signalp.tar test-data/rxlr_win_et_al_2007.fasta
tar -rf tmhmm_and_signalp.tar test-data/rxlr_win_et_al_2007.tabular
tar -rf tmhmm_and_signalp.tar test-data/rxlr_win_et_al_2007_sp3.tabular
tar -rf tmhmm_and_signalp.tar test-data/Adenovirus.fasta
tar -rf tmhmm_and_signalp.tar test-data/Adenovirus.promoter2.tabular

#Compress the tar file
gzip tmhmm_and_signalp.tar

#Check the output
echo "Expect a tar-ball 34 files, have:"
tar -tzf tmhmm_and_signalp.tar.gz | wc -l
