########download the data
le download.list |perl -ne 'chomp;print "ascp -P33001 -O33001 -QT -L- -l 1000M -i /home/suisha/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@$_ .\n"'>download.ebi.sh
le filereport_read_run_PRJEB23500_tsv.txt |perl -ne 'chomp;@t=split/\t/; $md=$t[11];$file=$t[12];@f=split/;/,$file;@ff1=split/\//,$f[0];@ff2=split/\//,$f[1];@mdt=split/;/,$md;print "$mdt[0]\t$ff1[-1]\n$mdt[1]\t$ff2[-1]\n";'>md5.ebi
le filereport_read_run_PRJEB23500_tsv.txt |perl -ne 'chomp;@t=split/\t/;print "$t[-2],/mnt/microbiota/suisha/Re27/fastq/$t[5]_1.fastq.gz,forward\n$t[-2],/mnt/microbiota/suisha/Re27/fastq/$t[5]_2.fastq.gz,reverse\n"' > path.list

## trim
cd Analysis
le path.list|perl -e '<>;while(<>){chomp;@t=split/,/; if(exists $h{$t[0]}){}else{$h{$t[0]}=1;print "cutadapt -Z -g CCTACGGGNGGCWGCAG -a GACTACHVGGGTATCTAATCC -o 01.trim/$t[0].1.fastq.gz -p 01.trim/$t[0].2.fastq.gz $t[1] "; $t[1] =~ s/_1.fastq.gz/_2.fastq.gz/g; print "$t[1]\n";}}'>trim.sh
sh trim.sh
perl -e 'open I,"path.list";$h= <I>;print $h;while(<I>){chomp;@t=split/,/;if(exists $h{$t[0]}){}else{$h{$t[0]}=1;print "$t[0],/mnt/microbiota/suisha/Re27/Analysis/01.trim/$t[0].1.fastq.gz,forward\n$t[0],/mnt/microbiota/suisha/Re27/Analysis/01.trim/$t[0].2.fastq.gz,reverse\n";}}' >path_trim.list
##import the data into qiime2
qiime tools import --input-path path_trim.list --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33 --type 'SampleData[PairedEndSequencesWithQuality]'
qiime demux summarize --i-data paired-end-demux.qza --o-visualization demux-summary.qzv

##dada2
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 268 --p-trunc-len-r 203 --o-table table.dada.268.203.qza --o-representative-sequences rep-seqs.dada.268.203.qza --o-denoising-stats denoising-stats.dada.268.203.qza --p-n-threads 0

qiime feature-table tabulate-seqs --i-data rep-seqs.dada.268.203.qza --o-visualization rep-seqs.dada.268.203.qzv
qiime feature-table summarize   --i-table table.dada.268.203.qza   --o-visualization table.dada.268.203.qzv
qiime metadata tabulate   --m-input-file denoising-stats.dada.268.203.qza --o-visualization denoising-stats.dada.268.203.qzv

## alpha diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.dada.268.203.qza --o-alignment aligned-rep-seqs.dada.qza --o-masked-alignment masked-aligned-rep-seqs.dada.qza --o-tree unrooted-tree.dada.qza --o-rooted-tree rooted-tree.dada.qza
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.dada.qza   --i-table table.dada.268.203.qza --p-sampling-depth 6040 --m-metadata-file metadata.txt --output-dir core-metrics-results_dada
qiime diversity alpha-rarefaction --i-table core-metrics-results_dada/rarefied_table.qza --i-phylogeny rooted-tree.dada.qza --m-metadata-file metadata.txt --o-visualization rare.dada.qzv --p-max-depth 6040 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}
qiime diversity alpha --i-table core-metrics-results_dada/rarefied_table.qza --p-metric 'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results_dada/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza

##taxonomy
qiime feature-classifier classify-sklearn   --i-classifier Re27.ref-seqs.classifier.qza --i-reads rep-seqs.deblur317.qza --o-classification taxonomy.qza
qiime metadata tabulate   --m-input-file taxonomy.qza  --o-visualization taxonomy.qzv
qiime feature-classifier classify-sklearn   --i-classifier /home/suisha/database/silva/99/Re27.ref-seqs.classifier.qza --i-reads rep-seqs.dada.268.203.qza --o-classification taxonomy.dada.268.203.qza&

qiime taxa barplot   --i-table core-metrics-results_dada/rarefied_table.qza --i-taxonomy taxonomy.dada.268.203.qza --m-metadata-file metadata.txt --o-visualization taxarare-bar-plots.qzv

##pathway
#conda activate qiime2-2021.2
qiime picrust2 full-pipeline    --i-table   table.dada.268.203.qza  --i-seq  rep-seqs-dada2.qza   --output-dir rep-seqs.dada.268.203.qza    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime feature-table summarize   --i-table q2-picrust2_output/pathway_abundance.qza  --o-visualization q2-picrust2_output/pathway_abundance.qzv
 qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 372680 --m-metadata-file metadata.txt --output-dir pathabun_core_metrics_out --p-n-jobs 12
qiime tools export    --input-path pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
 biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

qiime diversity core-metrics    --i-table q2-picrust2_output/ko_metagenome.qza --p-sampling-depth 3524278 --m-metadata-file metadata.txt --output-dir KO_core_metrics_out --p-n-jobs 36
 qiime tools export    --input-path KO_core_metrics_out/rarefied_table.qza    --output-path ko_exported
 biom convert    -i  ko_exported/feature-table.biom -o pathabun_exported/ko_metagenome.txt --to-tsv

