#######data filter and download
le filereport_read_run_PRJEB11419.ref23_tsv.txt |perl -ne 'chomp;@t=split/\t/;if($t[7]=="feces metagenome" || $t[7]=="gut metagenome" || $t[7]=="human gut metagenome"){print $t[12],"\n";}'> download.list
le download.list |perl -ne 'chomp;print "ascp -P33001 -O33001 -QT -L- -l 1000M -i /home/suisha/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@$_ .\n"'>download.ebi.sh
le filereport_read_run_PRJEB11419.ref23_tsv.txt |perl -ne 'chomp;@t=split/\t/;if($t[7]=="feces metagenome" || $t[7]=="gut metagenome" || $t[7]=="human gut metagenome"){@f=split/\//,$t[12];print $t[11],"\t",$f[-1],"\n"}' >md5.ebi

#metadata
le 10317_20200724-214900.txt|perl -ne 'chomp;@t=split/\t/;$n=@t;if($n == 667){print $_,"\n";}'>full.record.list
 grep "feces" full.record.list>feces.full.record.list
le feces.full.record.list|perl -ne 'chomp;@t=split/\t/;print "$t[0]\t$t[7]\t$t[41]\t$t[71]\t$t[91]\t$t[332]\t$t[360]\t$t[219]\n"'>feces.full.record.list.pick
le match_sample_metadata.txt |perl -e '%h;open I ,"match_sample_metadata.txt";<I>;while(<I>){chomp;@t=split/\t/;print `grep $t[0] study.txt `}'>pick.ebi.list


#only download the pick match sample
le match_sample_metadata.txt |perl -ne 'chomp;@t=split/\t/;print `grep $t[0] filereport_read_run_PRJEB11419_tsv.txt`'>pick.download.list

le pick.download.list |perl -ne '$i="ascp -P33001 -O33001 -QT -L- -l 1000M -i /home/suisha/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp\@";chomp;@t=split/\t/;if($t[8]=~/;/){@tt=split/;/,$t[8]; print $i, $tt[0], " .\n", $i, $tt[1], " .\n"}else{print $i, $t[8]," .\n";}' >download.matchsample.sh

le match_sample_metadata.txt |perl -ne 'chomp;@t=split/\t/;print `grep $t[0] filereport_read_run_PRJEB11419_tsv_md5.txt `'>pick.download.md5.list

########import the data into qiime2
find final_fastq/*|perl -e 'print "sample-id\tabsolute-filepath\n";while(<>){@t=split/\/|\./;print "$t[6].$t[7]\t$_"}' >import.list
le import.list |perl -e '$h=<>;chomp $h; $h=~s/\t/,/g;print "$h,direction\n";while(<>){chomp;$_=~s/\t/,/g;print "$_,forward\n"}'>import.file

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path import.file --input-format SingleEndFastqManifestPhred33 --output-path single-end-demux.qza

########dada2
nohup qiime dada2 denoise-single   --i-demultiplexed-seqs single-end-demux.qza   --p-trim-left 0   --p-trunc-len 151   --o-representative-sequences rep-seqs-dada2.qza   --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza&
qiime dada2 denoise-single   --i-demultiplexed-seqs single-end-demux.qza   --p-trim-left 0   --p-trunc-len 151   --o-representative-sequences rep-seqs-dada2.qza   --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza --p-n-threads 0
qiime demux summarize  --i-data single-end-demux.qza --o-visualization demux-summary.qzv --p-n 4615429
qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs-dada2.qzv
qiime feature-table summarize   --i-table table-dada2.qza --o-visualization table-dada2.qzv
qiime metadata tabulate   --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs2.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

########alpha diversity
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.dada.qza   --p-sampling-depth 8050 --m-metadata-file match_sample_metadata.mvtitle.txt --output-dir core-metrics-results

qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file match_sample_metadata.mvtitle.txt --o-visualization rare.qzv --p-max-depth 8050 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}

qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric  'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza

########taxonomy
nohup qiime feature-classifier classify-sklearn   --i-classifier Re23.ref-seqs.classifier.qza   --i-reads rep-seqs-dada2.qza   --o-classification taxonomy.qza
qiime taxa barplot   --i-table core-metrics-results/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxarare-bar-plots.qzv

#######pathway
#conda activate qiime2-2021.2
qiime picrust2 full-pipeline    --i-table  table-dada2.qza   --i-seq  rep-seqs-dada2.qza   --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime feature-table summarize   --i-table q2-picrust2_output/pathway_abundance.qza  --o-visualization q2-picrust2_output/pathway_abundance.qzv
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 673777 --m-metadata-file match_sample_metadata.mvtitle.txt --output-dir pathabun_core_metrics_out --p-n-jobs 12
qiime tools export    --input-path pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

