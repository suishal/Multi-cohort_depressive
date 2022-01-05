
########trim the adaptor
 le path.list |perl -e '<>;while(<>){chomp;@t=split/,/; if(exists $h{$t[0]}){}else{$h{$t[0]}=1;print "cutadapt -Z -g AGAGTTTGATCCTGGCTCAG -a CTGCTGCCTYCCGTA -o $t[0].1.fastq.gz -p $t[0].2.fastq.gz $t[1] "; $t[1] =~ s/_1.fastq.gz/_2.fastq.gz/g; print "$t[1]\n";}}'>trim.sh
 sh trim.sh
########import the data into qiime2
 perl -e 'open I,"../sample_id_from_metadata.txt";print "sample-id,absolute-filepath,direction\n";while(<I>){chomp;@t=split/,/;print "$t[0],$t[0].1.fastq.gz,forward\n$t[0],$t[0].2.fastq.gz,reverse\n"}'>import.list
 qiime tools import --input-path import.list --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33 --type 'SampleData[PairedEndSequencesWithQuality]'
########dada2
qiime dada2 denoise-paired   --i-demultiplexed-seqs paired-end-demux.qza  --p-trim-left-f 0   --p-trim-left-r 0   --p-trunc-len-f 273   --p-trunc-len-r 236   --o-table table.dada.qza   --o-representative-sequences rep-seqs.dada.qza   --o-denoising-stats denoising-stats.dada.qza --p-n-threads 0 &
qiime feature-table tabulate-seqs --i-data rep-seqs.dada.qza --o-visualization rep-seqs.dada.qzv
 qiime feature-table summarize   --i-table table.dada.qza   --o-visualization table.dada.qzv
qiime metadata tabulate   --m-input-file denoising-stats.dada.qza --o-visualization denoising-stats.dada.qzv
########alpha diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.dada.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.dada.qza   --p-sampling-depth 12154   --m-metadata-file metadata.txt --output-dir core-metrics-results
qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file metadata.txt --o-visualization rare.qzv --p-max-depth 12154 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric  'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza

########taxonomy
qiime feature-classifier classify-sklearn   --i-classifier Re1.ref-seqs.classifier.qza   --i-reads rep-seqs.dada.qza   --o-classification taxonomy.qza
qiime taxa barplot   --i-table core-metrics-results/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxarare-bar-plots.qzv

qiime metadata tabulate   --m-input-file taxonomy.qza  --o-visualization taxonomy.qzv

########pathway
conda activate qiime2-2021.2
qiime picrust2 full-pipeline    --i-table table.dada.qza    --i-seq rep-seqs.dada.qza    --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 1408013    --m-metadata-file metadata.txt    --output-dir pathabun_core_metrics_out    --p-n-jobs 12
qiime tools export    --input-path  pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv


