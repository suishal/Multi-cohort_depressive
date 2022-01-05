###trim and flash
#example of trim
#cutadapt -Z -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o 01.trim/103.1.fastq.gz -p 01.trim/103.2.fastq.gz study_12382_raw_100520-200346/per_sample_FASTQ/88054/103_S1_L001_R1_001.fastq.gz Re77/study_12382_raw_100520-200346/per_sample_FASTQ/88054/103_S1_L001_R2_001.fastq.gz
#example of flash
##/home/suisha/bin/FLASH-1.2.11-Linux-x86_64/flash 103.1.fastq 103.2.fastq --allow-outies -o 103

###import data to qiime2
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path flash.import.list --input-format SingleEndFastqManifestPhred33 --output-path single-end-demux.qza

###dada2
nohup qiime dada2 denoise-single --i-demultiplexed-seqs single-end-demux.qza --p-trim-left 0 --p-trunc-len 401 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza --p-n-threads 0

##diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs2.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table-dada2.qza   --p-sampling-depth 8886 --m-metadata-file metadata.txt --output-dir core-metrics-results
qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file metadata.txt --o-visualization rare.qzv --p-max-depth 8886 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric  'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza

##taxonomy
nohup qiime feature-classifier classify-sklearn   --i-classifier  Ref27.classifier_qiime19.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza&
qiime taxa barplot   --i-table core-metrics-results/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxarare-bar-plots.qzv

##pathway
qiime picrust2 full-pipeline    --i-table table.dada.qza    --i-seq rep-seqs.dada.qza    --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 869885    --m-metadata-file metadata.txt    --output-dir pathabun_core_metrics_out    --p-n-jobs 12
qiime tools export    --input-path  pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

