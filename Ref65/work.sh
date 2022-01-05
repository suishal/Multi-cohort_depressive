 
####import data
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path import.file --input-format SingleEndFastqManifestPhred33 --output-path single-end-demux.qza
nohup qiime demux summarize --i-data single-end-demux.qza --o-visualization demux-summary.qzv --p-n 2562532

###dada2
nohup qiime dada2 denoise-single   --i-demultiplexed-seqs single-end-demux.qza   --p-trim-left 0   --p-trunc-len 281 --o-representative-sequences rep-seqs-dada2.281.qza   --o-table table-dada2.281.qza --o-denoising-stats stats-dada2.281.qza --p-n-threads 8&

qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs-dada2.qzv
qiime feature-table summarize   --i-table table-dada2.qza --o-visualization table-dada2.qzv
qiime metadata tabulate   --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv

###taxonomy
nohup qiime feature-classifier classify-sklearn   --i-classifier Re65.ref-seqs.classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza

###alpha diversity
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table-dada2.qza   --p-sampling-depth 4562 --m-metadata-file metadata.txt --output-dir core-metrics-results
qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file metadata.txt --o-visualization rare.qzv --p-max-depth 4562 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}

###pathway
qiime picrust2 full-pipeline    --i-table table.dada.qza    --i-seq rep-seqs.dada.qza    --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth  494271    --m-metadata-file metadata.txt    --output-dir pathabun_core_metrics_out    --p-n-jobs 12
qiime tools export    --input-path pathabun_core_metrics_out/rarefied_table.qza    --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

