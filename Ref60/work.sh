####import the data into qiime2
nohup qiime tools import --input-path import.list.match_id.csv --output-path paired-end-demux.match.qza --input-format PairedEndFastqManifestPhred33 --type 'SampleData[PairedEndSequencesWithQuality]'
qiime metadata tabulate --m-input-file match_sample.metadata.txt_rn --o-visualization match_sample.metadata.qzv
qiime demux summarize --i-data paired-end-demux.match.qza --o-visualization demux-summary.qzv

####dada2
nohup qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.match.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 220 --p-trunc-len-r 218 --o-table table.denoise.qza --o-representative-sequences rep-seqs.dada.qza --o-denoising-stats denoising-stats.dada.qza --p-n-threads 0
qiime feature-table tabulate-seqs --i-data rep-seqs.dada.qza --o-visualization rep-seqs.dada.qzv
qiime feature-table summarize   --i-table table.denoise.qza   --o-visualization table.denoise.qzv
qiime metadata tabulate   --m-input-file denoising-stats.dada.qza --o-visualization denoising-stats.dada.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.dada.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

####diversity
nohup qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza --i-table table.denoise.qza --p-sampling-depth 7718 --m-metadata-file match_sample.metadata.txt_rn --output-dir core-metrics-results_7718 --p-n-jobs 8&

qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file match_sample.metadata.txt_rn --o-visualization rare.dada.qzv --p-max-depth 7719 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}

qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza

####taxonomy
nohup qiime feature-classifier classify-sklearn   --i-classifier Re60.ref-seqs.classifier.qza --i-reads rep-seqs.dada.qza --o-classification taxonomy.qza

qiime taxa barplot   --i-table core-metrics-results_7718/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file match_sample.metadata.txt_rn --o-visualization taxarare-bar-plots.qzv

####pathway
seqtk seq -r d47700ac-1993-4558-91b1-50b7b558a3d3/data/dna-sequences.fasta >rev.dna-sequences.fasta
picrust2_pipeline.py -s rev.dna-sequences.fasta -i 935c95c9-6382-49d0-86d4-558694d685ac/data/feature-table.biom -o picrust2_out_pipeline -p 32
picrust2 full-pipeline    --i-table  table.denoise.qza   --i-seq  rev.dna-sequences.qza   --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 1808387    --m-metadata-file ../match_sample.metadata.txt_rn    --output-dir pathabun_core_metrics_out    --p-n-jobs 24


