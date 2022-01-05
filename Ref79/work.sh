###import the data into qiime2

split_libraries.py -m ANMAP.txt -f rgsbt.an.fna -q cizhw.an.qual -b 10 -d --record_qual_scores -o Split_Library_Output/
convert_fastaqual_fastq.py -f Split_Library_Output/seqs.fna -q Split_Library_Output/seqs_filtered.qual -o fastq_files/
split_sequence_file_on_sample_ids.py -i fastq_files/seqs.fastq --file_type fastq -o split_fastq_files
find split_fastq_files/*|perl -e 'print "sample-id,absolute-filepath,direction\n";while(<>){chomp;@t=split/\/|\./;print "$t[-3],$_,forward\n"}'>import.list

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path import.list --input-format SingleEndFastqManifestPhred33 --output-path single-end-demux.qza
qiime demux summarize  --i-data single-end-demux.qza --o-visualization demux-summary.qzv

####dada2
qiime dada2 denoise-single --i-demultiplexed-seqs single-end-demux.qza --p-trim-left 0 --p-trunc-len 404 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza

qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs-dada2.qzv
qiime feature-table summarize   --i-table table-dada2.qza --o-visualization table-dada2.qzv
qiime metadata tabulate   --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv

####alpha diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs2.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
le ANMAP.txt|perl -ne 'chomp; @t=split/\t/; print $t[0],"\t",substr($t[0],0,3),"\t",$t[4],"\n"' >metadata.txt
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table-dada2.qza --p-sampling-depth 4000 --m-metadata-file metadata.txt --output-dir core-metrics-results
qiime diversity alpha-rarefaction --i-table core-metrics-results/rarefied_table.qza --i-phylogeny rooted-tree.qza --m-metadata-file metadata.txt --o-visualization rare.qzv --p-max-depth 4000 --p-metrics {'simpson_e','chao1','observed_otus','faith_pd','shannon'}

qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.qza
qiime diversity alpha --i-table core-metrics-results/rarefied_table.qza --p-metric  'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza

####taxonomy
qiime feature-classifier classify-sklearn   --i-classifier Re79.ref-seqs.classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza
qiime taxa barplot   --i-table core-metrics-results/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization taxarare-bar-plots.qzv

#pathway
qiime picrust2 full-pipeline    --i-table table-dada2.qza    --i-seq rep-seqs-dada2.qza    --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 869885    --m-metadata-file metadata.txt    --output-dir pathabun_core_metrics_out    --p-n-jobs 12
qiime tools export    --input-path  pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

