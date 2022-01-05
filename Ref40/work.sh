####import data to qiime2

le filereport_read_run_PRJNA375065_tsv_sampleid.txt |perl -ne 'chomp;@t=split/\t/;print "cutadapt -Z -g CCTACGGGNGGCWGCAG -a GACTACHVGGGTATCTAATCC -o 01.trim/$t[-1].1.fastq.gz -p 01.trim/$t[-1].2.fastq.gz /mnt/microbiota/suisha/Re40/fastq/$t[3]_1.fastq.gz /mnt/microbiota/suisha/Re40/fastq/$t[3]_2.fastq.gz\n"'>trim.sh
sh trim.sh
 perl -e 'open I,"metadata.txt";print "sample-id,absolute-filepath,direction\n";<I>;while(<I>){chomp;@t=split/\t/;print "$t[0],/mnt/microbiota/suisha/Re40/01.trim/$t[0].1.fastq.gz,forward\n$t[0],/mnt/microbiota/suisha/Re40/01.trim/$t[0].2..fastq.gz,reverse\n"}'>import.list

qiime tools import --input-path import.list --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33 --type 'SampleData[PairedEndSequencesWithQuality]'
qiime demux summarize --i-data paired-end-demux.qza --o-visualization demux-summary.qzv

####dada2
qiime dada2 denoise-paired   --i-demultiplexed-seqs paired-end-demux.qza  --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 233 --p-trunc-len-r 214 --o-table table.dada.qza --o-representative-sequences rep-seqs.dada.qza --o-denoising-stats denoising-stats.dada.qza --p-n-threads 0&
qiime feature-table filter-features --i-table table.dada.qza --p-min-frequency 2 --o-filtered-table table_nosingletons.data.qza
qiime vsearch join-pairs  --i-demultiplexed-seqs paired-end-demux.qza --o-joined-sequences demux-paired-join.qza
qiime metadata tabulate   --m-input-file denoising-stats.dada.qza --o-visualization denoising-stats.dada.qzv

###diversity
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.dada.qza --o-alignment aligned-rep-seqs.dada.qza --o-masked-alignment masked-aligned-rep-seqs.dada.qza --o-tree unrooted-tree.dada.qza --o-rooted-tree rooted-tree.dada.qza
qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.dada.qza   --i-table table.dada.qza   --p-sampling-depth 10670 --m-metadata-file metadata.reformat.txt --output-dir core-metrics-results-dada

qiime diversity alpha --i-table core-metrics-results-dada/rarefied_table.qza --p-metric 'chao1' --o-alpha-diversity chao1_alpha_diverisyt.dada.qza
qiime diversity alpha --i-table core-metrics-results-dada/rarefied_table.qza --p-metric 'simpson' --o-alpha-diversity simpson_alpha_diverisyt.qza

###taxonomy
qiime feature-classifier classify-sklearn   --i-classifier Re40.ref-seqs.classifier.qza   --i-reads rep-seqs.dada.qza   --o-classification taxonomy.qza
qiime taxa barplot   --i-table core-metrics-results-dada/rarefied_table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.reformat.txt --o-visualization taxarare-bar-dada.plots.qzv

###pathway
conda activate qiime2-2021.2
qiime picrust2 full-pipeline    --i-table table.dada.qza    --i-seq rep-seqs.dada.qza    --output-dir q2-picrust2_output    --p-placement-tool epa-ng   --p-threads 12    --p-hsp-method mp    --p-max-nsti 2 --verbose
qiime diversity core-metrics    --i-table q2-picrust2_output/pathway_abundance.qza    --p-sampling-depth 1209732    --m-metadata-file metadata.reformat.txt    --output-dir pathabun_core_metrics_out    --p-n-jobs 12
qiime tools export    --input-path  pathabun_core_metrics_out/rarefied_table.qza --output-path pathabun_exported
biom convert    -i pathabun_exported/feature-table.biom    -o pathabun_exported/feature-table.biom.tsv    --to-tsv

