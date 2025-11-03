# Create clust tables gene expression


clust_events_pladb_exons<-na.omit(data.frame(ID=pladb_exons$PSI$EVENT, control=rowMeans(pladb_exons$PSI[,7:9]), pladb_1mm=rowMeans(pladb_exons$PSI[,10:12]),pladb_10mm=rowMeans(pladb_exons$PSI[,13:15])))
write_tsv(clust_events_pladb_exons, "clust_events_splicing_pladb_onlyexons.txt")

clust_events_pladb_introns<-na.omit(data.frame(ID=pladb_introns$PSI$EVENT, control=rowMeans(pladb_introns$PSI[,7:9]), pladb_1mm=rowMeans(pladb_introns$PSI[,10:12]),pladb_10mm=rowMeans(pladb_introns$PSI[,13:15])))
write_tsv(clust_events_pladb_introns, "clust_events_splicing_pladb_onlyintrons.txt")
# Read clust output


# Read the TSV file
# Skip the first row (for descriptive cluster information) and specify the second row as column names
clust_out <- read.delim("Clusters_Objects.tsv", header = FALSE, skip = 1, stringsAsFactors = FALSE)[-1,]

# Set proper column names using the first row of the file
column_names <- read.delim("Clusters_Objects.tsv", header = FALSE, nrows = 1, stringsAsFactors = FALSE)
colnames(clust_out) <- column_names
clust_out[clust_out == ""] <- NA

clustered_exons<-read_tsv("results/clust_results_tao_exons/Clusters_Objects.tsv")[-1,]