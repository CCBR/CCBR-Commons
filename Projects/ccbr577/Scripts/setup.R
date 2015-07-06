# setup.R
# Convert files, and generally set up data for densityScan
# Randall Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


######### read in mm9 genes and convert to bed format #########
tmp <- read.csv('../Builds/mm9G37_genes.csv', stringsAsFactors = FALSE, na.strings = '-')

tmp$chrNumber <- with(tmp, ifelse(Chromosome.Name == 'X', 20,
                           ifelse(Chromosome.Name == 'Y', 21,
                           ifelse(Chromosome.Name == 'MT', 22, as.numeric(Chromosome.Name)))))

tmp$Chromosome.Name <- paste('chr', tmp$Chromosome.Name, sep = '')

tmp <- subset(tmp, !is.na(tmp$chrNumber))

tmp <- tmp[order(tmp$chrNumber, tmp$Gene.Start..bp.),]

tmp <- unique(tmp)
write.table(unique(subset(tmp, select = c(Chromosome.Name, Gene.Start..bp., Gene.End..bp.,Ensembl.Gene.ID))),
            file = '../Builds/mm9G37_genes.bed', sep = '\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE)


######### Add Ensembl ID / GO ID Translation table #########

cat("Ensembl\tGO\tMGI\n", file = '../Builds/mm9G37_genes.txt')
write.table(unique(subset(tmp, select = c(Ensembl.Gene.ID, GO.Term.Accession, MGI.symbol))),
            file = '../Builds/mm9G37_genes.txt', sep = '\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
