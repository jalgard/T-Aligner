library("edgeR")
dest_files <- list.files('./', pattern='.dest$')
print(dest_files)
dest_samples <- strsplit(dest_files, '.', fixed = TRUE)
dest_samples <- matrix(unlist(dest_samples),ncol=2,byrow=TRUE)[,1]
args <- commandArgs(trailingOnly = TRUE)
conditions <- args[1:length(args)-1]

raw_data <- read.table(dest_files[1], sep='\t', header = FALSE, quote = "", row.names= 1 )

for(i in 2:length(dest_files))
{
    next_row <- read.table(dest_files[i], sep='\t', header = FALSE, quote = "")
    raw_data[,i] = next_row[,2]
}

y <- DGEList(counts = raw_data, group = conditions)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
#results_edgeR <- topTags(et, n = nrow(raw_data), sort.by = "none")
print(topTags(et))
write.table(topTags(et, n = nrow(raw_data), sort.by = "none"), file=args[length(args)], sep='\t', quote=FALSE, col.names=FALSE)

pdf(file = 'DETEST.pdf', width = 5, height = 5)
edge_data <- na.omit(topTags(et, n = nrow(raw_data), sort.by = "none"))
#print(as.vector(edge_data[,2]))
#plot(as.numeric(unlist(edge_data[,1])), -1*log(as.numeric(unlist(edge_data[,4]))))
plot(as.numeric(unlist(edge_data[,1])), as.numeric(unlist(edge_data[,2])))
dev.off()
