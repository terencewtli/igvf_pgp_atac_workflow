library(ArchR)


proj<- loadArchRProject(path="/pathway/ArchR_Project_directory/")

cell.metadata <- as.data.frame(Object_data@meta.data)
cell.metadata$Cell <-row.names(cell.metadata)
proj@cellColData$Cell<- gsub("#","_",row.names(proj@cellColData) )
proj@cellColData$Name <- cell.metadata[match( proj@cellColData$Cell,cell.metadata$Cell),]$Cluster

proj <- proj[proj@cellColData$Cell %in% cell.metadata$Cell,]

library(BSgenome.Hsapiens.UCSC.hg38)

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Name")
pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
                                 ArchRProj = proj,
                                   groupBy = "Name",
                                     pathToMacs2 = pathToMacs2
                                 )

proj <- addPeakMatrix(proj)
proj@cellColData$Name <- as.character(proj@cellColData$Name )

markersPeaks <- getMarkerFeatures(
                                    ArchRProj = proj,
                                      useMatrix = "PeakMatrix",
                                        groupBy = "Name",
                                          bias = c("TSSEnrichment", "log10(nFrags)"),
                                            testMethod = "wilcoxon"
                                    )

peak.info<-as.data.frame(markersPeaks@assays@data)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

##Collect diff peaks
gene.peak.info<-NULL
for(Names in names(markerList)){
      print(Names)
  peak.data<-as.data.frame(markerList[[Names]])
    if(nrow(peak.data)==0){next}
    peak.data$chr <- peak.data$seqnames
      peak.data$Sample <- Names
      gene.peak.info<-rbind(gene.peak.info, peak.data)
}


write.table(gene.peak.info, paste0(outdir, Exp.Name, ".FullArchR_NamedCluster_DiffPeaks.table.txt"),
                        sep="\t", quote=F)
