  ##################################################################
  #generate_coding_genes_annotation_df
  #no arguments 
  #collects all the detected genes from the counts file and connects to biomart to find out the coding genes
  ##################################################################
  generate_coding_genes_annotation_df = function(){
    source("/home/cidr/Documents/work/programs/myUtils.R")
    
    #*****************************************************
    #collect all ensembl ids
    #hardcoded cross files path
    #*****************************************************
    igra_negative_cross_counts_dir = "/home/cidr/Documents/work/fresh/automatedCounts/cross_sectional_igra_negative"
    igra_negative_cross_counts_files = list.files(path = igra_negative_cross_counts_dir , pattern = "[.]csv$")
    
    #*****************************************************************************
    #ensembl_gene_ids will contain all the ensemble ids from all the counts files
    #*****************************************************************************
    ensembl_gene_ids = c()
    
    #the below loop picks up all the ensemble ids from the igra_negative cross ection files
    for(file in igra_negative_cross_counts_files)
    {
      filepath = getPath(filename = file, directory = igra_negative_cross_counts_dir)
      counts_df = read.csv(file = filepath,stringsAsFactors = FALSE)
      ensembl_gene_ids = c(ensembl_gene_ids,counts_df[,1])
    }
    
    #*****************************************************************************
    #the below loop picks up all the ensembl ids from the tb cross section files
    #hardcoded cross files path
    #*****************************************************************************
    tb_cross_counts_dir = "/home/cidr/Documents/work/fresh/automatedCounts/cross_sectional_tb/"
    tb_cross_counts_files = list.files(path = tb_cross_counts_dir , pattern = "[.]csv$")
    
    for(file in tb_cross_counts_files)
    {
      filepath = getPath(filename = file, directory = tb_cross_counts_dir)
      counts_df = read.csv(file = filepath,stringsAsFactors = FALSE)
      ensembl_gene_ids = c(ensembl_gene_ids,counts_df[,1])
    }
    
    
    ensembl_gene_ids = unique(ensembl_gene_ids)
    
    #*****************************************************************************
    #cleaning the ensembl ids removes the numbers after the point in the id
    #example ENSG12345678.10 will me made ENSG12345678
    #*****************************************************************************
    clean_ensembl_gene_ids = getCleanEnsembleIds(ensembl_gene_ids)
    clean_ensembl_gene_ids = unique(clean_ensembl_gene_ids)
    
    #*****************************************************************************
    #connecting to biomart for getting annotation and saving in df
    #*****************************************************************************
    all_protein_coding_genes_df = getAnnotationFromBiomart(ids = clean_ensembl_gene_ids, filter = "ensembl_gene_id", columns = c("ensembl_gene_id","hgnc_symbol","description","gene_biotype"))
    
    #*****************************************************************************
    #saving the annotation data
    #hardcoded annotations of all detected genes file path
    #*****************************************************************************
    write.csv(all_protein_coding_genes_df, file = "/home/cidr/Documents/work/fresh/annotation_info/all_gene_details.csv",row.names = FALSE)
    
    #*****************************************************************************
    #hardcoded filepath of all protein coding genes
    #*****************************************************************************
    all_protein_coding_genes_df = filter_annotation_df_for_coding_genes(all_protein_coding_genes_df)
    write.csv(all_protein_coding_genes_df, file = "/home/cidr/Documents/work/fresh/annotation_info/all_coding_genes_details.csv",row.names = FALSE)
  }
  
  
  ##########################################################################
  # remove_non_coding_genes
  # load a dataframe containing all the protein coding genes
  ##########################################################################
  remove_non_coding_genes = function(df, protein_coding_genes_csv_filepath = "/home/cidr/Documents/work/fresh/annotation_info/all_coding_genes_details.csv")
  {
    if(nrow(df)>0)
    {
      protein_coding_genes = read.csv(file = protein_coding_genes_csv_filepath, stringsAsFactors = FALSE)
      coding_genes = protein_coding_genes$hgnc_symbol
      
      df = df[df$hgnc_symbol %in% coding_genes,]
    }
    else
    {
      stop("the data frame is empty")
    }
    return(df)
  }
  ##########################################################################
  # remove_non_coding_genes_from_deseq_result
  # load a dataframe containing all the protein coding genes
  ##########################################################################
  remove_non_coding_genes_from_deseq_result = function(df, protein_coding_genes_csv_filepath = "/home/cidr/Documents/work/fresh/annotation_info/all_coding_genes_details.csv")
  {
    if(nrow(df)>0)
    {
      protein_coding_genes = read.csv(file = protein_coding_genes_csv_filepath, stringsAsFactors = FALSE)
      coding_genes = protein_coding_genes$ensembl_gene_id
      
      df = df[getCleanEnsembleIds(rownames(df)) %in% coding_genes,]
    }
    else
    {
      stop("the data frame is empty")
    }
    return(df)
  }
  
  ##########################################################################
  #filter_annotation_df_for_coding_genes
  #this is a supplementary function and is used to construct a dataframe containing the codind genes annotation
  ##########################################################################
  filter_annotation_df_for_coding_genes=function(annotation_df)
  {
    library("dplyr")
    
    #*****************************************************************************
    #below are the allowed biotypes which represent the coding genes type
    #*****************************************************************************
    allowed_biotypes=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
    
    annotation_df=annotation_df %>% dplyr::filter(gene_biotype %in% allowed_biotypes)
    #sum(duplicated(annotation_df$ensembl_gene_id))
    
    return(annotation_df)
  }
  
  ##########################################################################
  #getAnnotationFromBiomart
  #a very critical function for connecting to Biomart and fetching annotation info.
  #filter is the type of id (hgnc_symbol, ensembl_gene_id) you are supplying
  #
  ##########################################################################
  getAnnotationFromBiomart=function(ids, filter, columns, mirror = c("useast", "uswest", "asia", "www"))
  {
    library(biomaRt)
    mirror = match.arg(mirror)
    mart = useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl",mirror = mirror)
    
    ensembl = useDataset(dataset="hsapiens_gene_ensembl",mart=mart)
    #biomartResults=getBM(attributes = columns, filters ="ensembl_gene_id",values = id,mart = mart)
    biomartResults=getBM(attributes = columns, filters =filter,values = ids,mart = mart)
    return(biomartResults)
  }
  
  ##########################################################################
  #get gene ontology annotation for a gene
  #a very critical function for connecting to Biomart and fetching annotation info.
  #filter is the type of id (hgnc_symbol, ensembl_gene_id) you are supplying
  #
  ##########################################################################
  getGeneOntologyAnnotationFromBiomart=function(ids)
  {
    library(biomaRt)
    #mirror = match.arg(mirror)
    filter = "hgnc_symbol"
    columns = c("hgnc_symbol","go_id")
    #mirror = c("useast", "uswest", "asia", "www")
    mirror = "useast"
    mart = useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl",mirror = mirror)
    
    ensembl = useDataset(dataset="hsapiens_gene_ensembl",mart=mart)
    #biomartResults=getBM(attributes = columns, filters ="ensembl_gene_id",values = id,mart = mart)
    biomartResults=getBM(attributes = columns, filters =filter,values = ids,mart = ensembl)
    print(head(biomartResults))
    library(GO.db)
    go_ids = biomartResults$go_id
    #go_annotation = Term(GOTERM)[go_ids[Ontology(go_ids)=="BP"]]
    go_annotation = Term(GOTERM)[go_ids]
    biomartResults$go_annotation = go_annotation
    biomartResults = biomartResults %>% dplyr::filter(!is.na(go_annotation))
    getAnnotationString = function(id)
    {
      annot = biomartResults %>% dplyr::filter(hgnc_symbol == id) %>% dplyr::select(go_annotation) %>% pull() %>% paste(collapse = ",")
      return(annot)
    }
    go_annotation = sapply(ids,getAnnotationString)
    names(go_annotation) = ids
    print(go_annotation)
    #library(GO.db)
    #go_ids = biomartResults$go_id
    #go_annotation = Term(GOTERM)[go_ids[Ontology(go_ids)=="BP"]]
    #print(head(go_annotation))
    #print(head(go_annotation))
    return(go_annotation)
  }
  
  ##########################################################################
  #getCleanEnsembleIds
  #takes in ensmble ids with the character '.' in it and removes the '.' and everything after it
  ##########################################################################
  getCleanEnsembleIds=function(geneIds)
  {
    geneIds=gsub("[.][0-9]*","",geneIds)
    geneIds=gsub("_[A-Z]*_[A-Z]*","",geneIds)
    return(geneIds)
  }
  
  ##########################################################################
  #getPath
  #return a combined path that is the full path of the file
  ##########################################################################
  getPath=function(filename,directory)
  {
    library(stringr)
    #print(directory)
    #print(filename)
    if(str_detect(directory,"/$"))
    {
      #print("in ")
      directory=str_trunc(directory,nchar(directory)-1,ellipsis="")
    }
    #print(directory)
    return(file.path(directory,filename))
  }
  
  ##########################################################################
  #filterByPvalue
  #filters a differential expression dataframe by p value 
  ##########################################################################
    filterByPvalue = function(df, cutoff = 0.05)
  {
    df=df[!is.na(df$padj),]
    df=df[df$pvalue<cutoff,]
    return(df)
  }
  
  ##########################################################################
  #filterByPadjValue
  #filters a differential expression dataframe by padjusted value 
  ##########################################################################
  filterByPadjValue = function(df, cutoff = 0.05)
  {
    df=df[!is.na(df$padj),]
    df=df[df$padj<cutoff,]
    return(df)
  }
  
  ##########################################################################
  #filterByFoldChange
  #filter differential expression dataframes based on fold change values supplied
  ##########################################################################
  filterByFoldChange = function(df, cutoff)
  {
    df=df[abs(df$log2FoldChange)>=cutoff,]
    return(df)
  }
  
  ##########################################################################
  #filter 
  #a general purpose filter function for filtering differential expression dataframes
  # method specifies based on what you want to filter the dataframe
  #methods takes values (padj, pval)
  ##########################################################################
  filter = function(df, method, cutoff)
  {
    if(method=="padj")
    {
      df=filterByPadjValue(df,cutoff)
    }
    if(method=="pval")
    {
      df=filterByPvalue(df,cutoff)
    }else
    {
      print("wrong method chosen. nothign filtered")
    }
    return(df)
  }
  
  ##########################################################################
  # filterCodingGenesDf
  # the mode specifies whose differential expression data it is 
  # mode takes in 2 values (s, d) 
  # s for shubha
  #d for durbar
  # this was done because the two have generated two different types of output files
  ##########################################################################
  filterCodingGenesDf=function(df,mode="")
  {
    library(biomaRt)
    filterBiotypes = c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
    if(mode=="d")
    {
      genes=as.character(df$gene_symbol)
      if(nrow(df)>0)
      {
        biomartdf=getAnnotationFromBiomart(ids=genes,filter="hgnc_symbol",columns=c("hgnc_symbol","gene_biotype"))
        selectedGenes=biomartdf$hgnc_symbol[biomartdf$gene_biotype %in% filterBiotypes]
        df=df[df$gene_symbol %in% selectedGenes,]
        #print(biomartdf)
        df=df[,c("gene_symbol","description","log2FoldChange","padj","pvalue")]
      }
    }
    else if(mode=="s")
    {
      genes=as.character(df$Symbol)
      if(nrow(df)>0)
      {
        biomartdf=getAnnotationFromBiomart(ids=genes,filter="hgnc_symbol",columns=c("hgnc_symbol","gene_biotype"))
        selectedGenes=biomartdf$hgnc_symbol[biomartdf$gene_biotype %in% filterBiotypes]
        df=df[df$Symbol %in% selectedGenes,]
        #print(biomartdf)
        df=df[,c("Symbol","description","log2FoldChange","padj","pvalue","type_of_gene")]
      }
    }
    else
    {
      stop("supply proper mode. tell whose data this is.")
    }
    return(df)
  }
  
  ##########################################################################
  #filterCodingGenes
  #given a vector of gene ids this function can return the gene ids that are coding
  #the type argument specifies what type of ids these are
  ##########################################################################
  filterCodingGenes=function(genes,gene_id_type = c("ensembl_gene_id","hgnc_symbol"))
  {
    library(biomaRt)
    #*****************************************************************************
    #filterBiotypes store the allowed gene types
    #*****************************************************************************
    gene_id_type = match.arg(gene_id_type)
    filterBiotypes=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
    
    selectedGenes=c()
    if(length(genes)>0)
    {
      #*****************************************************************************
      #needs internet
      #*****************************************************************************
      biomartdf=getAnnotationFromBiomart(ids=genes,filter=gene_id_type,columns=c(gene_id_type,"gene_biotype"))
      selectedGenes=biomartdf[,gene_id_type][biomartdf$gene_biotype %in% filterBiotypes]
    }
    else
    {
      print("no genes supplied. cant figure out why you called me. such a waste of time. bye")
    }
    print(selectedGenes)
    return(selectedGenes)
  }
  
  ##########################################################################
  # writePathwaysToSheet
  # this function needs fixes in the way it presents the columns.. the style of columns
  # tobefixed
  ##########################################################################
  writePathwaysToSheet=function(wb,sheetName,pathwayresult)
  {
    sheet <- createSheet(wb, sheetName = sheetName)
    s=CellStyle(wb)+Alignment(horizontal = "ALIGN_CENTER",wrapText=T)
    requiredData=as.data.frame(pathwayresult)
    if(nrow(requiredData)>0)
    {
      if(str_detect(sheetName,"ora_kegg") | str_detect(sheetName,"enricher"))
      {
        requiredData=subset(requiredData,select=c("ID","Description","geneID","pvalue","p.adjust"))
      }
      else if(str_detect(sheetName,"gse_GO"))
      {
        head(requiredData)
      }
      else
      {
        requiredData=subset(requiredData,select=c("ID","Description","core_enrichment","pvalue","p.adjust"))
      }
      requiredData = add_how_many_genes_enriched(requiredData)
      setColumnWidth(sheet, colIndex=c(2,3), colWidth=40)
      setColumnWidth(sheet, colIndex=c(1,4:ncol(requiredData)), colWidth=20)
      addDataFrame(requiredData,sheet = sheet,row.names = F)
    }
    else
    {
      print("cant add anything to sheet")
    }
  }
  
  ##########################################################################
  #compare_and_get_degs
  #this function compares two differential expression dataframes 
  #how the data is to be compared is decided by the "onlyIn" argument
  #"set1" will return the degs which are present only in set 1
  #"set2" will return the degs which are present only in set 2
  #"common" returns the common degs 
  ##########################################################################
  compare_and_get_degs = function(degs1_df, degs2_df, onlyIn = c("set1","set2","common"))
  {
    df = data.frame()
    #setting the set of choice
    onlyIn = match.arg(onlyIn)
    
    if(onlyIn == "set1"){
      
      selectedRows = !(degs1_df$gene_symbol %in% degs2_df$gene_symbol)
      df = degs1_df[selectedRows,]
      
    }else if(onlyIn == "set2"){
      
      selectedRows = !(degs2_df$gene_symbol %in% degs1_df$gene_symbol)
      df = degs2_df[selectedRows,]
      
    }else if(onlyIn == "common"){
      
      selectedRows = (degs1_df$gene_symbol %in% degs2_df$gene_symbol)
      df = degs1_df[selectedRows,]
      
    }else{
      stop("wrong value supplied for arg onlyIn")
    }
    
    return(df)
  }
  
  
  ##########################################################################
  #getKeggPathwaysOnlyIn
  #takes in two sets of differential expression results
  #return a set of pathways 
  #the results depend on the "onlyIn" argument
  #"set1" gets pathways only because of set1 degs
  #"set2" gets pathways only because of set2 degs
  #"common" gets common pathways 
  ##########################################################################
  getKeggPathwaysOnlyIn = function(set1_df,set2_df,onlyIn = c("set1","set2","common"))
  {
    library(stringr)
    df = data.frame()
    #setting the set of choice
    onlyIn = match.arg(onlyIn)
    print("in getKeggPathwaysOnlyIn")
    id_column_name = "ID"

    rows_set1_df = nrow(set1_df)
    rows_set2_df = nrow(set2_df)
    
    if(is.null(rows_set1_df))
    {
      rows_set1_df = 0
    }
    if(is.null(rows_set2_df))
    {
      rows_set2_df = 0
    }
    if(rows_set1_df == 0 & rows_set2_df == 0)
    {
      return(df)
    }
    else if(rows_set1_df > 0 & rows_set2_df > 0)
    {
      print("rows_set1_df > 0 & rows_set2_df > 0")
      print(head(set1_df))
      print(head(set2_df))
      print(id_column_name)
      pathways_in_set1 = set1_df[,id_column_name]
      pathways_in_set2 = set2_df[,id_column_name]
    }
    else if(rows_set1_df > 0 & rows_set2_df == 0){
      print("rows_set1_df > 0 & rows_set2_df == 0")
      pathways_in_set1 = set1_df[,id_column_name]
      pathways_in_set2 = c()
    }
    else if(rows_set1_df == 0 & rows_set2_df > 0){
      print("rows_set1_df == 0 & rows_set2_df > 0")
      print(set2_df)
      print(nrow(set2_df))
      pathways_in_set1 = c()
      pathways_in_set2 = set2_df[,id_column_name]
    }
    print(onlyIn)
    if(onlyIn == "set1")  
    {
      if(rows_set1_df == 0)
      {
        return(df)
      }
      selected_rows = !(pathways_in_set1 %in% pathways_in_set2)
      print("enter onlyIn == \"set1\"")
      df = set1_df[selected_rows,]
      print("exit onlyIn == \"set1\"")
    }else if(onlyIn == "set2")
    {
      if(rows_set2_df == 0)
      {
        return(df)
      }
      selected_rows = !(pathways_in_set2 %in% pathways_in_set1)
      
      df = set2_df[selected_rows,]
    }else if(onlyIn == "common"){
      if(rows_set1_df == 0 | rows_set2_df == 0)
      {
        return(df)
      }
      selected_rows = pathways_in_set1 %in% pathways_in_set2
      df = set1_df[selected_rows,]
    }else{
      stop("wrong argument provided. you have to mention which set do you want ? set1 or set2")
    }
    print("exiting function getKeggPathwaysOnlyIn")
    return(as.data.frame(df))
  }
  
  
  ##########################################################################
  #get_common_pathways
  #takes in two sets of pathways results and returns a set of pathways that are common in both the set
  #note that this function doesnt take the differential expression results as input
  ##########################################################################
  get_common_pathways = function(set1_df,set2_df)
  {
    df = data.frame()
    print("in get_common_pathways")
    id_column_name = "ID"
    
    pathways_in_set1 = set1_df[,id_column_name]
    pathways_in_set2 = set2_df[,id_column_name]
    
    selected_rows = pathways_in_set1 %in% pathways_in_set2
    
    df = set1_df[selected_rows,]
    
    return(df)
  }
  
  ########################
  #get gene ontologies 
  ########################
  getGeneOntologies=function(genes, level = 0, ontology_type = "BP")
  {
    library(clusterProfiler)
    library(org.Hs.eg.db)
    #gmtFilePath="/home/cidr/Documents/work/programs/wikipathways-20190710-gmt-Homo_sapiens.gmt"
    ###############################################################
    # something before analysis
    ###############################################################
    #wp2gene = clusterProfiler::read.gmt(gmtfile=gmtFilePath)
    #wp2gene = wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    #wpid2gene = wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
    #wpid2name = wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
    #**********************************************************
    #performing gene set enrichment analysis
    #**********************************************************
    #entrzId=clusterProfiler::bitr(names(geneList),fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    set.seed(1234)
    #geneList=sort(geneList,decreasing = T)
    #**********************************************************
    #performing gse in kegg
    #**********************************************************
    res= groupGO(gene     = genes, OrgDb    = org.Hs.eg.db, ont      = ontology_type, level    = level, readable = TRUE)
    
    return(as.data.frame(res))
  }
  
  ########################
  #get kegg pathwaysfrom entrez ids 
  ########################
  getKeggPathways=function(geneList, pathwaysType = c("gsea","enricher"))
  {
    #print(length(geneList))
    #print(geneList)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    print("in getKeggPathways")
    pathwaysType = match.arg(pathwaysType)
    #gmtFilePath="/home/cidr/Documents/work/programs/wikipathways-20190710-gmt-Homo_sapiens.gmt"
    ###############################################################
    # something before analysis
    ###############################################################
    #wp2gene = clusterProfiler::read.gmt(gmtfile=gmtFilePath)
    #wp2gene = wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    #wpid2gene = wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
    #wpid2name = wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
    #**********************************************************
    #performing gene set enrichment analysis
    #**********************************************************
    #entrzId=clusterProfiler::bitr(names(geneList),fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    set.seed(1234)
    geneList=sort(geneList,decreasing = T)
    #**********************************************************
    #performing gse in kegg
    #**********************************************************
    if(pathwaysType == "gsea")
    {
      res= gseKEGG(gene= geneList,organism     = 'hsa',pvalueCutoff = 0.05,nPerm = 10000,seed = T)
      res=setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
      res = as.data.frame(res)
      if(nrow(res)>0){
        res = add_how_many_genes_enriched(res)
        res = res %>% dplyr::select(ID, Description, n, core_enrichment, "p.adjust",  enrichmentScore)
        #head(res)
      }
    }else if(pathwaysType == "enricher"){
      res= enrichKEGG(gene =  names(geneList),organism     = 'hsa',pvalueCutoff = 0.05)
      #print(typeof(res))
      res=setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
      res = as.data.frame(res)
      res = res %>% dplyr::select(ID, Description, Count, geneID, "p.adjust")
    }
    #print(head(res))
    if(is.null(nrow(res)))
    {
      res = data.frame()
    }
    return(res)
  }
  
  ###############################
  #get the reactome pathways just from gene symbols
  ###############################
  getReactomePathways = function(gene_symbols)
  {
    library(org.Hs.eg.db)
    library(dplyr)
    library(ReactomePA)
    entrez_ids = unlist(mapIds(org.Hs.eg.db, gene_symbols, "ENTREZID", "SYMBOL"))
    reactome_df = enrichPathway(gene = entrez_ids)
    reactome_df = as.data.frame(setReadable(reactome_df, org.Hs.eg.db, keyType = "ENTREZID")) %>% dplyr::select(ID, Description, Count, geneID, "p.adjust")
    return(reactome_df)
  }
  
  #get kegg pathways from gene symbols 
  
  get_kegg_pathways_from_gene_symbols = function(gene_symbols)
  {
    library(org.Hs.eg.db)
    library(clusterProfiler)
    entrez_ids = unlist(mapIds(org.Hs.eg.db, gene_symbols, "ENTREZID", "SYMBOL"))
    res= enrichKEGG(gene =  entrez_ids,organism     = 'hsa',pvalueCutoff = 0.05)
    #print(typeof(res))
    res=setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
    res = as.data.frame(res)
    res = res %>% dplyr::select(ID, Description, Count, geneID, "p.adjust")
    return(res)
  }
  
  add_how_many_genes_enriched = function(pathways_df)
  {
    library(dplyr)
    library(stringr)
    print("in add_how_many_genes_enriched")
    print(colnames(pathways_df))
    print(head(pathways_df))
    pathways_df = pathways_df %>% dplyr::mutate(n = str_count(core_enrichment, "/")+1)
    return(pathways_df)
  }
  
  doGSEA=function(geneList,wb,sheetNamePrefix="",backgroundGeneNames)
  {
    print(length(geneList))
    print(geneList)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    gmtFilePath="/home/cidr/Documents/work/programs/wikipathways-20190710-gmt-Homo_sapiens.gmt"
    ###############################################################
    # something before analysis
    ###############################################################
    wp2gene = clusterProfiler::read.gmt(gmtfile=gmtFilePath)
    wp2gene = wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    wpid2gene = wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
    wpid2name = wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
    #**********************************************************
    #performing gene set enrichment analysis
    #**********************************************************
    #entrzId=clusterProfiler::bitr(names(geneList),fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    set.seed(1234)
    geneList=sort(geneList,decreasing = T)
    gsea= GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=T,nPerm=10000,seed = T)
    gsea=setReadable(gsea, org.Hs.eg.db, keyType = "ENTREZID")
    head(gsea)
    #gsea = add_how_many_genes_enriched(as.data.frame(gsea))
    #print(head(gsea))
    writePathwaysToSheet(wb=wb,sheetName=paste(sheetNamePrefix,"gsea",sep=""),gsea)
    
    
    #**********************************************************
    #performing gse in kegg
    #**********************************************************
    res= gseKEGG(gene= geneList,organism     = 'hsa',pvalueCutoff = 0.05,nPerm = 10000,seed = T)
    res=setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
    #res = add_how_many_genes_enriched(as.data.frame(res))
    #print(head(res))
    head(res)
    writePathwaysToSheet(wb=wb,sheetName=paste(sheetNamePrefix,"gse_kegg",sep=""),res)
    #gsea=as.data.frame(gsea)
  }
  
  getCommon=function(ds,sv)
  {
    common=as.character(ds$gene_symbol[ds$gene_symbol %in% sv$Symbol])
    return(common)
  }
  getOnlyInDs=function(ds,sv)
  {
    common=getCommon(ds=ds,sv=sv)
    onlyInDs=as.character(ds$gene_symbol[!(ds$gene_symbol %in% common)])
    return(onlyInDs)
  }
  getOnlyInSv=function(ds,sv)
  {
    common=getCommon(ds=ds,sv=sv)
    onlyInSv=as.character(sv$Symbol[!(sv$Symbol %in% common)])
    return(onlyInSv)
  }
  getCommonDf=function(ds,sv)
  {
    common=as.character(ds$gene_symbol[ds$gene_symbol %in% sv$Symbol])
    ds=ds[ds$gene_symbol %in% common,]
    ds=ds[order(ds$gene_symbol),c("gene_symbol","description")]
    return(ds)
  }
  getOnlyInDsDf=function(ds,sv)
  {
    common=getCommon(ds=ds,sv=sv)
    onlyInDsDf=ds[!(ds$gene_symbol %in% common),]
    onlyInDsDf=onlyInDsDf[order(onlyInDsDf$gene_symbol),]
    return(onlyInDsDf)
  }
  getOnlyInSvDf=function(ds,sv)
  {
    common=getCommon(ds=ds,sv=sv)
    onlyInSvDf=sv[!(sv$Symbol %in% common),]
    onlyInSvDf=onlyInSvDf[order(onlyInSvDf$Symbol),]
    return(onlyInSvDf)
  }
  
  
  #**************************************************************
  # addRowToDetailsDf 
  # function used to add rows to a df ( mostly for creating thge details df of venn diagrams)
  #**************************************************************
  addRowToDetailsDf=function(df,detail,value)
  {
    rownames(df)=NULL
    df=rbind(df,cbind(detail,value))
    return(df)
  }
  
  #**************************************************************
  # getDetailsDf 
  # function used to get the details df. requires just the 
  #**************************************************************
  getDetailsDf=function(set1,set2,names=c())
  {
    common=set1[set1 %in% set2]
    total_common=length(common)
    
    total_set1=length(set1)
    total_set2=length(set2)
    
    only_set1=set1[!(set1 %in% common)]
    total_only_set1=length(only_set1)
    only_set2=set2[!(set2 %in% common)]
    total_only_set2=length(only_set2)
    
    overlap_perc_set1=format(round((total_common/total_set1)*100,2),nsmall = 2)
    overlap_perc_set2=format(round((total_common/total_set2)*100,2),nsmall = 2)
    
    df=data.frame(stringsAsFactors = F)
    df=addRowToDetailsDf(df,detail = paste("total in",names[1]),value = total_set1)
    df=addRowToDetailsDf(df,detail = paste("total in",names[2]),value = total_set2)
    
    df=addRowToDetailsDf(df,detail = paste("common"),value = total_common)
    df=addRowToDetailsDf(df,detail = paste("only in",names[1]),value = total_only_set1)
    df=addRowToDetailsDf(df,detail = paste("only in",names[2]),value = total_only_set2)
    
    df=addRowToDetailsDf(df,detail = paste("overlap of",names[1]),value = paste(overlap_perc_set1,"%"))
    df=addRowToDetailsDf(df,detail = paste("overlap of",names[2]),value = paste(overlap_perc_set2,"%"))
    return(df)
  }
  
  getGenesListForVenn=function(ds,sv)
  {
    common=getCommon
  }
  
  getMeansDf=function(scr_0hr_df,scr_2hr_df,scr_24hr_df,b3_0hr_df,b3_2hr_df,b3_24hr_df)
  {
    mean_scr_0hr=rowMeans(scr_0hr_df)
    mean_scr_2hr=rowMeans(scr_2hr_df)
    mean_scr_24hr=rowMeans(scr_24hr_df)
    
    mean_b3_0hr=rowMeans(b3_0hr_df)
    mean_b3_2hr=rowMeans(b3_2hr_df)
    mean_b3_24hr=rowMeans(b3_24hr_df)
    
    mean_df=cbind(mean_scr_0hr,mean_scr_2hr,mean_scr_24hr,mean_b3_0hr,mean_b3_2hr,mean_b3_24hr)
    return(mean_df)
  }
  
  get_counts_matrix_from_annotated_counts_df = function(annotated_counts_df)
  {
    selectedColumns=str_detect(colnames(annotated_counts_df),pattern="[.]sam$")
    return(as.matrix(annotated_counts_df[,selectedColumns]))
  }
  
  getGeneCountsDf=function(cross_section_csv_dir,genes)
  {
    x=c(0,2,24)
    files=list.files(path=cross_section_csv_dir,pattern = "[.]csv$")
    tempdf=data.frame()
    for(i in x)
    {
      file=files[str_detect(string = files,pattern = paste(i,"hr",sep=""))]
      print(file)
      df=read.csv(file=getPath(filename = file,directory = cross_section_csv_dir),header = T)
      df=df[df$entrezgene_id!="6349",]
      test=df[as.character(df$hgnc_symbol) %in% genes,]
      selectedColumns=str_detect(columns,pattern="[.]sam$")
      test=test[,selectedColumns]
      colnames(test)=c(paste("scr",1:3,sep="_"),paste("b3",1:3,sep="_"))
      tempdf=rbind(tempdf,test)
    }
    rownames(tempdf)=x
    print(tempdf)
  }
  getGeneWiseCountsReport=function(scr_0hr_df,scr_2hr_df,scr_24hr_df,b3_0hr_df,b3_2hr_df,b3_24hr_df,outputFilePath)
  {
    wb=createWorkbook()
    x=c(0,2,24)
    
    genes=rownames(scr_0hr_df)
    genes=sort(genes)
    i=1
    for(i in 1:length(genes))
    {
      gene=genes[i]
      print(gene)
      df=data.frame()
      sheet=createSheet(wb,sheetName = gene)
      row_scr_0hr=scr_0hr_df[gene,]
      row_scr_2hr=scr_2hr_df[gene,]
      row_scr_24hr=scr_24hr_df[gene,]
      
      row_scr_df=rbind(row_scr_0hr,row_scr_2hr,row_scr_24hr)
      rownames(row_scr_df)=x
      
      row_b3_0hr=b3_0hr_df[gene,]
      row_b3_2hr=b3_2hr_df[gene,]
      row_b3_24hr=b3_24hr_df[gene,]
      
      row_b3_df=rbind(row_b3_0hr,row_b3_2hr,row_b3_24hr)
      rownames(row_b3_df)=x
      
      df=cbind(row_scr_df,row_b3_df)
      addDataFrame(as.data.frame(df),sheet = sheet,row.names = F)
      i=i+1
    }
    saveWorkbook(wb,file = outputFilePath)
  }
  
  getGeneExpressionInfoDf=function(gene,degsDirectoryPath)
  {
    library(dplyr)
    files=list.files(path=degsDirectoryPath,pattern="[.]csv$")
    info_df=data.frame()
    for(file in files)
    {
      df=read.csv(getPath(filename=file,directory =degsDirectoryPath),header=T,stringsAsFactors = F)
      gene_df=df %>% select(gene_symbol,description,log2FoldChange,padj,pvalue) %>% dplyr::filter(gene_symbol==gene)
      if(nrow(gene_df)>0)
      {
        # print("********************************************************************************")
        # print(paste("gene     :",gene))
        # print(paste("found in :",file))
        # print(gene_df)
        # print("********************************************************************************")
        gene_df=gene_df %>% mutate(fileName=file) %>% select(fileName,gene_symbol,description,log2FoldChange,padj,pvalue)
        info_df=rbind(info_df,gene_df)
      }
    }
    return(info_df)
  }
  
  getGeneExpressionExcelReport=function(genes,degsDirectoryPath,outputDirectoryPath,outputFileName)
  {
    library(xlsx)
    wb=createWorkbook()
    print(length(genes))
    i=1
    for(i in 1:length(genes))
    {
      gene=genes[i]
      print("processing gene ")
      print(gene)
      sheet=xlsx::createSheet(wb,sheetName = gene)
      print("adding data frame")
      df=getGeneExpressionInfoDf(gene,degsDirectoryPath = degsDirectoryPath)
      if(nrow(df)>0)
      {
        xlsx::addDataFrame(df,sheet = sheet,row.names = F)
      }
    }
    xlsx::saveWorkbook(wb,file=getPath(filename = outputFileName,directory = outputDirectoryPath))
    gc()
  }
  
  getGeneCountsExcelReport=function(genes,cross_section_csv_dir,outputDirectoryPath,outputFileName,filteredEntrezGeneIds="")
  {
    library(xlsx)
    wb=createWorkbook()
    print(length(genes))
    i=1
    for(i in 1:length(genes))
    {
      gene=genes[i]
      print("processing gene ")
      print(gene)
      sheet=xlsx::createSheet(wb,sheetName = gene)
      print("adding data frame")
      df=getGeneCountsDf(gene = gene,cross_section_csv_dir = cross_section_csv_dir,filteredEntrezGeneIds = filteredEntrezGeneIds)
      if(nrow(df)>0)
      {
        xlsx::addDataFrame(df,sheet = sheet,row.names = T)
      }
    }
    xlsx::saveWorkbook(wb,file=getPath(filename = outputFileName,directory = outputDirectoryPath))
    gc()
  }
  
  getGeneCountsDf=function(cross_section_csv_dir,gene,filteredEntrezGeneIds="")
  {
    library(dplyr)
    x=c(0,2,24)
    files=list.files(path=cross_section_csv_dir,pattern = "[.]csv$")
    tempdf=data.frame()
    for(i in x)
    {
      file=files[str_detect(string = files,pattern = paste(i,"hr",sep=""))]
      print(file)
      df=read.csv(file=getPath(filename = file,directory = cross_section_csv_dir),header = T, stringsAsFactor=F)
      df=df %>% dplyr::filter(!(entrezgene_id %in% filteredEntrezGeneIds))
      test=df %>% dplyr::filter(hgnc_symbol == gene)
      columns=colnames(test)
      selectedColumns=columns[str_detect(columns,pattern="[.]sam$")]
      test=test %>% select(selectedColumns)
      colnames(test)=c(paste("scr",1:3,sep="_"),paste("b3",1:3,sep="_"))
      tempdf=rbind(tempdf,test)
    }
    if(nrow(tempdf)>0){
      rownames(tempdf)=x
    }
    return(tempdf)
  }
  
  ####################################################
  # checkCountsFileForSanity 
  # can be called by providing the file name . be careful the program must be running in the directory where the file is
  ####################################################
  getAnalysisType=function(file)
  {
    counts_detail=str_split(file,pattern = "[.]")
    counts_detail=counts_detail[[1]][1]
    time_points=str_split(string=counts_detail,pattern = "_")[[1]]
    time_points=time_points[str_detect(time_points,pattern="hr")]
    
    base_time_point=time_points[1]
    treatment_time_point=time_points[2]
    
    if(base_time_point==treatment_time_point)
    {
      return("cross_sectional")
    }
    else
    {
      return("longitudinal")
    }
  }
  
  ####################################################
  # checkCountsFileForSanity 
  # can be called by providing the file name . be careful the program must be running in the directory where the file is
  ####################################################
  countsSanityCheck=function(counts_dir)
  {
    wd=getwd()
    setwd(counts_dir)
    
    counts_files=list.files(path=counts_dir, pattern="[.]csv$")
    analysis_type=c()
    qa_status=c()
    for(file in counts_files)
    {
      temp_analysis_type=getAnalysisType(file = file)
      
      print("*********************************************************************************")
      print(file)
      print(temp_analysis_type)
      qa_status=c(qa_status,checkCountsFileForSanity(file=file))
      print("*********************************************************************************")
      analysis_type=c(analysis_type,temp_analysis_type)
    }
    counts_analysis_df=as.data.frame(cbind(counts_files,analysis_type,qa_status),stringsAsFactors=F)
    setwd(wd)
    return(counts_analysis_df)
  }
  
  ####################################################
  # checkCountsFileForSanity 
  # can be called by providing the file name . be careful the program must be running in the directory where the file is
  ####################################################
  checkCountsFileForSanity=function(file="")
  {
    tests_passed=""
    if(file=="")
    {
      print("please provide a file name")
      return()
    }
    else
    {
      
      wd=getwd()
      setwd(counts_dir)
      
      df=read.csv(file = file,header = T,stringsAsFactors = F)
      counts_detail=str_split(file,pattern = "[.]")
      counts_detail=counts_detail[[1]][1]
      time_points=str_split(string=counts_detail,pattern = "_")[[1]]
      time_points=time_points[str_detect(time_points,pattern="hr")]
      
      base_time_point=time_points[1]
      treatment_time_point=time_points[2]
      
      columns=colnames(df)
      columns=tolower(columns)
      
      total_columns=length(columns[2:length(columns)])
      base_time_point_found_in=sum(str_detect(string=columns,pattern = as.character(base_time_point)))
      if(base_time_point_found_in > 0)
      {
        print("base time point check OKAY ")
        tests_passed=paste(tests_passed,"base_time_point_check_okay")
      }else
      {
        print(paste("sample : ",base_time_point,treatment_time_point))
        print("FAILED base time point check ")
        print(base_time_point_found_in)
        tests_passed=paste(tests_passed,"FAILED_base_time_point_check")
      }
      treatment_time_point_found_in=sum(str_detect(pattern = treatment_time_point,string = columns),na.rm = T)
      if(((treatment_time_point_found_in == (total_columns/2)) & (base_time_point != treatment_time_point)) | ((treatment_time_point_found_in == total_columns) & (treatment_time_point == base_time_point)))
      {
        print("treatment_time_point OKAY ")
        tests_passed=paste(tests_passed,"treatment_time_point_check_okay",sep="+")
      }else
      {
        print(paste("sample : ",base_time_point,treatment_time_point))
        if(treatment_time_point=="0hr" & treatment_time_point_found_in >0)
        {
          tests_passed=paste(tests_passed,"treatment_time_point_check_okay",sep="+")
        }else{
          print(treatment_time_point_found_in)
          tests_passed=paste(tests_passed,"FAILED_treatment_time_point_check",sep="+")
        }
      }
      
    }
    setwd(wd)
    return(tests_passed)
  }
  
  get_normalized_counts_df=function(counts_df)
  {
    library(stringr)
    library(DESeq2)
    
    
    print(head(counts_df))
    print(paste("nrow(counts_df[rowSums(counts_df)<10,])",nrow(counts_df[rowSums(counts_df)<10,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<6,])",nrow(counts_df[rowSums(counts_df)<6,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<3,])",nrow(counts_df[rowSums(counts_df)<2,])))
    
    
    counts_df=counts_df[rowSums(counts_df)>2,]
    
    columns=colnames(counts_df)
    
    #preparing the details dataframe
    sam_file_details=as.data.frame(cbind(columns[1:(ncol(counts_df)/2)],columns[((ncol(counts_df)/2)+1): ncol(counts_df)]))
    colnames(sam_file_details)=c("control","treatment")
    
    
    #***********************************************************
    #conducting differential expression analysis
    #***********************************************************
    colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","treat", "treat", "treat")))
    dds = DESeqDataSetFromMatrix(counts_df, colData, formula(~ condition))
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized=TRUE)
    
    return(as.data.frame(normalized_counts))
  }
  
  get_unfiltered_differential_expression_results_df=function(counts_df)
  {
    library(stringr)
    library(DESeq2)
    
    
    
    print(paste("nrow(counts_df[rowSums(counts_df)<10,])",nrow(counts_df[rowSums(counts_df)<10,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<6,])",nrow(counts_df[rowSums(counts_df)<6,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<3,])",nrow(counts_df[rowSums(counts_df)<2,])))
    
    
    counts_df=counts_df[rowSums(counts_df)>2,]
    
    columns=colnames(counts_df)
    
    #preparing the details dataframe
    sam_file_details=as.data.frame(cbind(columns[1:(ncol(counts_df)/2)],columns[((ncol(counts_df)/2)+1): ncol(counts_df)]))
    colnames(sam_file_details)=c("control","treatment")
    
    
    #***********************************************************
    #conducting differential expression analysis
    #***********************************************************
    colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","treat", "treat", "treat")))
    dds = DESeqDataSetFromMatrix(counts_df, colData, formula(~ condition))
    
    # run DEseq
    res= DESeq(dds)
    
    #plotMA(res)
    
    #diffrentiallyExpressed 
    deGenes = results(res,alpha=0.05)
    
    # order by BH adjusted p-value
    deGenesOrdered = deGenes[order(deGenes$padj),]
    print(nrow(deGenesOrdered))
    deGenesOrdered=filterByPvalue(df=deGenesOrdered,cutoff = 0.05)
    print(nrow(deGenesOrdered))
    return(as.data.frame(deGenesOrdered))
  }
  
    # get_differential_expression_results_df=function(counts_df)
    # {
    #   library(stringr)
    #   library(DESeq2)
    #   
    #   
    #   
    #   print(paste("nrow(counts_df[rowSums(counts_df)<10,])",nrow(counts_df[rowSums(counts_df)<10,])))
    #   print(paste("nrow(counts_df[rowSums(counts_df)<6,])",nrow(counts_df[rowSums(counts_df)<6,])))
    #   print(paste("nrow(counts_df[rowSums(counts_df)<3,])",nrow(counts_df[rowSums(counts_df)<2,])))
    #   
    #   
    #   counts_df=counts_df[rowSums(counts_df)>2,]
    #   
    #   columns=colnames(counts_df)
    #   
    #   #preparing the details dataframe
    #   sam_file_details=as.data.frame(cbind(columns[1:(ncol(counts_df)/2)],columns[((ncol(counts_df)/2)+1): ncol(counts_df)]))
    #   colnames(sam_file_details)=c("control","treatment")
    #   
    #   
    #   #***********************************************************
    #   #conducting differential expression analysis
    #   #***********************************************************
    #   colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","treat", "treat", "treat")), replicate = c("1", "2", "3", "1", "2" ,"3"))
    #   
    #   dds = DESeqDataSetFromMatrix(counts_df, colData, formula(~ condition))
    #   
    #   # run DEseq
    #   res= DESeq(dds)
    #   
    #   #plotMA(res)
    #   
    #   #diffrentiallyExpressed 
    #   deGenes = results(res,alpha=0.05)
    #   
    #   # order by BH adjusted p-value
    #   deGenesOrdered = deGenes[order(deGenes$padj),]
    #   #deGenesOrdered=filterByPadjValue(df=deGenesOrdered,cutoff = 0.05)
    #   
    #   return(as.data.frame(deGenesOrdered))
    # }
    # 
  ################################################################################
  # get_annotated_deg_df 
  # returns an annotated differential expression results dataframe
  ################################################################################
  get_annotated_deg_df = function(deg_file_path,degs_df=NA, filter_by_padj = TRUE)
  {
    library(dplyr)
    library(stringr)
    if(nrow(degs_df)==0)
    {
      return(degs_df)
    }
    if(is.na(degs_df))
    {
      if(str_ends(deg_file_path,"csv"))
      {
        degs_df=read.csv(file = deg_file_path,header = T,stringsAsFactors = F)
      }
      if(str_ends(deg_file_path,"xlsx") | str_ends(deg_file_path,"xls"))
      {
        degs_df=read_excel(path =deg_file_path,sheet = 1)
      }
    }
    if(filter_by_padj){
      degs_df=filterByPadjValue(degs_df,cutoff = 0.05)
    }
    print(head(degs_df))
    # columns=colnames(degs_df)
    # columns[1]="ensembl_gene_id"
    # colnames(degs_df)=columns
    # print(head(degs_df))
    #names(degs_df)
    #rownames(degs_df)=degs_df$ensembl_gene_id
    my_ensemble_ids=rownames(degs_df)
    clean_ensemble_ids=getCleanEnsembleIds(my_ensemble_ids)
    degs_df$clean_ensemble_id=clean_ensemble_ids
    degs_df$ensemble_id=my_ensemble_ids
    #length(clean_ensemble_ids)
    #length(unique(clean_ensemble_ids))
    
    #**************************************************
    # fetching annotation from biomart
    #**************************************************
    columns=c("ensembl_gene_id","entrezgene_id","hgnc_symbol","description","gene_biotype")
    annotation_df=getAnnotationFromBiomart(ids = clean_ensemble_ids,filter = "ensembl_gene_id",columns = columns)
    nrow(annotation_df)
    
    
    allowed_biotypes=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
    
    annotation_df=annotation_df %>% dplyr::filter(gene_biotype %in% allowed_biotypes)
    #sum(duplicated(annotation_df$ensembl_gene_id))
    
    
    annotation_summary=as.data.frame(table(annotation_df$ensembl_gene_id))
    #names(annotation_summary)
    #annotation_summary %>% dplyr::filter(Freq>1)
    
    # duplicated_genes=annotation_df[duplicated(annotation_df$ensembl_gene_id),1]
    # annotation_df[annotation_df$ensembl_gene_id %in% duplicated_genes,]
    # degs_df[degs_df$clean_ensemble_id %in% duplicated_genes,]
    
    
    #little experimet
    #degs_df[degs_df$clean_ensemble_id=="ENSG00000143226",] %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id"))
    
    #nrow(degs_df %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id")))
    print(head(degs_df))
    print(head(annotation_df))
    print(class(degs_df))
    print(class(annotation_df))
    degs_df=degs_df %>% dplyr::inner_join(annotation_df,by=c("clean_ensemble_id"="ensembl_gene_id"))
    names(degs_df)
    print(head(degs_df))
    degs_df=degs_df %>% dplyr::select(hgnc_symbol,description,log2FoldChange,padj,pvalue,ensemble_id,clean_ensemble_id,entrezgene_id) %>% dplyr::rename(gene_symbol=hgnc_symbol) %>% arrange(padj)
    
    return(degs_df)
  }
  
  
  generate_counts_data_excel=function(human_readable_counts_df,output_filepath)
  {
    library(openxlsx)
    head(human_readable_counts_df)
    
    all_data=human_readable_counts_df
    
    data_0hr_df= all_data %>% dplyr::filter(time==0)
    data_2hr_df= all_data %>% dplyr::filter(time==2)
    data_24hr_df=  all_data %>% dplyr::filter(time==24)
    
    wb=createWorkbook()
    addWorksheet(wb,sheetName = "all_data")
    addWorksheet(wb,sheetName = "0hr")
    addWorksheet(wb,sheetName = "2hr")
    addWorksheet(wb,sheetName = "24hr")
    
    writeData(wb=wb,sheet="all_data",x=all_data)
    writeData(wb=wb,sheet="0hr",x=data_0hr_df)
    writeData(wb=wb,sheet="2hr",x=data_2hr_df)
    writeData(wb=wb,sheet="24hr",x=data_24hr_df)
    
    saveWorkbook(wb,file = output_filepath, overwrite = TRUE)
    detach("package:openxlsx", unload = TRUE)
  }
  
  
  ###########################################################
  #get_count_matrix
  #this function return a count matrix for the raw counts df supplied to it
  ###########################################################
  
  get_count_matrix = function(counts_df)
  {
    if(nrow(counts_df)!=0)
    {
      rownames(counts_df) = counts_df[,1]
      counts_df = counts_df[,-1]
    }
    return(counts_df)
  }
  
  ###########################################################
  #get_count_matrix
  #this function return a count matrix for the raw counts df supplied to it
  ###########################################################
  
  get_count_matrix_from_csv = function(counts_file_path)
  {
    counts_df = read.csv(file = counts_file_path, stringsAsFactors = FALSE)
    if(nrow(counts_df)!=0)
    {
      rownames(counts_df) = counts_df[,1]
      counts_df = counts_df[,-1]
    }
    return(counts_df)
  }
  
  #############################################
  # modify_count_matrix_and_save_as_csv
  # takes an existing csv file as input and returns a modified csv file
  #############################################
  modify_count_matrix_and_save_as_csv=function(input_file_path="",output_file_path="",to_modify=0)
  {
    df=read.csv(file =input_file_path,header = T,stringsAsFactors = F)
    df[,3:8]=get_modified_counts_df(df[,3:8],to_modify=to_modify)
    write.csv(df,file = output_file_path,row.names = F)
    return(df)
  }
  #############################################
  # do_venn_analysis
  # does the venn analysis by using two standard differential expression result files
  #############################################
  do_venn_analysis=function(differential_expr_csv_filepath1,differential_expr_csv_filepath2,set1_name,set2_name,excel_report_output_dir)
  {
    library(VennDiagram)
    library(xlsx)
    
    quick_summary = get_quick_summary(degs1_filepath = differential_expr_csv_filepath1,degs2_filepath = differential_expr_csv_filepath2,set1_name = set1_name, set2_name = set2_name)
    #setting the venn filepath
    venn_file_name=paste("venn_analysis_modified_",set1_name,set2_name,".png",sep="")
    venn_diagram_file_path=getPath(filename = venn_file_name,directory = excel_report_output_dir)
    
    #setting the title of  the venn diagram
    main_title=paste("Venn (modified counts)",set1_name,set2_name)
    
    #loading the required dataframes
    set1_degs_df=read.csv(file=differential_expr_csv_filepath1,header = T,stringsAsFactors = F)
    print(head(set1_degs_df))
    set2_degs_df=read.csv(file=differential_expr_csv_filepath2,header = T,stringsAsFactors = F)
    print(head(set2_degs_df))
    #generating the venn diagram
    set_list=list(set1_degs_df$gene_symbol,set2_degs_df$gene_symbol)
    names(set_list)=c(set1_name,set2_name)
    venn.diagram(set_list,filename = venn_diagram_file_path,fill=c("blue","red"),imagetype = "png",scaled=F,cat.pos=c(0,0))
    
    #generating the details dataframe
    details_df=getDetailsDf(set1=set1_degs_df$gene_symbol,set2=set2_degs_df$gene_symbol,names = c(set1_name,set2_name))
    
    #generating the commons dataframe
    common_df=set1_degs_df %>% inner_join(set2_degs_df,c("gene_symbol"="gene_symbol"),suffix=paste("_",c(set1_name,set2_name),sep="")) 
    print(head(common_df))
    only_in_set1_df= set1_degs_df %>% dplyr::filter(!(gene_symbol %in% common_df$gene_symbol))
    print(head(only_in_set1_df))
    only_in_set2_df= set2_degs_df %>% dplyr::filter(!(gene_symbol %in% common_df$gene_symbol))
    print(head(only_in_set2_df))
    #setting the excel workbook name
    excel_file_name=paste("excelReport_venn_",set1_name,"_",set2_name,".xlsx",sep = "")
    excel_file_path=getPath(filename = excel_file_name,dir=excel_report_output_dir)
    #generating the workbook
    wb=createWorkbook()
    
    #add details sheets
    sheet=createSheet(wb,sheetName = "details")
    addDataFrame(x=as.data.frame(details_df),row.names = F,sheet = sheet)
    addDataFrame(x=as.data.frame(quick_summary),row.names = T,sheet = sheet, startRow = nrow(details_df) + 3)
    addPicture(file=venn_diagram_file_path,sheet=sheet,startRow = 4,startColumn = 6)
    
    #add commons sheet
    sheet=createSheet(sheetName = "common",wb=wb)
    addDataFrame(x=as.data.frame(common_df),sheet = sheet,row.names = F)
    
    #********************************
    #picking up gene names from set1. the log fold change will also be from set 1. more can be done in this to incorporate the set 2 genes as well 
    #********************************
    temp_common_df=set1_degs_df %>% dplyr::filter(gene_symbol %in% common_df$gene_symbol)
    common_genelist_for_gsea=get_genelist_for_gsea(temp_common_df)
    try(doGSEA(geneList = common_genelist_for_gsea,wb=wb,sheetNamePrefix = "pthwys_common_genes_"),silent=T)
    
    #add onlyInSet1 sheet
    sheet=createSheet(sheetName = paste("onlyIn_",set1_name,sep=""),wb=wb)
    addDataFrame(x=as.data.frame(only_in_set1_df),sheet = sheet,row.names = F)
    onlyinSet1_genelist_for_gsea=get_genelist_for_gsea(only_in_set1_df)
    try(doGSEA(geneList = onlyinSet1_genelist_for_gsea,wb=wb,sheetNamePrefix = paste("pthwys_onlyIn_",set1_name)),silent=T)
    
    #add onlyInSet2 sheet
    sheet=createSheet(sheetName = paste("onlyIn_",set2_name,sep=""),wb=wb)
    addDataFrame(x=as.data.frame(only_in_set2_df),sheet = sheet,row.names = F)
    onlyinSet2_genelist_for_gsea=get_genelist_for_gsea(only_in_set2_df)
    try(doGSEA(geneList = onlyinSet2_genelist_for_gsea,wb=wb,sheetNamePrefix = paste("pthwys_onlyIn_",set2_name)),silent=T)
    
    saveWorkbook(wb,file=excel_file_path)
  }
  
  
  #############################################
  # get_genelist_for_gsea
  #this function return a geneList which is used by clusterprofiler for gsea
  #############################################
  get_genelist_for_gsea=function(my_degs_df)
  {
    library(dplyr)
    my_degs_df=my_degs_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    #print(head(my_degs_df))
    genelist=my_degs_df$log2FoldChange
    #print(head(my_degs_df$gene_symbol))
    names(genelist)=my_degs_df$entrezgene_id
    genelist=sort(genelist,decreasing = T)
    return(genelist)
  }
  #############################################
  # addRowToDetailsDf
  # adds a row to out details dataframe (2 columns, detail and value)
  #############################################
  addRowToDetailsDf=function(df,detail,value)
  {
    rownames(df)=NULL
    df=rbind(df,cbind(detail,value))
    return(df)
  }
  
  #############################################
  # getDetailsDf
  # by taking two sets of gene lists ( gene symbol vectors) returns a dataframe with info on overlap (2 columns, detail and value)
  #############################################
  getDetailsDf=function(set1,set2,names=c())
  {
    common=set1[set1 %in% set2]
    total_common=length(common)
    
    total_set1=length(set1)
    total_set2=length(set2)
    
    only_set1=set1[!(set1 %in% common)]
    total_only_set1=length(only_set1)
    only_set2=set2[!(set2 %in% common)]
    total_only_set2=length(only_set2)
    
    overlap_perc_set1=format(round((total_common/total_set1)*100,2),nsmall = 2)
    overlap_perc_set2=format(round((total_common/total_set2)*100,2),nsmall = 2)
    
    df=data.frame(stringsAsFactors = F)
    df=addRowToDetailsDf(df,detail = paste("total in",names[1]),value = total_set1)
    df=addRowToDetailsDf(df,detail = paste("total in",names[2]),value = total_set2)
    
    df=addRowToDetailsDf(df,detail = paste("common"),value = total_common)
    df=addRowToDetailsDf(df,detail = paste("only in",names[1]),value = total_only_set1)
    df=addRowToDetailsDf(df,detail = paste("only in",names[2]),value = total_only_set2)
    
    df=addRowToDetailsDf(df,detail = paste("overlap of",names[1]),value = paste(overlap_perc_set1,"%"))
    df=addRowToDetailsDf(df,detail = paste("overlap of",names[2]),value = paste(overlap_perc_set2,"%"))
    return(df)
  }
  
  
  #############################################
  # get_quick_summary
  # the get_quick_summary function takes in two degs files and return valuable summary from the two files which can be directly used ( and is used) in presentations
  #############################################
  get_quick_summary = function(degs1_filepath, degs2_filepath, set1_name,set2_name, filter_by_padj = TRUE)
  {
    library(dplyr)
    set1_df = read.csv(file = degs1_filepath, stringsAsFactors = FALSE) 
    set2_df = read.csv(file = degs2_filepath, stringsAsFactors = FALSE)
    set1_df = set1_df%>% dplyr::filter(!is.na(padj))
    set2_df = set2_df%>% dplyr::filter(!is.na(padj))
    
    set1_df = set1_df %>% dplyr::filter(!is.na(padj), gene_symbol!="")
    set2_df = set2_df %>% dplyr::filter(!is.na(padj), gene_symbol!="")
    set1_df = set1_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    set2_df = set2_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    
    if(filter_by_padj)
    {
      set1_df = set1_df%>% dplyr::filter(padj<0.05)
      set2_df = set2_df%>% dplyr::filter(padj<0.05)
    }
    
    set1_total_genes = nrow(set1_df)
    set2_total_genes = nrow(set2_df)
    
    col_names = c(set1_name,set2_name)
    row_names = c("total","up (>= 2 fold change)", "down (>= 2 fold change)", paste("up (exclusive)"), paste("down (exclusive)"))
    df=data.frame(row.names = row_names, stringsAsFactors = F)
    
    print(rownames(df))
    set1_up_genes = sum(set1_df$log2FoldChange >= 1)
    set2_up_genes = sum(set2_df$log2FoldChange >= 1)
    
    set1_down_genes = sum(set1_df$log2FoldChange <= -1)
    set2_down_genes = sum(set2_df$log2FoldChange <= -1)
    
    set1_up_gene_symbols = set1_df$gene_symbol[set1_df$log2FoldChange >= 1]
    set2_up_gene_symbols = set2_df$gene_symbol[set2_df$log2FoldChange >= 1]
    
    print(length(set1_up_gene_symbols))
    print(length(set2_up_gene_symbols))
    
    set1_down_gene_symbols = set1_df$gene_symbol[set1_df$log2FoldChange <= -1]
    set2_down_gene_symbols = set2_df$gene_symbol[set2_df$log2FoldChange <= -1]
    
    print(length(set1_down_gene_symbols))
    print(length(set2_down_gene_symbols))
    
    only_set1_up_genes = sum(!(set1_up_gene_symbols %in% set2_up_gene_symbols))
    only_set2_up_genes = sum(!(set2_up_gene_symbols %in% set1_up_gene_symbols))
    
    only_set1_down_genes = sum(!(set1_down_gene_symbols %in% set2_down_gene_symbols))
    only_set2_down_genes = sum(!(set2_down_gene_symbols %in% set1_down_gene_symbols))
    
    df[,set1_name] = c(set1_total_genes,set1_up_genes,set1_down_genes,only_set1_up_genes,only_set1_down_genes)
    df[,set2_name] = c(set2_total_genes,set2_up_genes,set2_down_genes,only_set2_up_genes,only_set2_down_genes)
    
    # set_list=list(set1_up_gene_symbols,set2_up_gene_symbols)
    # names(set_list)=c(set1_name,set2_name)
    # 
    # up_venn_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/plots/up_test_venn.png"
    # venn.diagram(set_list,filename = up_venn_filepath,fill=c("blue","red"),imagetype = "png",scaled=F,cat.pos=c(0,0))
    # 
    # set_list=list(set1_down_gene_symbols,set2_down_gene_symbols)
    # names(set_list)=c(set1_name,set2_name)
    # 
    # down_venn_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/plots/down_test_venn.png"
    # venn.diagram(set_list,filename = down_venn_filepath,fill=c("blue","red"),imagetype = "png",scaled=F,cat.pos=c(0,0))
    
    return(df)
  }
  ####################################################
  # myZeroHourOutlierDetector 
  # seraches for outlier in the upper half of the sample vector
  # expects a sample size of 6 (3 n each group)
  ####################################################
  myZeroHourOutlierDetector=function(x)
  {
    b3=x[1:3]
    #print(b3)
    max_b3=max(b3)
    #print(max_b3)
    max_index=which(b3==max_b3)[1]
    #print(max_index)
    my_mean_b3=sum(b3[-max_index],na.rm = T)
    #print(my_mean_b3)
    if(my_mean_b3==0)
    {
      my_mean_b3=1
    }
    if(max_b3>(2*my_mean_b3))
    {
      # if(max_b3>100)
      # {
      #   return(T)
      # }else{
      #   return(F)
      # }
      return(T)
    }
    return(F)
  }
  
  ####################################################
  # myTwoHourOutlierDetector 
  # seraches for outlier in the lower half of the sample vector
  # expects a sample size of 6 (3 n each group)
  ####################################################
  myTwoHourOutlierDetector=function(x)
  {
    b3=x[4:6]
    #print(b3)
    max_b3=max(b3)
    #print(max_b3)
    max_index=which(b3==max_b3)[1]
    #print(max_index)
    my_mean_b3=sum(b3[-max_index],na.rm = T)
    #print(my_mean_b3)
    if(my_mean_b3==0)
    {
      my_mean_b3=1
    }
    if(max_b3>(111*my_mean_b3))
    {
      # if(max_b3>100)
      # {
      #   return(T)
      # }else{
      #   return(F)
      # }
      return(T)
    }
    return(F)
  }
  
  ####################################################
  # myRowMean
  # calculats the mean after excluding the outlier and returns the same
  ####################################################
  myRowMean=function(x)
  {
    b3=x[1:3]
    #print(b3)
    max_b3=max(b3)
    #print(max_b3)
    max_index=which(b3==max_b3)[1]
    #print(max_index)
    my_mean_b3=sum(b3[-max_index],na.rm = T)/length(b3[-max_index])
    #print(my_mean_b3)
    return(round(my_mean_b3))
  }
  
  ####################################################
  # getOutlierIndex
  # returns the index of the outlier
  ####################################################
  getOutlierIndex=function(x,to_modify=0)
  {
    #print(paste("in get outlier index : to modify - ",to_modify))
    if(to_modify==0)
    {
      b3=x[1:3]
    }
    if(to_modify==2)
    {
      b3=x[4:6]
    }
    #print(b3)
    max_b3=max(b3)
    #print(max_b3)
    #cat("B3 = ",b3,"\n")
    #cat("max b3 ",max_b3,"\n")
    #cat("x = ",x,"\n")
    #cat("which(x==max_b3) = ",which(x==max_b3),"\n")
    max_index=which(x==max_b3)[1]
    #print(max_index)
    #print(names(max_index))
    return(max_index)
  }
  
  ####################################################
  # getOutlierName
  # returns the name of the outlier sample ( columname)
  ####################################################
  getOutlierName=function(x)
  {
    b3=x[1:3]
    #print(b3)
    max_b3=max(b3)
    #print(max_b3)
    max_index=which(b3==max_b3)[1]
    #print(names(max_index))
    return(names(max_index))
  }
  
  ###################################################
  # differentialExpressionAnalysis
  # writes the differential expression excel report
  ###################################################
  differentialExpressionAnalysis=function(counts_dir,counts_filename,diffential_expression_output_dir,diffential_expression_output_filename)
  {
    library(xlsx)
    library(stringr)
    library(DESeq2)
    wb=createWorkbook()
    
    wd=getwd()
    setwd(counts_dir)
    #***********************************************************
    #reading raw reads from counts file
    #***********************************************************
    countsFilePath=getPath(filename = counts_filename,directory = counts_dir)
    
    df=read.csv(file=countsFilePath,stringsAsFactors = F,header=T)
    rownames(df)=df[,1]
    df=df[,2:ncol(df)]
    
    # rMedian=apply(df,1,myRowMedian)
    # print(paste("sum(rMedian==0)",sum(rMedian==0)))
    print(paste("nrow(df[rowSums(df)<10,])",nrow(df[rowSums(df)<10,])))
    print(paste("nrow(df[rowSums(df)<6,])",nrow(df[rowSums(df)<6,])))
    print(paste("nrow(df[rowSums(df)<3,])",nrow(df[rowSums(df)<2,])))
    df=df[rowSums(df)>2,]
    
    #print(paste("counts data : ",countsFilePath))
    #write.csv(df,file = countsFilePath,row.names = T,col.names = T)
    columns=colnames(df)
    
    sam_file_details=as.data.frame(cbind(columns[1:(ncol(df)/2)],columns[((ncol(df)/2)+1): ncol(df)]))
    colnames(sam_file_details)=c("control","treatment")
    #print(sam_file_details)
    
    
    
    #***********************************************************
    #conducting differential expression analysis
    #***********************************************************
    colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","treat", "treat", "treat")))
    dds = DESeqDataSetFromMatrix(df, colData, formula(~ condition))
    
    # run DEseq
    res= DESeq(dds)
    
    #plotMA(res)
    
    #diffrentiallyExpressed 
    deGenes <- results(res,alpha=0.05)
    
    # order by BH adjusted p-value
    deGenesOrdered <- deGenes[order(deGenes$padj),]
    #deGenesOrdered=deGenesOrdered[!is.na(deGenesOrdered$padj),]
    #deGenesOrdered <- deGenesOrdered[deGenesOrdered$padj<0.05,]
    
    
    #######################################
    #setting up file paths for differential expression
    #######################################
    differential_expression_report_path=getPath(filename = diffential_expression_output_filename,directory = diffential_expression_output_dir)
    
    
    #print(paste("diffex csv path : ",differentialExpressionResultsFilePath))
    print(paste("diffex excel path : ",differential_expression_report_path))
    
    
    #######################################
    #adding results to excel sheets
    #######################################
    sheet=createSheet(wb,sheetName = "diff_ex_result")
    addDataFrame(deGenesOrdered,sheet = sheet,row.names = T)
    
    sheet=createSheet(wb,sheetName = "sam_file_details")
    addDataFrame(sam_file_details,sheet = sheet,row.names = F)
    
    #######################################
    #saving files
    #######################################
    #write.csv(deGenesOrdered,file=differentialExpressionResultsFilePath,row.names = T,col.names = T)
    saveWorkbook(wb,file=differential_expression_report_path)
    setwd(wd)
    # top of ordered matrix
    #head(deGenesOrdered)
    gc()
  }
  
  #############################################################################
  #get_mega_report
  # this function generates a master report for 2 differential expression files 
  # this is the most complete report one can get for the analysis. 
  # this is more like a venn analysis but detailed
  #############################################################################
  get_mega_report = function(degs1_filepath, degs2_filepath, set1_name, set2_name, excel_report_output_dir, pathwaysType = c("gsea","enricher"))
  {
    
    #loading the required packages for making an excel report and generating venn diagrams
    library(VennDiagram)
    library(xlsx)
    library(dplyr)
    pathwaysType = match.arg(pathwaysType)
    #reading the differential expression files into dataframes
    #the expression files should be in the standard format 
    degs1_df = read.csv(file = degs1_filepath, stringsAsFactors = FALSE)
    degs2_df = read.csv(file = degs2_filepath, stringsAsFactors = FALSE)
    degs1_df = degs1_df %>% dplyr::filter(!is.na(padj), gene_symbol!="")
    degs2_df = degs2_df %>% dplyr::filter(!is.na(padj), gene_symbol!="")
    degs1_df = degs1_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    degs2_df = degs2_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    degs1_df = as.data.frame(degs1_df)
    degs2_df = as.data.frame(degs2_df)
    #reordering the differential expression data based on the log2FoldChange , ready for analysis
    degs1_df=degs1_df[order(as.numeric(degs1_df$log2FoldChange),decreasing = TRUE),]
    degs2_df=degs2_df[order(as.numeric(degs2_df$log2FoldChange),decreasing = TRUE),]
    
    #filtering based padj value cutoff 
    padj_cutoff = 0.05
    degs1_df = degs1_df[degs1_df$padj < padj_cutoff,]
    degs2_df = degs2_df[degs2_df$padj < padj_cutoff,]
    
    #getting a quick summary of the degs
    #this quick summary will be used in the reoprt
    #the quick summary contains the upregulated and downregulated gene numbers
    print("fetching quick summary for mega report")
    quick_summary = get_quick_summary(degs1_filepath = degs1_filepath, degs2_filepath = degs2_filepath, set1_name = set1_name, set2_name = set2_name)
    
    #generating the details dataframe
    #this also contains the overlap proportions
    print("fetching deteails like %overlap between the two sets")
    details_df=getDetailsDf(set1=degs1_df$gene_symbol,set2=degs2_df$gene_symbol,names = c(set1_name,set2_name))
    
    #genelists are created for doing gene set enrichment analysis
    print("fetching the genelist for set 1 from degs1_df")
    all_set1_genelist = get_genelist_for_gsea(degs1_df)
    print("fetching the genelist for set 2 from degs2_df")
    all_set2_genelist = get_genelist_for_gsea(degs2_df)
    
    #this filter can be adjusted based on the idea of important genes to consider
    #here 1 and -1 represent 2 and 1/2 fold changes 
    upregulated_filter = 1 #2 fold upregulated
    downregulated_filter = -1 #2 fold downregulation filter
    
    #upregulated genes from set1
    set1_upregulated_degs_df = degs1_df[degs1_df$log2FoldChange >= upregulated_filter,]
    #upregulated genes from set2
    set2_upregulated_degs_df = degs2_df[degs2_df$log2FoldChange >= upregulated_filter,]
    
    #downregulated genes from set1
    set1_downregulated_degs_df = degs1_df[degs1_df$log2FoldChange <= downregulated_filter,]
    set1_downregulated_degs_df = set1_downregulated_degs_df[order(as.numeric(set1_downregulated_degs_df$log2FoldChange), decreasing = FALSE),]
    #downregulated genes from set2
    set2_downregulated_degs_df = degs2_df[degs2_df$log2FoldChange <= downregulated_filter,]
    set2_downregulated_degs_df = set2_downregulated_degs_df[order(as.numeric(set2_downregulated_degs_df$log2FoldChange), decreasing = FALSE),]
    #********************************
    #info regarding all degs
    #********************************
    #retrieving the common gene set from the complete data
    print("fetchign common degs")
    common_degs_df = compare_and_get_degs(degs1_df = degs1_df, degs2_df = degs2_df, onlyIn = "common")
    #retrieving the genes only in set1 from the entire set
    print("fetching onlyIn_set1_degs_df")
    onlyIn_set1_degs_df = compare_and_get_degs(degs1_df = degs1_df, degs2_df = degs2_df, onlyIn = "set1")
    #retreiving the genes only in set2 from the entire set
    print("fetching onlyIn_set2_degs_df")
    onlyIn_set2_degs_df = compare_and_get_degs(degs1_df = degs1_df, degs2_df = degs2_df, onlyIn = "set2")
    
    #genelists are created for unique genes doing gene set enrichment analysis
    print("fetching onlyIn_set1_genelist")
    onlyIn_set1_genelist = get_genelist_for_gsea(onlyIn_set1_degs_df)
    print("fetching onlyIn_set2_genelist")
    onlyIn_set2_genelist = get_genelist_for_gsea(onlyIn_set2_degs_df)
    
    #********************************
    #info regarding upregulated degs
    #********************************
    #retrieving common upregulated genes
    print("fetching common_up_degs_df")
    common_up_degs_df = compare_and_get_degs(degs1_df = set1_upregulated_degs_df, degs2_df = set2_upregulated_degs_df, onlyIn = "common")
    #retreiving upregulated genes in set1 only
    print("fetching onlyIn_set1_upregulated_degs_df")
    onlyIn_set1_upregulated_degs_df = compare_and_get_degs(degs1_df = set1_upregulated_degs_df, degs2_df = set2_upregulated_degs_df, onlyIn = "set1")
    #retrieving upregulated genes in set2 only
    print("fetching onlyIn_set2_upregulated_degs_df")
    onlyIn_set2_upregulated_degs_df = compare_and_get_degs(degs1_df = set1_upregulated_degs_df, degs2_df = set2_upregulated_degs_df, onlyIn = "set2")
    
    #********************************
    #info regarding downregulated degs
    #********************************
    #getting common downregulated genes
    print("fetching common_down_degs_df")
    common_down_degs_df = compare_and_get_degs(degs1_df = set1_downregulated_degs_df, degs2_df = set2_downregulated_degs_df, onlyIn = "common")
    #getting doownregulated genes only in set1
    print("fetching onlyIn_set1_downregulated_degs_df")
    onlyIn_set1_downregulated_degs_df = compare_and_get_degs(degs1_df = set1_downregulated_degs_df, degs2_df = set2_downregulated_degs_df, onlyIn = "set1")
    #getting downregulated genes only in set2
    print("fetching onlyIn_set2_downregulated_degs_df")
    onlyIn_set2_downregulated_degs_df = compare_and_get_degs(degs1_df = set1_downregulated_degs_df, degs2_df = set2_downregulated_degs_df, onlyIn = "set2")
    
    #getting the gene list for gsea of set1 upregulated genes
    print("fetching up_set1_genelist")
    up_set1_genelist = get_genelist_for_gsea(set1_upregulated_degs_df)
    #getting the gene list for gsea of set2 upregulated genes
    print("fetching up_set2_genelist")
    up_set2_genelist = get_genelist_for_gsea(set2_upregulated_degs_df)
    
    #getting the gene list for gsea of set1 downregulated genes
    print("fetching down_set1_genelist")
    down_set1_genelist = get_genelist_for_gsea(set1_downregulated_degs_df)
    #getting the gene list for gsea of set2 downregulated genes
    print("fetching down_set2_genelist")
    down_set2_genelist = get_genelist_for_gsea(set2_downregulated_degs_df)
    
    #*****************************************************
    #getting pathways from whole gene set
    #*****************************************************
    print("fetching all_set1_pathways_df")
    all_set1_pathways_df = try(getKeggPathways(all_set1_genelist, pathwaysType = "gsea"), silent = TRUE)
    print("fetching all_set2_pathways_df")
    all_set2_pathways_df = try(getKeggPathways(all_set2_genelist, pathwaysType = "gsea"), silent = TRUE)
    
    print("fetching set1_only_pathways_df")
    #set1_only_pathways_df = getKeggPathwaysOnlyIn(set1_df = all_set1_pathways_df, set2_df = all_set2_pathways_df,onlyIn = "set1")
    set1_only_pathways_df = getKeggPathwaysOnlyIn(set1_df = all_set1_pathways_df, set2_df = all_set2_pathways_df,onlyIn = "set1")
    #print(head(set1_only_pathways_df))
    print("fetching set2_only_pathways_df")
    set2_only_pathways_df = getKeggPathwaysOnlyIn(set1_df = all_set1_pathways_df, set2_df = all_set2_pathways_df,onlyIn = "set2")
    #print(head(set2_only_pathways_df))
    print("fetching common_pathways_df")
    common_pathways_df = getKeggPathwaysOnlyIn(set1_df = all_set1_pathways_df, set2_df = all_set2_pathways_df, onlyIn = "common")
    if(nrow(common_pathways_df)==0)
    {
      common_pathways_df   = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set1_only_pathways_df) == 0)
    {
      set1_only_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set2_only_pathways_df) == 0)
    {
      set2_only_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(all_set1_pathways_df)))
    {
      all_set1_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(all_set2_pathways_df)))
    {
      all_set2_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    #*****************************************************
    #getting unique genes pathways from whole gene set
    #*****************************************************
    print("fetching pathways_from_unique_genes_in_set1_df")
    pathways_from_unique_genes_in_set1_df = try(getKeggPathways(onlyIn_set1_genelist, pathwaysType = "enricher"), silent = TRUE)
    print("fetching pathways_from_unique_genes_in_set2_df")
    pathways_from_unique_genes_in_set2_df = try(getKeggPathways(onlyIn_set2_genelist, pathwaysType = "enricher"), silent = TRUE) 
    
    if(is.null(nrow(pathways_from_unique_genes_in_set1_df)))
    {
      pathways_from_unique_genes_in_set1_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(pathways_from_unique_genes_in_set2_df)))
    {
      pathways_from_unique_genes_in_set2_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    #*************************************
    #collecting pathways from upregulated genes in both the datasets
    #*************************************
    print("fetching up_set1_pathways_df")
    up_set1_pathways_df = try(getKeggPathways(up_set1_genelist, pathwaysType = "enricher"), silent = TRUE)
    print("fetching up_set2_pathways_df")
    up_set2_pathways_df = try(getKeggPathways(up_set2_genelist, pathwaysType = "enricher"), silent = TRUE)
    
    print("fetching common_up_pathways_df")
    common_up_pathways_df = getKeggPathwaysOnlyIn(set1_df = up_set1_pathways_df, set2_df = up_set2_pathways_df,onlyIn = "common")
    print("fetching set1_only_up_pathways_df")
    set1_only_up_pathways_df = getKeggPathwaysOnlyIn(set1_df = up_set1_pathways_df, set2_df = up_set2_pathways_df,onlyIn = "set1")
    print("fetching set2_only_up_pathways_df")
    set2_only_up_pathways_df = getKeggPathwaysOnlyIn(set1_df = up_set1_pathways_df, set2_df = up_set2_pathways_df,onlyIn = "set2")
    
    if(nrow(common_up_pathways_df)==0)
    {
      common_up_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set1_only_up_pathways_df)==0)
    {
      set1_only_up_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set2_only_up_pathways_df)==0)
    {
      set2_only_up_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(up_set1_pathways_df)))
    {
      up_set1_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(up_set2_pathways_df)))
    {
      up_set2_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    #*************************************
    #collecting pathways from downregulated genes in both the datasets
    #*************************************
    print("fetching down_set1_pathways_df")
    down_set1_pathways_df = try(getKeggPathways(down_set1_genelist, pathwaysType = "enricher"), silent = TRUE)
    print("fetching down_set2_pathways_df")
    down_set2_pathways_df = try(getKeggPathways(down_set2_genelist, pathwaysType = "enricher"), silent = TRUE)
    
    print("fetching set1_only_down_pathways_df")
    common_down_pathways_df = getKeggPathwaysOnlyIn(set1_df = down_set1_pathways_df, set2_df = down_set2_pathways_df,onlyIn = "common")
    print("fetching set1_only_down_pathways_df")
    set1_only_down_pathways_df = getKeggPathwaysOnlyIn(set1_df = down_set1_pathways_df, set2_df = down_set2_pathways_df,onlyIn = "set1")
    print("fetching set2_only_down_pathways_df")
    set2_only_down_pathways_df = getKeggPathwaysOnlyIn(set1_df = down_set1_pathways_df, set2_df = down_set2_pathways_df,onlyIn = "set2")
    if(nrow(common_down_pathways_df)==0)
    {
      common_down_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set1_only_down_pathways_df)==0)
    {
      set1_only_down_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(nrow(set2_only_down_pathways_df)==0)
    {
      set2_only_down_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    
    if(is.null(nrow(down_set1_pathways_df)))
    {
      down_set1_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    if(is.null(nrow(down_set2_pathways_df)))
    {
      down_set2_pathways_df = data.frame(c("no pathways enriched"), stringsAsFactors = FALSE)
    }
    #*************************************
    #creating venn diagram
    #*************************************
    venn_file_name=paste("venn_analysis_modified_",set1_name,set2_name,".png",sep="")
    venn_diagram_file_path=getPath(filename = venn_file_name,directory = excel_report_output_dir)
    #generating the venn diagram
    set_list=list(degs1_df$gene_symbol, degs2_df$gene_symbol)
    names(set_list)=c(set1_name, set2_name)
    print("generating venn diagram")
    venn.diagram(set_list,filename = venn_diagram_file_path,fill=c("blue","red"),imagetype = "png",scaled=F,cat.pos=c(0,0))
    #*************************************
    #creating and adding data to workbook
    #*************************************
    wb=createWorkbook()
    #*************************************
    #adding details and venn diagram to the workbook
    #*************************************
    #add details sheets
    print("creating details sheet")
    sheet=createSheet(wb,sheetName = "details")
    print("adding details_df")
    addDataFrame(x=as.data.frame(details_df),row.names = FALSE,sheet = sheet)
    print("adding quick_summary")
    addDataFrame(x=as.data.frame(quick_summary),row.names = T,sheet = sheet, startRow = nrow(details_df) + 3)
    print("adding picture")
    addPicture(file=venn_diagram_file_path,sheet=sheet,startRow = 4,startColumn = 6)
    
    #*************************************
    #adding existing degs knowledge data to workbook
    #*************************************
    TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
    
    sheet = createSheet(wb,sheetName=paste(set1_name,"degs", sep = "_"))
    addDataFrame(colnamesStyle = TABLE_COLNAMES_STYLE, x=degs1_df,sheet=sheet,row.names = FALSE)
    setColumnWidth(sheet = sheet,colIndex = 2, colWidth = 50)
    setColumnWidth(sheet = sheet,colIndex = 1, colWidth = 15)
    sheet = createSheet(wb,sheetName=paste("pthwys",set1_name, sep = "_"))
    print("all_set1_pathways_df")
    print(all_set1_pathways_df)
    addDataFrame(x=all_set1_pathways_df,sheet=sheet,row.names = FALSE)
    print("added all_set1_pathways_df")
    
    
    sheet = createSheet(wb,sheetName=paste(set2_name,"degs", sep = "_"))
    print("adding set2_df degs")
    addDataFrame(x=degs2_df,sheet=sheet,row.names = FALSE)
    
    sheet = createSheet(wb,sheetName=paste("pthwys",set2_name, sep = "_"))
    print("adding all_set1_pathways_df")
    addDataFrame(x=all_set2_pathways_df,sheet=sheet,row.names = FALSE)
    print("all_set2_pathways_df")
    #*************************************
    #adding normal data to workbook
    #*************************************
    #adding common genes from set1 and set2
    sheet = createSheet(wb,sheetName="common_degs")
    addDataFrame(x=common_degs_df,sheet=sheet,row.names = FALSE)
    
    #adding unique genes from set1 
    sheet = createSheet(wb,sheetName=paste("only",set1_name,"degs",sep = "_"))
    addDataFrame(x = onlyIn_set1_degs_df, sheet=sheet,row.names = F)
    
    #*************************************
    #adding pathways which are obtained from unique genesof set1
    #*************************************
    sheet = createSheet(wb,sheetName=paste("uniq_genes_pthwys",set1_name,sep = "_"))
    addDataFrame(x=pathways_from_unique_genes_in_set1_df,sheet=sheet,row.names = F)
    print("pathways_from_unique_genes_in_set1_df")
    #adding unique genes from set2
    sheet = createSheet(wb,sheetName=paste("only",set2_name,"degs",sep = "_"))
    addDataFrame(x = onlyIn_set2_degs_df, sheet=sheet,row.names = F)
    #*************************************
    #adding pathways which are obtained from unique genesof set2
    #*************************************
    sheet = createSheet(wb,sheetName=paste("uniq_genes_pthwys",set2_name,sep = "_"))
    addDataFrame(x=pathways_from_unique_genes_in_set2_df,sheet=sheet,row.names = F)
    print("pathways_from_unique_genes_in_set2_df")
    #*************************************
    #adding pathways data to workbook
    #*************************************
    sheet = createSheet(wb,sheetName="common_pathways")
    addDataFrame(x=common_pathways_df,sheet=sheet,row.names = F)
    print("common_pathways_df")
    sheet = createSheet(wb,sheetName=paste("only",set1_name,"keggPthwy",sep = "_"))
    addDataFrame(x=set1_only_pathways_df,sheet=sheet,row.names = F)
    print("set1_only_pathways_df")
    sheet = createSheet(wb,sheetName=paste("only",set2_name,"keggPthwy",sep = "_"))
    addDataFrame(x=set2_only_pathways_df,sheet=sheet,row.names = F)
    print("set2_only_pathways_df")
    
    #*************************************
    #adding upregulated data to workbook
    #*************************************
    
    sheet = createSheet(wb,sheetName="common_upregulated")
    addDataFrame(x=common_up_degs_df,sheet=sheet,row.names = FALSE)
    
    sheet = createSheet(wb,sheetName="pthwys_common_upregulated")
    addDataFrame(x=common_up_pathways_df,sheet=sheet,row.names = FALSE)
    print("common_up_pathways_df")
    sheet = createSheet(wb,sheetName=paste("only",set1_name,"up_degs",sep = "_"))
    addDataFrame(x = onlyIn_set1_upregulated_degs_df, sheet=sheet,row.names = FALSE)
    
    sheet = createSheet(wb,sheetName=paste("pthwysOnly",set1_name,"up",sep = "_"))
    addDataFrame(x = set1_only_up_pathways_df, sheet=sheet,row.names = FALSE)
    print("set1_only_up_pathways_df")
    sheet = createSheet(wb,sheetName=paste("onlyIn",set2_name,"up_degs",sep = "_"))
    addDataFrame(x = onlyIn_set2_upregulated_degs_df, sheet=sheet,row.names = FALSE)
    
    sheet = createSheet(wb,sheetName=paste("pthwysOnly",set2_name,"up",sep = "_"))
    addDataFrame(x = set2_only_up_pathways_df, sheet=sheet,row.names = FALSE)
    print("set2_only_up_pathways_df")
    #*************************************
    #adding downregulated data to workbook
    #*************************************
    sheet = createSheet(wb,sheetName="common_downregulated")
    addDataFrame(x=common_down_degs_df,sheet=sheet,row.names = F)
    
    sheet = createSheet(wb,sheetName="pthwys_common_downregulated")
    addDataFrame(x=common_down_pathways_df,sheet=sheet,row.names = FALSE)
    print("common_down_pathways_df")
    
    sheet = createSheet(wb,sheetName=paste("only",set1_name,"down_degs",sep = "_"))
    addDataFrame(x = onlyIn_set1_downregulated_degs_df, sheet=sheet,row.names = F)
    
    sheet = createSheet(wb,sheetName=paste("pthwysOnly",set1_name,"down",sep = "_"))
    addDataFrame(x = set1_only_down_pathways_df, sheet=sheet,row.names = F)
    print("set1_only_down_pathways_df")
    
    sheet = createSheet(wb,sheetName=paste("only",set2_name,"down_degs",sep = "_"))
    addDataFrame(x = onlyIn_set2_downregulated_degs_df, sheet=sheet,row.names = F)
    
    
    sheet = createSheet(wb,sheetName=paste("pthwysOnly",set2_name,"down",sep = "_"))
    addDataFrame(x = set2_only_down_pathways_df, sheet=sheet,row.names = FALSE)
    print("set2_only_down_pathways_df")
    #*************************************
    #setting the excel workbook name
    #*************************************
    excel_file_name=paste("mega_excelReport_venn_",set1_name,"_",set2_name,".xlsx",sep = "")
    excel_file_path=getPath(filename = excel_file_name,dir=excel_report_output_dir)
    
    saveWorkbook(wb,file=excel_file_path)
  }
  #######################################
  # get_degs_with_pathways
  # takes in a degs csv and generates an excel report with up regulated and downregulated results sparately along with pathways
  #######################################
  get_degs_with_pathways=function(degs_csv_path,excel_report_output_dir,excel_report_filename)
  {
    library(dplyr)
    library(xlsx)
    my_degs_df=read.csv(file=degs_csv_path,stringsAsFactors = F,header = T)
    my_degs_df=my_degs_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    my_degs_df = my_degs_df %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::filter(padj < 0.05)
    my_degs_df = my_degs_df %>% dplyr::filter(!is.na(padj), gene_symbol!="")
    my_degs_df = my_degs_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
    
    #print(head(my_degs_df))
    up_degs_df=my_degs_df %>% dplyr::filter(log2FoldChange >= 1)
    down_degs_df=my_degs_df %>% dplyr::filter(log2FoldChange <= -1)
    #getting the details dataframe
    up_down_details_df=get_up_down_details_df(my_degs_df)
    
    #preparing the gene list as required for the pathway enrichment analysis
    for_gsea_df=my_degs_df
    for_gsea_df=for_gsea_df %>% dplyr::arrange(log2FoldChange)
    #for_gsea_annotation_df=getAnnotationFromBiomart(ids=for_gsea_df$gene_symbol,filter = "hgnc_symbol",columns = c("hgnc_symbol","entrezgene_id"))
    #for_gsea_df=for_gsea_df %>% inner_join(for_gsea_annotation_df,c("gene_symbol"="hgnc_symbol"))
    for_gsea_gene_list=for_gsea_df$log2FoldChange
    names(for_gsea_gene_list)=for_gsea_df$entrezgene_id
    up_gene_list=for_gsea_gene_list[for_gsea_gene_list>0]
    down_gene_list=for_gsea_gene_list[for_gsea_gene_list<0]
    
    #adding results to workbook
    wb=createWorkbook()
    
    sheet=createSheet(wb,sheetName="details")
    addDataFrame(x=as.data.frame(up_down_details_df),sheet=sheet,row.names = F)
    
    sheet=createSheet(wb,sheetName = "all_degs")
    addDataFrame(x=as.data.frame(my_degs_df),sheet=sheet,row.names = F)
    try(doGSEA(geneList = for_gsea_gene_list,wb=wb,sheetNamePrefix = "pthwys_all_genes_"),silent=T)
    
    sheet=createSheet(wb,sheetName = "up_degs")
    addDataFrame(x=as.data.frame(up_degs_df),sheet=sheet,row.names = F)
    try(doGSEA(geneList = up_gene_list,wb=wb,sheetNamePrefix = "pthwys_up_gene_"),silent=T)
    
    sheet=createSheet(wb,sheetName = "down_degs")
    addDataFrame(x=as.data.frame(down_degs_df),sheet=sheet,row.names = F)
    try(doGSEA(geneList = down_gene_list,wb=wb,sheetNamePrefix = "pthwys_down_gene_"),silent=T)
    
    excel_report_filepath=getPath(filename = excel_report_filename,directory = excel_report_output_dir)
    saveWorkbook(wb,file=excel_report_filepath)
    gc()
  }
  
  #############################################
  # get_up_down_details_df
  # this function return the details of upregulated and downregulated genes as a dataframe
  # takes the actual degs scv as input
  #############################################
  get_up_down_details_df=function(degs_df)
  {
    all=nrow(degs_df)
    upregulated=nrow(degs_df[degs_df$log2FoldChange >= 1,])
    downregulated=nrow(degs_df[degs_df$log2FoldChange <= -1,])
    up_down_details_df=data.frame()
    up_down_details_df=addRowToDetailsDf(df = up_down_details_df,detail = "total differentially expressed genes",value =all )
    up_down_details_df=addRowToDetailsDf(df = up_down_details_df,detail = "total upregulated genes",value =upregulated )
    up_down_details_df=addRowToDetailsDf(df = up_down_details_df,detail = "total downregulated genes",value = downregulated )
    
    return(up_down_details_df)
  }
  
  ################################################################################
  # get_annotated_deg_df 
  # returns an annotated differential expression results dataframe
  ################################################################################
  get_annotated_deg_df=function(deg_file_path,degs_df=NA)
  {
    library(dplyr)
    library(stringr)
    if(nrow(degs_df)==0)
    {
      return(degs_df)
    }
    if(is.na(degs_df))
    {
      if(str_ends(deg_file_path,"csv"))
      {
        degs_df=read.csv(file = deg_file_path,header = T,stringsAsFactors = F)
      }
      if(str_ends(deg_file_path,"xlsx") | str_ends(deg_file_path,"xls"))
      {
        degs_df=read_excel(path =deg_file_path,sheet = 1)
      }
    }
    degs_df=filterByPadjValue(degs_df,cutoff = 0.05)
    print(head(degs_df))
    # columns=colnames(degs_df)
    # columns[1]="ensembl_gene_id"
    # colnames(degs_df)=columns
    # print(head(degs_df))
    #names(degs_df)
    #rownames(degs_df)=degs_df$ensembl_gene_id
    my_ensemble_ids=rownames(degs_df)
    clean_ensemble_ids=getCleanEnsembleIds(my_ensemble_ids)
    degs_df$clean_ensemble_id=clean_ensemble_ids
    degs_df$ensemble_id=my_ensemble_ids
    #length(clean_ensemble_ids)
    #length(unique(clean_ensemble_ids))
    
    #**************************************************
    # fetching annotation from biomart
    #**************************************************
    columns=c("ensembl_gene_id","entrezgene_id","hgnc_symbol","description","gene_biotype")
    annotation_df=getAnnotationFromBiomart(ids = clean_ensemble_ids,filter = "ensembl_gene_id",columns = columns)
    nrow(annotation_df)
    
    
    allowed_biotypes=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
    
    annotation_df=annotation_df %>% dplyr::filter(gene_biotype %in% allowed_biotypes)
    #sum(duplicated(annotation_df$ensembl_gene_id))
    
    
    annotation_summary=as.data.frame(table(annotation_df$ensembl_gene_id))
    #names(annotation_summary)
    #annotation_summary %>% dplyr::filter(Freq>1)
    
    # duplicated_genes=annotation_df[duplicated(annotation_df$ensembl_gene_id),1]
    # annotation_df[annotation_df$ensembl_gene_id %in% duplicated_genes,]
    # degs_df[degs_df$clean_ensemble_id %in% duplicated_genes,]
    
    
    #little experimet
    #degs_df[degs_df$clean_ensemble_id=="ENSG00000143226",] %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id"))
    
    #nrow(degs_df %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id")))
    print(head(degs_df))
    print(head(annotation_df))
    print(class(degs_df))
    print(class(annotation_df))
    degs_df=degs_df %>% dplyr::inner_join(annotation_df,by=c("clean_ensemble_id"="ensembl_gene_id"))
    names(degs_df)
    print(head(degs_df))
    degs_df=degs_df %>% dplyr::select(hgnc_symbol,description,log2FoldChange,padj,pvalue,ensemble_id,clean_ensemble_id,entrezgene_id) %>% dplyr::rename(gene_symbol=hgnc_symbol) %>% arrange(padj)
    
    return(degs_df)
  }
  
  ################################################################################
  # get_modified_counts_df 
  # modified the counts matrix by replacing the outlier with the mean of the other two
  # to_modify is the hour which has to be modified. default is set to 0 hr
  ################################################################################
  get_modified_counts_df=function(counts_df,to_modify=0)
  {
    print("**************************modifying the counts data**************************")
    #*********************************************
    #finding out all the outlier gene locations at 0 hours i.e the 1st half
    #*********************************************
    print("original row order")
    print(rownames(counts_df))
    zero_hour_outlier_genes=apply(counts_df,1,myZeroHourOutlierDetector)
    #print(zero_hour_outlier_genes)
    total_zero_hour_outliers=sum(zero_hour_outlier_genes,na.rm = T)
    print(paste("total_zero_hour_outliers",total_zero_hour_outliers))
    
    #*********************************************
    #finding out all the outlier gene locations at 2 hours i.e the 2nd half
    #*********************************************
    two_hour_outlier_genes=apply(counts_df,1,myTwoHourOutlierDetector)
    total_two_hour_outliers=sum(two_hour_outlier_genes,na.rm = T)
    print(paste("total_two_hour_outliers",total_two_hour_outliers))
    
    #*********************************************
    #constructing the dataframe of outlier genes at 0 hrs only
    #*********************************************
    at_zero_hours_only_df=counts_df[zero_hour_outlier_genes & !two_hour_outlier_genes,]
    cat("rownames(at_zero_hours_only_df)\n")
    print(rownames(at_zero_hours_only_df))
    #cat("at_zero_hours_only_df")
    #print(at_zero_hours_only_df)
    at_two_hours_only_df=counts_df[!zero_hour_outlier_genes & two_hour_outlier_genes,]
    cat("rownames(at_two_hours_only_df)\n")
    print(rownames(at_two_hours_only_df))
    #cat("at_two_hours_only_df")
    #print(at_two_hours_only_df)
    #*********************************************
    #constructing the dataframe of outlier genes at 0 hrs only
    #*********************************************
    not_at_zero_hours_only_df=counts_df[!(zero_hour_outlier_genes & !two_hour_outlier_genes),]
    not_at_two_hours_only_df=counts_df[!(!zero_hour_outlier_genes & two_hour_outlier_genes),]
    #*********************************************
    #constructing modified count matrix data
    #*********************************************
    modified_counts_df=counts_df
    if(to_modify==0 & nrow(at_zero_hours_only_df)>0)
    {
      #print(paste("to_modify",to_modify))
      modified_at_zero_hours_only_df=at_zero_hours_only_df
      #cat("modified_at_zero_hours_only_df\n")
      #print(modified_at_zero_hours_only_df)
      #print(nrow(modified_at_zero_hours_only_df))
      for(i in 1:nrow(at_zero_hours_only_df))
      {
        modified_at_zero_hours_only_df[i,apply(modified_at_zero_hours_only_df[i,],1,getOutlierIndex,to_modify)]=apply(modified_at_zero_hours_only_df[i,],1,myRowMean)
      }
      modified_counts_df=as.data.frame(rbind(modified_at_zero_hours_only_df,not_at_zero_hours_only_df))
      
    }
    
    if(to_modify==2 & nrow(at_two_hours_only_df)>0)
    {
      #print(paste("to_modify",to_modify))
      #cat("from get_modified_counts_df \n")
      modified_at_two_hours_only_df=at_two_hours_only_df
      #print(modified_at_two_hours_only_df)
      for(i in 1:nrow(at_two_hours_only_df))
      {
        modified_at_two_hours_only_df[i,apply(modified_at_two_hours_only_df[i,],1,getOutlierIndex,to_modify=to_modify)]=apply(modified_at_two_hours_only_df[i,],1,myRowMean)
      }
      
      modified_counts_df=as.data.frame(rbind(modified_at_two_hours_only_df,not_at_two_hours_only_df))
    }
    
    print(paste("nrow(counts_df)",nrow(counts_df)))
    print(paste("nrow(modified_counts_df)",nrow(modified_counts_df)))
    modified_counts_df=modified_counts_df[order(as.numeric(rownames(modified_counts_df))),]
    return(modified_counts_df)
  }
  
  #####################################################
  # get_differential_expression_results_df
  # takes in a counts dataframe and returns the differential expression results as a dataframe
  #####################################################
  get_differential_expression_results_df=function(counts_df, filter_by_padj = TRUE)
  {
    library(stringr)
    library(DESeq2)
    
    
    
    print(paste("nrow(counts_df[rowSums(counts_df)<10,])",nrow(counts_df[rowSums(counts_df)<10,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<6,])",nrow(counts_df[rowSums(counts_df)<6,])))
    print(paste("nrow(counts_df[rowSums(counts_df)<3,])",nrow(counts_df[rowSums(counts_df)<2,])))
    
    
    #counts_df=counts_df[rowSums(counts_df)>rowsum_cutoff,]
    
    columns=colnames(counts_df)
    
    #preparing the details dataframe
    sam_file_details=as.data.frame(cbind(columns[1:(ncol(counts_df)/2)],columns[((ncol(counts_df)/2)+1): ncol(counts_df)]))
    colnames(sam_file_details)=c("control","treatment")
    
    
    #***********************************************************
    #conducting differential expression analysis
    #***********************************************************
    colData = DataFrame(condition=factor(c(rep("control", ncol(counts_df)/2),rep("treatment", ncol(counts_df)/2))))
    dds = DESeqDataSetFromMatrix(counts_df, colData, formula(~ condition))
    
    # run DEseq
    res= DESeq(dds, minReplicatesForReplace = 3)
    
    #plotMA(res)
    
    #diffrentiallyExpressed 
    deGenes = results(res,alpha=0.05)
    
    # order by BH adjusted p-value
    deGenesOrdered = deGenes[order(deGenes$padj),]
    if(filter_by_padj){
      deGenesOrdered=filterByPadjValue(df=deGenesOrdered,cutoff = 0.05)
    }
    return(as.data.frame(deGenesOrdered))
  }
  
  #####################################################
  # get_standardised_differential_expression_df
  # takes in a raw differential expression output result and returns standardised differential expression dataframe from csv file
  #####################################################
  get_standardised_differential_expression_df=function(degs_csv_filepath)
  {
    library(dplyr)
    #reading the degs file
    degs_df=read.csv(file = degs_csv_filepath,header = T,stringsAsFactors = F)
    degs_df = degs_df[order(degs_df$padj),]
    
    print(head(degs_df))
    degs_df=filterByPadjValue(df=degs_df,cutoff = 0.05)
    if(nrow(degs_df)>0)
    {
      rownames(degs_df)=degs_df$X
      degs_df=degs_df%>% dplyr::select(-1)
      print(degs_df)
      degs_df=get_annotated_deg_df(degs_df = degs_df)
      print(head(degs_df))
      
      #filtering the data
      degs_df=degs_df %>% dplyr::select(gene_symbol,description,log2FoldChange,padj,pvalue,entrezgene_id,ensemble_id,clean_ensemble_id)
    }
    return(as.data.frame(degs_df))
  }
  
  ########################################################
  # modifyAndGetDifferentialExpression
  # takes in the original csv file and returns the differential expression results of a modified counts data . you can specify which sample set to modify. 0 for 1st sample set and 2 for the other sample set
  ########################################################
  modifyAndGetDifferentialExpression=function(counts_dir,original_counts_csv_filename,differential_expression_results_dir,differential_expression_filename,to_modify=0)
  {
    #loading original data
    original_counts_csv_file_path=getPath(filename = original_counts_csv_filename,directory = counts_dir)
    original_counts_df=read.csv(file=original_counts_csv_file_path,header=T,stringsAsFactors = F)
    
    #setting rownames. required for the deseq analysis
    rownames(original_counts_df)=original_counts_df[,1]
    
    #removing the rowname which was part of the dataframe
    original_counts_df=original_counts_df %>% dplyr::select(-1)
    
    #modifying the counts data
    modified_counts_df=get_modified_counts_df(original_counts_df,to_modify = to_modify)
    
    #getting the differential expression
    modified_diff_expr_result_df=get_differential_expression_results_df(counts_df = modified_counts_df)
    
    #annotating the counts data
    print(head(modified_diff_expr_result_df))
    
    if(nrow(modified_diff_expr_result_df)>0)
    {
      modified_diff_expr_result_df=get_annotated_deg_df(degs_df = modified_diff_expr_result_df)
      head(modified_diff_expr_result_df)
      
      #filtering the data
      modified_diff_expr_result_df=modified_diff_expr_result_df %>% dplyr::select(gene_symbol,description,log2FoldChange,padj,pvalue,entrezgene_id,ensemble_id,clean_ensemble_id)
      
      #saving the resuts
      differential_expression_results_file_path=getPath(filename = differential_expression_filename,directory = differential_expression_results_dir)
    }
    write.csv(modified_diff_expr_result_df,file = differential_expression_results_file_path,row.names = F)
    return(modified_diff_expr_result_df)
  }
  
  ###################################################################################################################
  #write_gene_behavior_report helps in building up a report with counts and timewise behavior  of a gene
  # adds the counts df and the behavior plot to a workbook object
  ###################################################################################################################
  write_gene_behavior_report=function(gene_name, output_directory, wb, behaviorIn = c("both", "igra", "tb"),data = c("modified","original"))
  {
    print("entered write_gene_behavior_report")
    #**********************************************************************
    #loading required packages 
    #**********************************************************************
    print("loading required packages ")
    library(ggplot2)
    library(openxlsx)
    library(dplyr)
    #**********************************************************************
    #matching argument 
    #**********************************************************************
    print("matching argument ")
    behaviorIn = match.arg(behaviorIn)
    data = match.arg(data)
    #**********************************************************************
    #setting a few decision making variables
    #**********************************************************************
    load_igra_counts = FALSE
    load_tb_counts = FALSE
    add_igra_sheet = FALSE
    add_tb_sheet = FALSE
    #*******************************************************************************************
    #initialising variables
    #*******************************************************************************************
    modified_igra_negative_counts_df = NULL
    modified_tb_counts_df = NULL
    igra_negative_counts_df = NULL
    tb_counts_df = NULL
    igra_plot = NULL
    tb_plot = NULL
    both_plot = NULL
    normal_sample_wise_counts_plot_path = NULL
    
    #*******************************************************************************************
    # making decisions on which files to load
    #*******************************************************************************************
    if(behaviorIn == "both")
    {
      load_igra_counts = TRUE
      load_tb_counts = TRUE
      add_igra_sheet = TRUE
      add_tb_sheet = TRUE
    }else if(behaviorIn == "igra")
    {
      load_igra_counts = TRUE
      add_igra_sheet = TRUE
    }else if(behaviorIn == "tb")
    {
      load_tb_counts = TRUE
      add_tb_sheet = TRUE
    }else{
      stop("please provide a correct value for behaviorIn")
      return()
    }
    #**********************************************************************
    #creating log folders to save the folders
    #**********************************************************************
    logs_path = getPath(filename = "logs",directory = output_directory)
    dir.create(logs_path)
    
    #*******************************************************************************************
    #loading necessary data
    #*******************************************************************************************
    print("loading necessary data")
    if(load_igra_counts)
    {
      if(data == "original")
      {
        modified_igra_negative_counts_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/counts_data/original_unmodified_gene_counts.csv")
      }else if(data == "modified"){
        modified_igra_negative_counts_df= read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/counts_data/sco71_replaced_with_mean_gene_counts.csv")
      }
      
      igra_negative_counts_df = modified_igra_negative_counts_df %>% dplyr::filter(gene_symbol == gene_name) 
      print("igra_negative_counts")
      print(head(igra_negative_counts_df))
      # modified_igra_negative_counts_df= read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/counts_data/sco71_replaced_with_mean_gene_counts.csv")
      # igra_negative_counts_df = modified_igra_negative_counts_df %>% dplyr::filter(gene_symbol == gene_name) 
      # print(head(igra_negative_counts_df))
    }
    if(load_tb_counts)
    {
      modified_tb_counts_df= read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/counts_data/tb_counts_data.csv")
      tb_counts_df = modified_tb_counts_df %>% dplyr::filter(gene_symbol == gene_name)
      print(head(tb_counts_df))
    }
    
    #**********************************************************************
    #loading the combined counts data required for plotting the behavior of the genes
    #**********************************************************************
    print("loading the combined counts data required for plotting the behavior of the genes")
    combined_df = ""
    if(data == "original")
    {
      combined_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/original_unmodified_combined_counts_for_plotting.csv", stringsAsFactors = FALSE)
    }else if(data == "modified"){
      combined_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/combined_counts_for_plotting.csv", stringsAsFactors = FALSE)
    }
    #combined_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/combined_counts_for_plotting.csv", stringsAsFactors = FALSE)
    if(load_igra_counts & load_tb_counts)
    {
      #**********************************************************************
      #generating the combined plots
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name) %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type + sample_type)+ ggtitle(gene_name)
      normal_sample_wise_counts_plot_path = getPath(filename = paste(gene_name,"normal_sample_wise_counts.jpg",sep= "_"),directory = logs_path)
      ggsave(filename = normal_sample_wise_counts_plot_path,dpi = 400) 
      
      #**********************************************************************
      #adding igra negative sheet 
      #**********************************************************************
      sheet_name=paste(gene_name,"igra_negative_counts",sep = "_")
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb=wb,sheet=sheet_name,x=igra_negative_counts_df)
      insertImage(wb = wb,sheet = sheet_name, startCol = 3, startRow = 6, file = normal_sample_wise_counts_plot_path,width=30,height = 20, units = "cm" ,dpi = 400)
      
      #**********************************************************************
      #adding tb sheet
      #**********************************************************************
      sheet_name = paste(gene_name,"tb_counts",sep = "_")
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb=wb,sheet=sheet_name,x=tb_counts_df)
      insertImage(wb,sheet = sheet_name, startCol = 3, startRow = 6 , file = normal_sample_wise_counts_plot_path,width=30,height = 20, units = "cm",dpi = 400 )
      
    }else if(load_igra_counts){
      #**********************************************************************
      #generating plots for igra negative type only
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name, sample_type == "igra_negative") %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type)+ ggtitle(gene_name)
      normal_sample_wise_counts_plot_path = getPath(filename = paste(gene_name,"normal_sample_wise_counts.jpg",sep= "_"),directory = logs_path)
      ggsave(filename = normal_sample_wise_counts_plot_path,dpi = 400) 
      
      #**********************************************************************
      #adding igra negative counts and plots only
      #**********************************************************************
      sheet_name=paste(gene_name,"igra_negative_counts",sep = "_")
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb=wb,sheet=sheet_name,x=igra_negative_counts_df)
      insertImage(wb = wb,sheet = sheet_name, startCol = 3, startRow = 6, file = normal_sample_wise_counts_plot_path,width=30,height = 20, units = "cm" ,dpi = 400)
      
    }else if(load_tb_counts){
      #**********************************************************************
      #generating plot for tb type data only
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name, sample_type == "tb") %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type) + ggtitle(gene_name)
      normal_sample_wise_counts_plot_path = getPath(filename = paste(gene_name,"normal_sample_wise_counts.jpg",sep= "_"),directory = logs_path)
      ggsave(filename = normal_sample_wise_counts_plot_path,dpi = 400)
      
      #**********************************************************************
      #attaching the counts fata for tb type data only
      #**********************************************************************
      sheet_name = paste(gene_name,"tb_counts",sep = "_")
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb=wb,sheet=sheet_name,x=tb_counts_df)
      insertImage(wb,sheet = sheet_name, startCol = 3, startRow = 6 , file = normal_sample_wise_counts_plot_path,width=30,height = 20, units = "cm",dpi = 400 )
    }
    print("exiting write_gene_behavior_report")
    detach("package:openxlsx", unload = TRUE) 
    #gc()
  }
  
  ###################################################################################################################
  #show_gene_behavior helps in building up a report with counts and timewise behavior  of a gene
  # adds the counts df and the behavior plot to a workbook object
  ###################################################################################################################
  show_gene_behavior=function(gene_name, behaviorIn = c("both", "igra", "tb"), data = c("modified","original"))
  {
    #**********************************************************************
    #loading required packages 
    #**********************************************************************
    library(ggplot2)
    library(openxlsx)
    library(dplyr)
    #**********************************************************************
    #matching argument 
    #**********************************************************************
    behaviorIn = match.arg(behaviorIn)
    data=match.arg(data)
    #**********************************************************************
    #setting a few decision making variables
    #**********************************************************************
    load_igra_counts = FALSE
    load_tb_counts = FALSE
    add_igra_sheet = FALSE
    add_tb_sheet = FALSE
    #*******************************************************************************************
    #initialising variables
    #*******************************************************************************************
    modified_igra_negative_counts_df = NULL
    modified_tb_counts_df = NULL
    igra_negative_counts_df = NULL
    tb_counts_df = NULL
    igra_plot = NULL
    tb_plot = NULL
    both_plot = NULL
    normal_sample_wise_counts_plot_path = NULL
    
    #*******************************************************************************************
    # making decisions on which files to load
    #*******************************************************************************************
    if(behaviorIn == "both")
    {
      load_igra_counts = TRUE
      load_tb_counts = TRUE
      add_igra_sheet = TRUE
      add_tb_sheet = TRUE
    }else if(behaviorIn == "igra")
    {
      load_igra_counts = TRUE
      add_igra_sheet = TRUE
    }else if(behaviorIn == "tb")
    {
      load_tb_counts = TRUE
      add_tb_sheet = TRUE
    }else{
      stop("please provide a correct value for behaviorIn")
      return()
    }
  
    #*******************************************************************************************
    #loading necessary data
    #*******************************************************************************************
    if(load_igra_counts)
    {
      modified_igra_negative_counts_df=""
      if(data == "original")
      {
        modified_igra_negative_counts_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/counts_data/original_unmodified_gene_counts.csv")
      }else if(data == "modified"){
        modified_igra_negative_counts_df= read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/counts_data/sco71_replaced_with_mean_gene_counts.csv")
      }
      
      igra_negative_counts_df = modified_igra_negative_counts_df %>% dplyr::filter(gene_symbol == gene_name) 
      print("igra_negative_counts")
      print(head(igra_negative_counts_df))
    }
    if(load_tb_counts)
    {
      modified_tb_counts_df= read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/counts_data/tb_counts_data.csv")
      tb_counts_df = modified_tb_counts_df %>% dplyr::filter(gene_symbol == gene_name) 
      print("tb counts")
      print(head(tb_counts_df))
    }
    
    #**********************************************************************
    #loading the combined counts data required for plotting the behavior of the genes
    #**********************************************************************
    combined_df = ""
    if(data == "original")
    {
      combined_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/original_unmodified_combined_counts_for_plotting.csv", stringsAsFactors = FALSE)
    }else if(data == "modified"){
      combined_df = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/combined_counts_for_plotting.csv", stringsAsFactors = FALSE)
    }
    if(load_igra_counts & load_tb_counts)
    {
      #**********************************************************************
      #generating the combined plots
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name) %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type + sample_type)+ ggtitle(gene_name)
      print(normal_sample_wise_counts_plot)
      
    }else if(load_igra_counts){
      #**********************************************************************
      #generating plots for igra negative type only
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name, sample_type == "igra_negative") %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type)+ ggtitle(gene_name)
      print(normal_sample_wise_counts_plot)
      
    }else if(load_tb_counts){
      #**********************************************************************
      #generating plot for tb type data only
      #**********************************************************************
      normal_sample_wise_counts_plot = combined_df %>% dplyr::filter(gene_symbol == gene_name, sample_type == "tb") %>% ggplot(aes(x = factor(time), y = counts, color = donor_id)) + geom_point(alpha = 0.4) +geom_line(aes(group = sample_id)) + facet_wrap(~treatment_type) + ggtitle(gene_name)
      print(normal_sample_wise_counts_plot)
    }
    
    detach("package:openxlsx", unload = TRUE) 
    gc()
  }
  
  #**************************************************************
  #fix_annotated_counts_df_column_names
  #fixes the columns of annotated counts matrix 
  #**************************************************************
  fix_igra_negative_annotated_counts_df_column_names = function(counts_df,hour)
  {
    rownames(counts_df) = counts_df[,"hgnc_symbol"]
    counts_df = counts_df %>% dplyr::select(4:9)
    columns=colnames(counts_df)
    columns[str_detect(columns,pattern = "SCO149")]="SCO149"
    columns[str_detect(columns,pattern = "SCO68")]="SCO68"
    columns[str_detect(columns,pattern = "SCO71")]="SCO71"
    columns[1:3]=paste("scr_",hour,"hr_",columns[1:3],sep="")
    columns[4:6]=paste("b3_",hour,"hr_",columns[4:6],sep="")
    colnames(counts_df)=columns
    print(head(counts_df))
    counts_df = as.matrix(counts_df)
    return(counts_df)
  }

  fix_igra_negative_unannotated_counts_df_column_names = function(counts_df,hour)
  {
    rownames(counts_df) = counts_df[,1]
    counts_df = counts_df[,-1]
    columns=colnames(counts_df)
    columns[str_detect(columns,pattern = "SCO149")]="SCO149"
    columns[str_detect(columns,pattern = "SCO68")]="SCO68"
    columns[str_detect(columns,pattern = "SCO71")]="SCO71"
    columns[1:3]=paste("scr_",hour,"hr_",columns[1:3],sep="")
    columns[4:6]=paste("b3_",hour,"hr_",columns[4:6],sep="")
    colnames(counts_df)=columns
    print(head(counts_df))
    counts_df = as.matrix(counts_df)
    return(counts_df)
  }
  
  median_based_count_matrix_filtering = function(counts_matrix, expected_median = 10)
  {
    library(dplyr)
    library(ggplot2)
    
    control_indices = 1:(ncol(counts_matrix)/2)
    treatment_indices = ((ncol(counts_matrix)/2)+1):ncol(counts_matrix)
    median_control = apply(counts_matrix, 1, function(x){return(median(x[control_indices]))})
    median_treatment = apply(counts_matrix, 1, function(x){return(median(x[treatment_indices]))})
    temp_counts_df = cbind.data.frame(counts_matrix, median_control, median_treatment)
    print(temp_counts_df %>% ggplot(aes(x = log10(median_control+1), y = log10(median_treatment+1))) + geom_point(alpha = 0.3)+geom_vline(xintercept = 1, color = "red")+geom_hline(yintercept = 1, color = "blue") + ggtitle("Before removing low counts"))
    
    selected_rows = apply(counts_matrix,1,rowwise_median_filter,expected_median)
    print(nrow(counts_matrix))
    not_selected_count_matrix = counts_matrix[!selected_rows,]
    counts_matrix = counts_matrix[selected_rows,]
    print("not_selected_rows below")
    print(nrow(counts_matrix))
    print(head(not_selected_count_matrix))
    print(ncol(counts_matrix))
    median_control = apply(counts_matrix, 1, function(x){return(median(x[control_indices]))})
    median_treatment = apply(counts_matrix, 1, function(x){return(median(x[treatment_indices]))})
    temp_counts_df = cbind.data.frame(counts_matrix, median_control, median_treatment)
    print(temp_counts_df %>% ggplot(aes(x = log10(median_control+1), y = log10(median_treatment+1))) + geom_point(alpha = 0.3)+geom_vline(xintercept = 1, color = "red")+geom_hline(yintercept = 1, color = "blue") + ggtitle("After removing low counts"))
    return(counts_matrix)
  }
  
rowwise_median_filter = function(counts, expected_median = 10)
{
  control_indices = 1:(length(counts)/2)
  treatment_indices = ((length(counts)/2)+1):length(counts)
  return(median(counts[control_indices])>=expected_median | median(counts[treatment_indices])>=expected_median)
}
