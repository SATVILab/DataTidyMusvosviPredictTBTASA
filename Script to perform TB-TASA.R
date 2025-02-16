rm(list=ls())
graphics.off()
library(openCyto)
library(flowCore)
library(ggcyto)
library(magick)
library(CytoML)
library(rio)
library(PeacoQC)

fcs_location = "" #Absolute path of folder containing .fcs files and a .wsp file with the comp matrix
results_location = "" #Absolute path of folder to save the data.


  openCyto::registerPlugins(fun = cytoUtils:::.gate_tail, methodName = "gate_tail")

  fcs_files = list.files(fcs_location,full.names = T,pattern = ".fcs")
  fcs_files = fcs_files[!grepl("Compensation Controls",fcs_files,fixed = T)]
  fcs_files = fcs_files[grepl(".fcs",fcs_files,fixed = T)]
  wsp_location = list.files(fcs_location,full.names = T,pattern = ".wsp")

  if(length(wsp_location)>1){
    wsp_location = wsp_location[order(file.info(wsp_location)$mtime,decreasing = T)][1]#Selects the most up to date .wsp

  }

  fcs_temp = load_cytoset_from_fcs(fcs_files[1])
  markers_in_fcs_files = as.data.frame(markernames(fcs_temp))
  colnames(markers_in_fcs_files) = "marker"

  if(sum(grepl("UVID|Live/Dead|Dump|DUMP|Live Dead|Viability|Aqua",markers_in_fcs_files$marker),na.rm = T)>=1){
    temp_viability = markers_in_fcs_files$marker[grepl("UVID|Live/Dead|Live Dead|Viability|Aqua",markers_in_fcs_files$marker)]

    base_lineage_gating =data.frame(
      alias = c("boundary", "singlets", "SSCneg", "lymph","viable", "cd3.gate"),
      pop = c("+", "+", "-", "+","-", "+"),
      parent = c("root", "boundary", "singlets", "SSCneg", "lymph","viable"),
      dims = c("FSC-A,FSC-H", "FSC-A,FSC-H", "SSC-A", "FSC-A,SSC-A",temp_viability, "CD3"),
      gating_method = c("boundary", "singletGate", "gate_mindensity2", "gate_flowclust_2d", "gate_mindensity2","gate_mindensity2"),
      gating_args = c("max=c(256000,256000),min=c(20000,20000)", "", "peaks=2,min=30000,max=230000,gate_range = c(40000,200000)", "K=1,quantile=.9","min=1.2,max=4,gate_range=c(2,3.5)", "min=0.5,max=4,gate_range=c(1,4)"),
      collapseDataForGating = rep(TRUE, 6),
      groupBy = rep("pid", 6),
      preprocessing_method = c("ppmyGate", "ppmyGate", "ppmyGate", NA, "ppmyGate","ppmyGate"),
      preprocessing_args = rep(NA, 6),
      stringsAsFactors = FALSE)

    hladr = data.frame(
      alias = c("hladr"),
      pop = c("+"),
      parent = c("viable"),
      dims = c("HLA-DR"),
      gating_method = c("gate_mindensity2"),
      gating_args = c("peaks= 2, min=1, max=4,gate_range=c(1,3.5)"),
      collapseDataForGating = c(TRUE),
      groupBy = c("pid"),
      preprocessing_method = c("ppmyGate"),
      preprocessing_args = c(NA),
      stringsAsFactors = FALSE)

    activated_IFNGandTNF_cd3 = data.frame(
      alias = c("activated_IFNGandTNF_cd3"),
      pop = c("+"),
      parent = c("cd3.gate"),
      dims = c("HLA-DR,IFNG"),
      gating_method = c("boolGate"),
      gating_args = c("cd3.gate/cd3_IFNGandTNF&viable/hladr"),
      collapseDataForGating = c(TRUE),
      groupBy = c("pid"),
      preprocessing_method = c("ppmyGate"),
      preprocessing_args = c(NA),
      stringsAsFactors = FALSE)

  }else{
    base_lineage_gating =data.frame(
      alias = c("boundary", "singlets", "SSCneg", "lymph", "cd3.gate"),
      pop = c("+", "+", "-", "+", "+"),
      parent = c("root", "boundary", "singlets", "SSCneg", "lymph"),
      dims = c("FSC-A,FSC-H", "FSC-A,FSC-H", "SSC-A,", "FSC-A,SSC-A", "CD3"),
      gating_method = c("boundary", "singletGate", "gate_mindensity2", "gate_flowclust_2d", "gate_mindensity2"),
      gating_args = c("max=c(256000,256000),min=c(20000,20000)", "", "peaks=2,min=30000,max=230000,gate_range = c(40000,200000)", "K=1,quantile=.9", "min=0.5,max=4,gate_range=c(1,4)"),
      collapseDataForGating = rep(TRUE, 5),
      groupBy = rep("pid", 5),
      preprocessing_method = c("ppmyGate", "ppmyGate", "ppmyGate", NA, "ppmyGate"),
      preprocessing_args = rep(NA, 5),
      stringsAsFactors = FALSE)


    hladr = data.frame(
      alias = c("hladr"),
      pop = c("+"),
      parent = c("lymph"),
      dims = c("HLA-DR"),
      gating_method = c("gate_mindensity2"),
      gating_args = c("peaks= 2, min=1, max=4,gate_range=c(1,3.5)"),
      collapseDataForGating = c(TRUE),
      groupBy = c("pid"),
      preprocessing_method = c("ppmyGate"),
      preprocessing_args = c(NA),
      stringsAsFactors = FALSE)

    activated_IFNGandTNF_cd3 = data.frame(
      alias = c("activated_IFNGandTNF_cd3"),
      pop = c("+"),
      parent = c("cd3.gate"),
      dims = c("HLA-DR,IFNG"),
      gating_method = c("boolGate"),
      gating_args = c("cd3.gate/cd3_IFNGandTNF&lymph/hladr"),
      collapseDataForGating = c(TRUE),
      groupBy = c("pid"),
      preprocessing_method = c("ppmyGate"),
      preprocessing_args = c(NA),
      stringsAsFactors = FALSE)
  }

  cytokine_IFNG = data.frame(
    alias = c("cd3IFNG"),
    pop = c("+"),
    parent = c("cd3.gate"),
    dims = c("IFNG"),
    gating_method = c("gate_tail"),
    gating_args = c("auto_tol = TRUE, bias = 0.3"),
    collapseDataForGating = c(TRUE),
    groupBy = c("pid"),
    preprocessing_method = c("ppmyGate"),
    preprocessing_args = c(NA),
    stringsAsFactors = FALSE)

  cytokine_TNF = data.frame(
    alias = c("cd3TNF"),
    pop = c("+"),
    parent = c("cd3.gate"),
    dims = c("TNF"),
    gating_method = c("gate_tail"),
    gating_args = c("auto_tol = TRUE, bias = 0.3"),
    collapseDataForGating = c(TRUE),
    groupBy = c("pid"),
    preprocessing_method = c("ppmyGate"),
    preprocessing_args = c(NA),
    stringsAsFactors = FALSE)

  cd3_IFNGandTNF = data.frame(
    alias = c("cd3_IFNGandTNF"),
    pop = c("++"),
    parent = c("cd3.gate"),
    dims = c("IFNG,TNF"),
    gating_method = c("refGate"),
    gating_args = c("cd3.gate/cd3IFNG:cd3.gate/cd3TNF"),
    collapseDataForGating = c(TRUE),
    groupBy = c("pid"),
    preprocessing_method = c("ppmyGate"),
    preprocessing_args = c(NA),
    stringsAsFactors = FALSE)





  temp_gating_strategy = rbind.data.frame(base_lineage_gating,hladr,cytokine_IFNG,cytokine_TNF,cd3_IFNGandTNF,activated_IFNGandTNF_cd3)

  temp_file <- tempfile(fileext = ".csv")
  write.csv(temp_gating_strategy, temp_file, row.names = FALSE)

  gating_strategy <- gatingTemplate(temp_file)#use
  temp_gating <- read.csv(temp_file)

  metadata <- NULL
  for(i in 1:length(fcs_files)){
    temp_meta = as.data.frame(fcs_files[i])
    colnames(temp_meta) ="file_location"
    temp_keywords =keyword(read.FCS(fcs_files[i]))
    temp_meta$batch = gsub("/","_",temp_keywords$`$SRC`,fixed = T)
    temp_meta$stimulation = gsub("/","_",temp_keywords$`TUBE NAME`,fixed = T)
    temp_meta$stimulation = gsub("\\\\","_",temp_meta$stimulation,fixed = T)
    temp_meta$file_name = temp_keywords$`$FIL`
    metadata = rbind.data.frame(metadata,temp_meta)
  }


  batches_to_analyse <- unique(metadata$batch)

  .ppmyGate <- function(fs, gs, gm, channels=NA,groupBy=NA,isCollapse=NA, ...) {
    xChannel = channels[1]
    yChannel = channels[1]
    d <- c()
    for(i in c(1:length(fs))) {
      d <- c(d,rep.int(pData(fs[i])$control,nrow(exprs(fs[[i]]))))
    }
    return(as.logical(d))
  }
  register_plugins(fun=.ppmyGate, methodName='ppmyGate', dep=NA, "preprocessing")

  .polyGate <- function(fr, pp_res, channels, filterId="polygate", ...){
    args <- list(...)
    g <- data.frame(x=args$x, y=args$y)
    colnames(g) <- channels
    flowCore::polygonGate(.gate=g, filterId=filterId)
  }
  register_plugins(fun=.polyGate, methodName='polyGate', dep=NA)

  .myGate <- function(fr, pp_res, channels=NA, filterId="ppgate", ...){
    my_gate <- tailgate(fr[pp_res,],channel=channels, filter_id=filterId, ...)
    return(my_gate)
  }
  register_plugins(fun=.myGate,methodName='myGate',dep=NA)


  for(temp_sampleBatch in batches_to_analyse){

    if(!is.null(results_location)){
      results_folder_path =  file.path(results_location, paste(temp_sampleBatch,"results",sep = "_"))
    }else{
      results_folder_path = file.path(fcs_location, paste(temp_sampleBatch,"results",sep = "_"))
    }

    dir.create(results_folder_path, showWarnings = F)
    setwd(results_folder_path)

    summary_df = NULL
    res = tryCatch({
      temp_key <- subset(metadata, metadata$batch == temp_sampleBatch)
      sample_outputPath = file.path(results_folder_path,paste0("flow_results_",temp_sampleBatch))#Where I want the
      dir.create(sample_outputPath, showWarnings = F)
      for(i in 1:dim(temp_key)[1]){
        temp_key[i,]$file_location <- fcs_files[grepl(temp_key[i,]$file_name, fcs_files,fixed = T)]
      }

      ws = open_flowjo_xml(wsp_location)
      comp_mat = flowjo_to_gatingset(ws,name="Compensation",execute = F)[1]
      comp_mat = gs_get_compensations(comp_mat[[1]])
      comp_mat = as.data.frame(comp_mat[[1]]@spillover)
      comp_mat = as.matrix(comp_mat)

      for(i in 1:dim(comp_mat)[1]){
        colnames(comp_mat)[i] = unlist(strsplit(colnames(comp_mat)[i]," :",fixed = T))[1]
        rownames(comp_mat)[i] = unlist(strsplit(rownames(comp_mat)[i]," :",fixed = T))[1]
      }

      colnames(comp_mat) = gsub("Comp-","",colnames(comp_mat),fixed = T)
      rownames(comp_mat) = gsub("Comp-","",rownames(comp_mat),fixed = T)

      fcsFiles <- temp_key$file_location
      temp_fcs  <- load_cytoset_from_fcs(fcsFiles)
      temp_key$name <- temp_key$file_name

      temp_flowset <- cytoset_to_flowSet(temp_fcs)
      pData(temp_fcs)$pid <- temp_sampleBatch

      temp_stim_order <- NULL
      for (i in pData(temp_fcs)$name) {
        temp_stim_order <- c(temp_stim_order, subset(temp_key, temp_key$name == i)$stimulation)
      }

      pData(temp_fcs)$stimulation <- temp_stim_order
      pData(temp_fcs)$study_visit <- unique(subset(temp_key, temp_key$name == i)$study_visit)

      pData(temp_fcs)$control <- ifelse(pData(temp_fcs)$stimulation%in%c("UNS","unstim","Unstim"), "TRUE", "FALSE")

      if(!is.null(keyword(temp_fcs[[1]])$SPILL)){
        if(unique(grepl("/",colnames(keyword(temp_fcs[[1]])$SPILL),fixed = T))){
          colnames(comp_mat) = gsub("_","/",colnames(comp_mat),fixed = T)
          rownames(comp_mat) = gsub("_","/",rownames(comp_mat),fixed = T)
        }
      }else{
        if(unique(grepl("/",colnames(keyword(temp_fcs[[1]])$`$SPILLOVER`),fixed = T))){
          colnames(comp_mat) = gsub("_","/",colnames(comp_mat),fixed = T)
          rownames(comp_mat) = gsub("_","/",rownames(comp_mat),fixed = T)
        }

      }




      temp_fcs <- compensate(temp_fcs, comp_mat)

      chnls <- names(temp_fcs[[1]])
      chnls <- chnls[!grepl("FSC|SSC|Time", chnls)]
      markernames(temp_fcs) = gsub("ifng|IFNg|IFN-g|IFN-y|IFNy","IFNG",markernames(temp_fcs))
      markernames(temp_fcs) = gsub("tnf|TNFa|TNF-a|TNF-aphla","TNF",markernames(temp_fcs))
      markernames(temp_fcs) = gsub("hladr|hla-dr|HLADR","HLA-DR",markernames(temp_fcs))
      markernames(temp_fcs) = gsub("cd3","CD3",markernames(temp_fcs))
      temp_ncfs <- temp_fcs
      trans <- estimateLogicle(temp_ncfs[[1]], channels = chnls)
      temp_ncfs <- transform(temp_ncfs, trans)

      temp_ncfs_qc = list()
      for(i in 1:length(temp_ncfs)){
        temp_qced = PeacoQC(temp_ncfs[[i]],channels = colnames(temp_ncfs)[colnames(temp_ncfs)!="Time"],save_fcs = F,report=T,plot=T)
        temp_ncfs_qc[[i]] <- temp_ncfs[[i]][temp_qced$GoodCells, ]
      }
      names(temp_ncfs_qc) <- rownames(pData(temp_ncfs))
      temp_ncfs <- cytoset(temp_ncfs_qc)


      temp_ncfs <- GatingSet(temp_ncfs)
      gt_gating(gating_strategy, temp_ncfs)

      for(i in 1:length(temp_ncfs)) {
        temp_df = gh_pop_compare_stats(temp_ncfs[[i]])
        temp_df$name = pData(temp_ncfs[i])$name
        temp_df$pid = pData(temp_ncfs[i])$pid
        temp_df$stimulation = pData(temp_ncfs[i])$stimulation
        summary_df = rbind.data.frame(summary_df, temp_df)
      }


      tasa_summary_df = NULL
      temp_df = subset(summary_df,summary_df$node%in%c("cd3.gate","cd3_IFNGandTNF","activated_IFNGandTNF_cd3","cd3TNF","cd3IFNG"))
      for(i_1 in unique(temp_df$stimulation)){
        temp_meta = unique.data.frame(subset(temp_df,temp_df$stimulation==i_1,select = c("pid","stimulation")))

        temp_stim = subset(temp_df,temp_df$stimulation==i_1)
        temp_unstim = subset(temp_df,temp_df$stimulation=="UNS")

        temp_fisher = rbind(c(subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.count,(subset(temp_stim,temp_stim$node=="cd3.gate")$openCyto.count-subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.count)),
                            c(subset(temp_unstim,temp_unstim$node=="cd3_IFNGandTNF")$openCyto.count,(subset(temp_unstim,temp_unstim$node=="cd3.gate")$openCyto.count-subset(temp_unstim,temp_unstim$node=="cd3_IFNGandTNF")$openCyto.count)))
        temp_fisher = fisher.test(temp_fisher)

        temp_meta$ifng_tnf_freq = subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.freq
        temp_meta$p_value = temp_fisher$p.value
        temp_meta$fold_change = subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.freq/subset(temp_unstim,temp_unstim$node=="cd3_IFNGandTNF")$openCyto.freq
        temp_meta$responder = ifelse(temp_meta$p_value<=0.01&temp_meta$fold_change>=3,"yes","no")
        temp_meta$number_of_cells = subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.count

        if(temp_meta$responder=="yes"){
          temp_meta$tasa = subset(temp_stim,temp_stim$node=="activated_IFNGandTNF_cd3")$openCyto.freq/subset(temp_stim,temp_stim$node=="cd3_IFNGandTNF")$openCyto.freq*100
        }else{
          temp_meta$tasa = NA
        }
        tasa_summary_df = rbind.data.frame(tasa_summary_df,temp_meta)
      }


      if("viable"%in%temp_gating$alias){
        gate_population_key = rbind.data.frame(
          c("boundary","FSC-A","FSC-H","cell_subsets"),
          c("singlets","FSC-A","FSC-H","cell_subsets"),
          c("SSCneg","SSC-A","FSC-A","cell_subsets"),
          c("lymph","FSC-A","SSC-A","cell_subsets"),
          c("viable","CD3",temp_viability,"cell_subsets"),
          c("cd3.gate","IFNG","CD3","cell_subsets"),
          c("cd3_IFNGandTNF","IFNG","TNF","cell_subsets"),
          c("hladr","HLA-DR","CD3","cell_subsets"))
        colnames(gate_population_key) = c("gate","xlabel","ylabel","summary_plot")
      }else{
        gate_population_key = rbind.data.frame(
          c("boundary","FSC-A","FSC-H","cell_subsets"),
          c("singlets","FSC-A","FSC-H","cell_subsets"),
          c("SSCneg","SSC-A","FSC-A","cell_subsets"),
          c("lymph","FSC-A","SSC-A","cell_subsets"),
          c("cd3.gate","IFNG","CD3","cell_subsets"),
          c("cd3_IFNGandTNF","IFNG","TNF","cell_subsets"),
          c("hladr","HLA-DR","CD3","cell_subsets"))
        colnames(gate_population_key) = c("gate","xlabel","ylabel","summary_plot")
      }

      gatesToVisualize = gate_population_key$gate

      for (i in 1:length(temp_ncfs)) {
        outputPath1 = file.path(sample_outputPath, paste0("Flowplots_", pData(temp_ncfs[[i]])$pid, "_", pData(temp_ncfs[[i]])$stimulation))
        dir.create(outputPath1, showWarnings = F)
        for (i_1 in 1:length(gatesToVisualize)) {
          temp_gate = gate_population_key[i_1,]$gate

          Temp.png.file.name = paste(pData(temp_ncfs[[i]])$pid,"_",pData(temp_ncfs[[i]])$stimulation,"_",gsub("/","_",temp_gate,fixed = T),".png",sep = "")
          Temp.png.file.name = gsub("/","_",Temp.png.file.name)
          png(filename = paste(outputPath1,"/",Temp.png.file.name,sep = ""))

          if(temp_gate%in%c("boundary","singlets","SSCneg","lymph")){
            print(autoplot(temp_ncfs[[i]],bins=200,gate=temp_gate,x=subset(gate_population_key,gate_population_key$gate==temp_gate)$xlabel,y=subset(gate_population_key,gate_population_key$gate==temp_gate)$ylabel) + geom_density2d(colour = "black")+
                    geom_stats(size = 12)+

                    theme(
                      panel.background = element_rect(fill = "white",
                                                      colour = "black",
                                                      linewidth = 0.5, linetype = "solid"),
                      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                                      colour = "white"),
                      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                                      colour = "white"),# Increase font size for axis titles
                      axis.title = element_text(size = 30),
                      # Increase font size for axis text (numbers)
                      axis.text = element_text(size = 30)
                    ))
          }else{
            print(autoplot(temp_ncfs[[i]],bins=200,gate=temp_gate,x=subset(gate_population_key,gate_population_key$gate==temp_gate)$xlabel,y=subset(gate_population_key,gate_population_key$gate==temp_gate)$ylabel) + geom_density2d(colour = "black")+
                    geom_stats(size = 12)+
                    ggcyto_par_set(limits = list(x=c(-0.5,5),y=c(-0.5,5)))+
                    theme(
                      panel.background = element_rect(fill = "white",
                                                      colour = "black",
                                                      linewidth = 0.5, linetype = "solid"),
                      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                                      colour = "white"),
                      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                                      colour = "white"), # Increase font size for axis titles
                      axis.title = element_text(size = 30),
                      # Increase font size for axis text (numbers)
                      axis.text = element_text(size = 30)
                    ))
          }
          dev.off()
        }
      }

      for (i in 1:length(temp_ncfs)) {
        temp_plots = gate_population_key

        temp_folder = file.path(sample_outputPath,paste0("Flowplots_", pData(temp_ncfs[[i]])$pid, "_", pData(temp_ncfs[[i]])$stimulation))
        temp_folder = list.files(path=temp_folder,full.names = T)

        png(filename = paste(results_location,"/","GatingStrategy_",pData(temp_ncfs[[i]])$pid, "_", pData(temp_ncfs[[i]])$stimulation,".png",sep = ""),width = 2400,height = 1800)

        par(mfrow=c(3,3))
        gates_for_strategy = gate_population_key

        for(i_1 in gates_for_strategy$gate){
          temp.img = image_read(temp_folder[grepl(paste(pData(temp_ncfs[i])$stimulation,"_",gsub("/","_",i_1,fixed = T),".png",sep = ""),temp_folder)],depth = 16)
          par(mar=c(0,0,0,0))
          par(xpd=NA)
          plot(temp.img)

        }
        dev.off()
      }

      gs_cleanup_temp(temp_ncfs)
      setwd(results_folder_path)
      write.csv(summary_df, paste(temp_sampleBatch,"Results_DF.csv",sep = "_"),row.names = FALSE)
      write.csv(tasa_summary_df, paste(temp_sampleBatch,"tb_tasa.csv",sep = "_"),row.names = FALSE)

    }, error=function(e)  write.csv(summary_df,paste(temp_sampleBatch,"Error_in_run_Data_Not_collected.csv",sep = "_")))

  }


