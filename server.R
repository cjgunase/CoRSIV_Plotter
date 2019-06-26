library(shiny)
library(corrplot)
CoRSIV <- read.csv("./CoRSIV_9926_gencode28.csv")
function(input, output,session) {
  
  output$gene_symbol = renderText({
    paste0("Option2: Details of CoRSIVs in +/- 1Mb of  ", toupper(input$gene_sym))
  })
  
  output$coodinate = renderText({
    paste0("Option3: Details of CoRSIVs in +/- 1Mb of  ", input$ucsc_coord)
  })
  
  output$help <- renderText({
    paste0("--------------------------------------------------------------------------------")
  })

  
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("chr",input$select, "zip", sep=".")
    },

    content <- function(file) {
      file.copy(paste("chr",input$select, "zip", sep="."), file)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename <- function() {
      x = input$ucsc_coord
      paste("CoRSIVs flanking 1Mb of",x,".pdf",sep = "")
    },
    
    content <- function(file) {
      ###################################################################
      x = input$ucsc_coord
      #x = "1:1068100"
      chr_num <- as.numeric(strsplit(x = x,split=":")[[1]][1])
      mid <- as.numeric( as.numeric(gsub(",", "",strsplit(x = x,split=":")[[1]][2])))
      start <- mid-1000000
      end <- mid + 1000000
      CoRSIV_chr1 <- CoRSIV[CoRSIV$Chromosome==chr_num & (start < CoRSIV$BP) & (CoRSIV$BP < end),]
      validate(
        need(dim(CoRSIV_chr1)[1]!=0, "Incorrect region or No CoRSIVs found in that region")
      )
      
      
      annotated_blocks <- CoRSIV_chr1[order(CoRSIV_chr1$BP),]
      candidate_blocks <-annotated_blocks
      
      
      
      min_corr = 0
      bin_num = 1
      pos<-c()
      cor_val <-c()
      start=bin_num
      end = dim(candidate_blocks)[1]
      
      #start = 1
      #end = 207
      num_of_bins = dim(candidate_blocks)[1]
      test <- t(candidate_blocks[start:end,5:14])
      colnames(test) <-candidate_blocks$Bin.Name[start:end]
      col_names<-as.character(candidate_blocks$Bin.Name[start:end])
      M <-matrix(0,dim(test)[2],dim(test)[2])
      for(k in 1:(dim(test)[2]-1)){
        if(k>40){# increase this to plot columns in figure
          #M[i1,i2] <-0
          break
        }
        for(i in 1:(dim(test)[2]-k)){
          i1 <-i
          i2 <-i+k
          M[i1,i2] <- cor(test[,i1],test[,i2],use="pairwise.complete.obs")
          if(k==1){
            if(M[i1,i2] >= min_corr){
              #cor_pair_df<-rbind(cor_pair_df,c(as.character(col_names[i1]),as.character(col_names[i2]),as.numeric(M[i1,i2])))
              pos <- c(pos,as.character(col_names[i1]))
              cor_val <-c(cor_val,as.numeric(M[i1,i2]))
              pos <- c(pos,as.character(col_names[i2]))
              cor_val <-c(cor_val,as.numeric(M[i1,i2]))
              
            }
          }
          
        }
      }
      m <- M
      
      #paste(candidate_blocks$EMBL_ID,candidate_blocks$Bin.Name,sep = "__")
      
      colnames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      rownames(m)<-as.character(paste(paste(annotated_blocks$USCS_Coordinates_CoRSIV,annotated_blocks$GS1.Name,annotated_blocks$GS1.Distance,
                                            annotated_blocks$GS2.Name,annotated_blocks$GS2.Distance,
                                            annotated_blocks$GS3.Name,annotated_blocks$GS3.Distance,
                                            annotated_blocks$GE1.Name,annotated_blocks$GE1.Distance,
                                            annotated_blocks$GE2.Name,annotated_blocks$GE2.Distance,
                                            annotated_blocks$promoter....3kb.TSS..overlapping.genes,
                                            annotated_blocks$gene.body.overlapping.genes,
                                            annotated_blocks$three.prime.region.....3kb.TES..overlapping.genes,sep = " | "),candidate_blocks$Uniq_ID,candidate_blocks$Bin.Name,sep = " ------> "))[start:end]
      
      #colnames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      #rownames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      
      #install.packages("extrafont")
      #library("extrafont")
      #font_import()
      #library(extrafont)
      #loadfonts()
      if(num_of_bins<50){
        size = 50
      }else if(50 < num_of_bins & num_of_bins < 75){
        size = num_of_bins
      }else if(num_of_bins > 100){
        size = 100
      }
      
      pdf(file=paste("./CoRSIVs_map_of_1Mb_window_cerntered_on_given_position.pdf",sep = ""),family = "Courier",width = size,height = size)
      corrplot(m, type = "upper", tl.pos = "td",
               method = "square", tl.cex = 1, tl.col = 'black',
               diag = FALSE,tl.srt=45,addgrid.col = NA,is.corr = FALSE,addCoef.col = "black",number.cex = .25,cl.lim = c(-1, 1))
      
      dev.off()
      
      
      ##################################################################
      file.copy("CoRSIVs_map_of_1Mb_window_cerntered_on_given_position.pdf", file)
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename <- function() {
      x = toupper(input$gene_sym)
      paste("CoRSIVs_Overlapping_",x,"_and_CoRSIVs_in_+/-1Mb_region.pdf",sep = "")
    },
    
    content <- function(file) {
      
      #################################################################
      #CoRSIV <- read.csv("./CoRSIV_9926_gencode28.csv")
      
      x = toupper(input$gene_sym)
      
      CoRSIV_gene <- CoRSIV[(CoRSIV$promoter....3kb.TSS..overlapping.genes==x)|(CoRSIV$gene.body.overlapping.genes==x)|(CoRSIV$three.prime.region.....3kb.TES..overlapping.genes==x),]

      validate(
        need(dim(CoRSIV_gene)[1]!=0 & dim(CoRSIV_gene)[1] < 100, paste("CoRSIVs were found overlapping",x,"or incorrect gene symbol !",sep = "_"))
        
      )
      
      chr_num <- CoRSIV_gene$Chromosome[1]
      mid <- CoRSIV_gene$BP[1]
      
      start <- mid-1000000
      end <- mid + 1000000
      
      CoRSIV_chr <- CoRSIV[CoRSIV$Chromosome==chr_num & (start < CoRSIV$BP) & (CoRSIV$BP < end),]
      
      
      
      annotated_blocks <- CoRSIV_chr[order(CoRSIV_chr$BP),]
      candidate_blocks <-annotated_blocks
      
      
      
      
      min_corr = 0
      library(corrplot)
      bin_num = 1
      pos<-c()
      cor_val <-c()
      start=bin_num
      end = dim(candidate_blocks)[1]
      
      #start = 1
      #end = 207
      num_of_bins = dim(candidate_blocks)[1]
      test <- t(candidate_blocks[start:end,5:14])
      colnames(test) <-candidate_blocks$Bin.Name[start:end]
      col_names<-as.character(candidate_blocks$Bin.Name[start:end])
      M <-matrix(0,dim(test)[2],dim(test)[2])
      for(k in 1:(dim(test)[2]-1)){
        if(k>40){# increase this to plot columns in figure
          #M[i1,i2] <-0
          break
        }
        for(i in 1:(dim(test)[2]-k)){
          i1 <-i
          i2 <-i+k
          M[i1,i2] <- cor(test[,i1],test[,i2],use="pairwise.complete.obs")
          if(k==1){
            if(M[i1,i2] >= min_corr){
              #cor_pair_df<-rbind(cor_pair_df,c(as.character(col_names[i1]),as.character(col_names[i2]),as.numeric(M[i1,i2])))
              pos <- c(pos,as.character(col_names[i1]))
              cor_val <-c(cor_val,as.numeric(M[i1,i2]))
              pos <- c(pos,as.character(col_names[i2]))
              cor_val <-c(cor_val,as.numeric(M[i1,i2]))
              
            }
          }
          
        }
      }
      m <- M
      
      #paste(candidate_blocks$EMBL_ID,candidate_blocks$Bin.Name,sep = "__")
      
      colnames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      rownames(m)<-as.character(paste(paste(annotated_blocks$USCS_Coordinates_CoRSIV,annotated_blocks$GS1.Name,annotated_blocks$GS1.Distance,
                                            annotated_blocks$GS2.Name,annotated_blocks$GS2.Distance,
                                            annotated_blocks$GS3.Name,annotated_blocks$GS3.Distance,
                                            annotated_blocks$GE1.Name,annotated_blocks$GE1.Distance,
                                            annotated_blocks$GE2.Name,annotated_blocks$GE2.Distance,
                                            annotated_blocks$promoter....3kb.TSS..overlapping.genes,
                                            annotated_blocks$gene.body.overlapping.genes,
                                            annotated_blocks$three.prime.region.....3kb.TES..overlapping.genes,sep = " | "),candidate_blocks$Uniq_ID,candidate_blocks$Bin.Name,sep = " ------> "))[start:end]
      
      #colnames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      #rownames(m)<-as.character(candidate_blocks$Bin.Name)[start:end]
      
      #install.packages("extrafont")
      #library("extrafont")
      #font_import()
      #library(extrafont)
      #loadfonts()
      if(num_of_bins<50){
        size = 50
      }else if(50 < num_of_bins & num_of_bins < 100){
        size = num_of_bins
      }else if(num_of_bins > 100){
        size = 100
      }
      
      pdf(file=paste("CoRSIVs_Overlapping_",x,"_and_CoRSIVs_in_1Mb_region.pdf",sep = ""),family = "Courier",width = size,height = size)
      corrplot(m, type = "upper", tl.pos = "td",
               method = "square", tl.cex = 1, tl.col = 'black',
               diag = FALSE,tl.srt=45,addgrid.col = NA,is.corr = FALSE,addCoef.col = "black",number.cex = .25,cl.lim = c(-1, 1))
      dev.off()
      
      ###############################################################################

      
      
      file.copy(paste("CoRSIVs_Overlapping_",x,"_and_CoRSIVs_in_1Mb_region.pdf",sep = ""), file)

    }

  )
  

  
  output$coolplot2 <-  renderTable({
    #CoRSIV <- read.csv("./CoRSIV_9926_gencode28.csv")
    
    x = toupper(input$gene_sym)
    
    
    CoRSIV_gene <- CoRSIV[(CoRSIV$promoter....3kb.TSS..overlapping.genes==x)|(CoRSIV$gene.body.overlapping.genes==x)|(CoRSIV$three.prime.region.....3kb.TES..overlapping.genes==x),]
    cat(file=stderr(), "dimension 1",dim(CoRSIV_gene)[1], "\n")
    
    validate(
      need(dim(CoRSIV_gene)[1]!=0 & dim(CoRSIV_gene)[1] < 100, paste("No CoRSIVs were found within +/- 1 Mb of the",x,"gene or incorrect gene symbol !",sep = " "))
      
    )
    
    chr_num <- CoRSIV_gene$Chromosome[1]
    mid <- CoRSIV_gene$BP[1]
    
    start <- mid-1000000
    end <- mid + 1000000
    
    CoRSIV_chr <- CoRSIV[CoRSIV$Chromosome==chr_num & (start < CoRSIV$BP) & (CoRSIV$BP < end),]
    

    
    annotated_blocks <- CoRSIV_chr[order(CoRSIV_chr$BP),]
    candidate_blocks <-annotated_blocks

    
    
    temp2 <- CoRSIV_gene[c("Uniq_ID","CorSIV_ID","genomic_range_of_CorSIV","CoRSIV_IIR","min.ITC","USCS_Coordinates_CoRSIV","promoter....3kb.TSS..overlapping.genes"           
                 ,"gene.body.overlapping.genes"                      
                 ,"three.prime.region.....3kb.TES..overlapping.genes")] 
    temp2 <- temp2[!duplicated(temp2$Uniq_ID),]
    temp2 <- temp2[order(temp2$CorSIV_ID),]
    temp2

    
    #plot(rnorm(100))
  })
  
  
  output$coolplot <-  renderTable({
    #CoRSIV <- read.csv("./CoRSIV_9926_gencode28.csv")
    x = input$ucsc_coord
    #x = "1:1068100"
    
    chr_num <- as.numeric(strsplit(x = x,split=":")[[1]][1])
    mid <- as.numeric( as.numeric(gsub(",", "",strsplit(x = x,split=":")[[1]][2])))
    
    start <- mid-1000000
    end <- mid + 1000000
    
    CoRSIV_chr1 <- CoRSIV[CoRSIV$Chromosome==chr_num & (start < CoRSIV$BP) & (CoRSIV$BP < end),]
    if(is.na(mid)){CoRSIV_chr1<-data.frame()}
    validate(
      need(dim(CoRSIV_chr1)[1]!=0, paste("No CoRSIVs were found within +/- 1Mb of the",x, "genomic coordinate or incorrect input format !",sep = " "))
    )
    
    
    annotated_blocks <- CoRSIV_chr1[order(CoRSIV_chr1$BP),]
    candidate_blocks <-annotated_blocks
    
    
    
    temp1 <- candidate_blocks[c("Uniq_ID","CorSIV_ID","genomic_range_of_CorSIV","CoRSIV_IIR","min.ITC","USCS_Coordinates_CoRSIV","promoter....3kb.TSS..overlapping.genes"           
                                ,"gene.body.overlapping.genes"                      
                                ,"three.prime.region.....3kb.TES..overlapping.genes")]                                      
    temp1 <- temp1[!duplicated(temp1$Uniq_ID),]
    temp1 <- temp1[order(temp1$CorSIV_ID),]
    temp1
    
    #plot(rnorm(100))
  })
  
}