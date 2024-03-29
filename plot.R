library(stringr)
library(ggplot2)
library(gggenes)
library(aplot)
library(ggeasy)
library(RColorBrewer)
library(labeling) #eventually not necessary
library(rebus)    #eventually not necessary
library(stringr)
library(shiny)
library(DT)
library(dplyr)
library(shinyWidgets)
library(shinycssloaders)
library(colourpicker)

col1 <- c("#228b22",#forestgreen
          "#186118", #dark green
          "#4f94cd", #steelblue3
          "#95bfe1", #light blue
          "#666666", #gray40
          "#dbdbdb", #gray86
          "#e12f0a", #red
          "#ed826c", #light red
          "#01e8fe", #aquamarine
          "#018b98", #dark aquamarine
          "#d601e1", #violet
          "#96019e", #dark violet
          "#fe8601", #orange
          "#be4e0f", #brown
          "#000000", #black
          "#ffffff"
) 

###program
RNAplot <- function(Data, Gff, start, end, alpha= 0.8, graph_size = 3,color=c(),
                    subgenes = FALSE, incom_genes ="_partial", max_read=NA, filter=TRUE, 
                    rename_graphs=NULL, line_visible=TRUE, arrow_body_height= 7, 
                    arrowhead_height=10, arrowhead_width=8, Gfill = TRUE, Gsize =1.2, 
                    msize= 4, Mwidth =7, ntlength =NA){
  if(is.null(alpha)) alpha <- 0.8
  if(is.null(graph_size)) graph_size <- 3
  if(is.null(subgenes)) subgenes <- TRUE
  if(is.null(line_visible)) line_visible <- TRUE
  if(is.null(incom_genes)) incom_genes <- "_partial"
  if(is.null(arrow_body_height)) arrow_body_height <- 7
  if(is.null(arrowhead_height)) arrowhead_height <- 10
  if(is.null(arrowhead_width)) arrowhead_width <- 8
  if(is.null(max_read)) max_read <- NA
  if(is.null(ntlength)) ntlength <- NA
  
  
  Data2 <- Data
  N <- str_sub(colnames(Data2),1, str_locate(colnames(Data2),or("fwd","rev"))[,1]-1)
  N1 <- unique(N[!is.na(N)])
  
  #rename
  if(length(rename_graphs)==0){rename_graphs <- N1}
  
  #filter
  q <- filter
  if(is.null(q)){q<- N}
  Data2 <- Data2[,which(N %in% q)]
  N <- N1[which(N1 %in% q)]
  rename_graphs <- rename_graphs[which(N1 %in% q)]
  
  
  #prepare data.frame for further processing
  n <- start:end
  No<- rep(NA, length(colnames(Data2))*length(n))
  D <- data.frame(Name= No, x=No, y=No)
  D$Name <- as.vector(sapply(colnames(Data2), function(i) rep(i,length(n))))
  D$x    <- rep(n,length(colnames(Data2)))
  D$y   <- unlist(Data2[n,])
  si <- 13
  if(length(Gff$color)==0){ Gff$color <- NA }
  Gff$color <- str_replace_all(Gff$color,SPC,"")
  if(length(Gff$from)==0){ Gff$from <- NA }
  if(length(Gff$to)==0){ Gff$to <- NA }
  if(length(Gff$subcolor)==0){ Gff$subcolor <- NA }
  if(length(Gff$label_type)==0){ Gff$label_type <- "bold.italic" 
  Gff$label_type[which(Gff$type!="gene")] <- "bold"}
  if(length(Gff$lane)==0){ Gff$lane <- 1 }
  
  if(!is.logical(Gff$orientation)){
    test <- c("+", "fwd", "fw", 1, "forward", "Forward", "f","TRUE")
    Gff$orientation <- Gff$orientation %in% test
  }
  
  
  
  
  
  
  
  z <- sort(sapply(unique(D$Name), function(i) sum(D$y[which(D$Name==i)])), decreasing =TRUE)
  D$Name <- factor(D$Name, levels = names(z))#sorts samples by total number of reads
  
  f <- which(str_detect(D$Name,"fwd"))              #which RNA-reads are forward
  
  Nf <- names(z)[which(str_detect(names(z),"fwd"))]
  Nf <- str_sub(Nf, 1, str_locate(Nf,"fwd")[,1]-1)  #get Names of RNA-reads, without "fwd"
  
  if(length(color)==0){
    col <- col1[1:length(Nf)]    #create colors for RNA-reads
    col <- as.character(sapply(Nf, function(i) col[which(N %in% i)]))
  }else{
    col <- as.character(sapply(Nf, function(i) color[which(N1 %in% i)]))
  }
  
  Nf1 <- rename_graphs[sapply(Nf, function(i) which(N==i))]
  
  if(!is.na(max_read)){
    D$y[D$y>max_read] <- max_read
    D$y[D$y<(-max_read)] <- -max_read
  }
  
  p1 <- ggplot(D[f,], aes(x=x, y=y))+
    theme_bw() +                                                         #deletes background
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x.bottom=element_line(size=1),
          axis.line.y.left=element_line(size=1),
          axis.text=element_text(size=si, face="bold"),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y.left = element_line(size=1),
          plot.background = element_rect(fill='transparent', color=NA))+ #cosmetic
    scale_y_continuous(expand = c(0, 0),limits = c(0,max_read))          #y-axis starts at 0
  
  
  
  if(Gfill){
    p1 <- p1 + geom_area(aes(fill=Name), alpha=alpha, 
                         position="identity",color="#000000", size=0.4)+         
      labs(fill = "")+                                                     #no Name for legend
      scale_fill_manual(values=col, labels = Nf1)
  }else{
    p1 <- p1 + geom_line(aes(color=Name), size=Gsize)+         
      labs(color = "")+                                                     #no Name for legend
      scale_color_manual(values=col, labels = Nf1)}
  
  
  
  D$Name <- factor(D$Name, levels = rev(names(z)))                      #reorder Names by size
  r <- which(str_detect(D$Name,"rev"))                                  #which RNArun is reverse
  Nr <- rev(names(z))[which(str_detect(rev(names(z)),"rev"))]
  Nr <- str_sub(Nr, 1, str_locate(Nr,"rev")[,1]-1)
  col2 <- as.character(sapply(Nr, function(i) col[which(Nf %in% i)]))   #resort colors for plot 2
  
  p2 <- ggplot(D[r,], aes(x=x, y=y)) +     #lower lineplot
    geom_hline(yintercept=0, color="#000000", size=1)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x.top=element_line(size=1),
          axis.line.y.left=element_line(size=1),
          axis.text=element_text(size=si, face="bold"),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y.left = element_line(size=1),
          plot.background = element_rect(fill='transparent', color=NA))+
    scale_x_continuous(position="top",limits = c(start,end), expand = c(0, 0))+                   #x-axis at top
    scale_y_continuous(expand = c(0, 0),limits = c(-max_read,0), labels=abs)     #y-axis starts at zero and shows absolut numbers
  
  if(Gfill){
    p2 <- p2 + geom_area(aes(fill=Name), alpha=alpha, show.legend = FALSE, 
                         position="identity",color="#000000", size=0.4)+         
      scale_fill_manual(values=col2, labels = Nr)           #reorder color
  }else{
    p2 <- p2 + geom_line(aes(color=Name), size=Gsize, show.legend = FALSE)+                      #no Name for legend
      scale_color_manual(values=col2, labels = Nr)
  }
  
  gf <- which(Gff$orientation  & ((Gff$start>start & Gff$start<end) |
                                    (Gff$end>start & Gff$end<end) |
                                    (Gff$start<start & Gff$end>start)|
                                    (Gff$start<end & Gff$end>end)))      #get all relevant Gff-elements
  gr <- which(!Gff$orientation & ((Gff$start>start & Gff$start<end) |
                                    (Gff$end>start & Gff$end<end) |
                                    (Gff$start<start & Gff$end>start)|
                                    (Gff$start<end & Gff$end>end)))
  
  lb <- arrow_body_height #width of gene-arrow
  lh <-arrowhead_height
  lw <-arrowhead_width
  
  if(length(c(gf, gr))>0){
    G <- Gff[c(gf, gr),]
    G$start <- sapply(c(gf, gr), function(i) min(Gff$start[i],Gff$end[i]))  #gives the lower value to column "start"
    G$end <- sapply(c(gf, gr), function(i) max(Gff$start[i],Gff$end[i]))    #gives the higher value to column "end"
    
    s <- which((G$start<start & !G$orientation) |
                 G$end>end & G$orientation)         #which Gff-elements would have an arrow and over the edge
    s1 <- which(G$start<start | G$start>end | 
                  G$end<start | G$end>end)             #which Gff-elements are too long
    
    if(length(s1)>0){ G$gene[s1] <- str_c(G$gene[s1], incom_genes) }  #rename elements that are to long
    
    
    G$start <- sapply(G$start, function(i) max(i, start))   #Gff-start is mininmal start of area of Interest
    G$end <- sapply(G$end, function(i) min(i, end))         #Gff-end is maxinmal end of area of Interest
    
    G$start1<-NA                  #get second category for Gff-elemnts with different arrow 
    G$end1<-NA
    G$startall <- G$start         #get category for all Gff-elements for labeling
    G$endall <- G$end
    if(length(s)>0){
      G$start1[s] <- G$start[s]   
      G$end1[s] <- G$end[s]
      G$start[s] <- NA
      G$end[s] <- NA
    }
    
    
    gf <- which(G$orientation) #which Gff-element is forward
    gr <- which(!G$orientation)#which Gff-element is reverse
  }else{G <- data.frame(matrix(NA,ncol=ncol(Gff)))
  colnames(G) <- colnames(Gff)}
  
  if(length(unique(G$lane))>1){
    g <- unique(G$lane)
    l1 <- length(unique(G$lane))
    r1 <-1:(l1*2)+nrow(G)
    gf <- c(gf,1:(l1)+nrow(G))
    gr <- c(gr,1:l1+nrow(G)+l1)
    G[r1,] <-NA
    G$orientation[r1] <- c(rep(TRUE,l1),rep(FALSE,l1))
    G$lane[r1] <- rep(g,2)
    G$molecule[r1] <- G$molecule[1]
    
  }
  if(is.na(ntlength)){
    G$gene[(G$endall -G$startall)/(end-start)*Mwidth < strwidth(G$gene, units="inches")*1.5 ] <-""
  }else{
    G$gene[abs(G$endall -G$startall) < ntlength ] <-""
  }
  
  if(length(gf)==0){     #in case there is no forward Gff-element...
    if(line_visible){
      p3 <- ggplot() +     #write empty plot
        geom_gene_arrow()+
        theme_genes()+
        scale_x_continuous(limits = c(start,end), expand = c(0, 0))+
        theme(axis.line.x.bottom=element_line(size=1),
              axis.text=element_text(size=si, face="bold"), #x-axis is plotted
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",
              plot.background = element_rect(fill='transparent', color=NA))
    }else{
      p3 <- ggplot() +     #write empty plot
        geom_gene_arrow()+
        theme_void()+
        scale_x_continuous(limits = c(start,end), expand = c(0, 0))+
        theme(axis.line.x.bottom=element_line(size=1),
              axis.text=element_text(size=si, face="bold"), #x-axis is plotted
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",
              plot.background = element_rect(fill='transparent', color=NA))
    }
    
  }else{
    col1 <- unique(G$color[gf])
    names(col1)<- unique(G$color[gf])
    p3 <- ggplot(G[gf,], aes(y = as.character(lane))) +
      geom_gene_arrow(aes(forward = orientation, fill=color,             #write Gene-map
                          xmin = as.integer(start), xmax = as.integer(end)),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(max(0,lw), "mm"), 
                      arrow_body_height = unit(lb, "mm"))+     #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start1), 
                          xmax = as.integer(end1), fill=color),
                      arrowhead_height = unit(lb*0.6, "mm"), 
                      arrowhead_width = unit(lb*0.3, "mm"), 
                      arrow_body_height = unit(lb, "mm"))+    #secound Gff-elements with ending outside area with irregular arrowheads
      scale_fill_manual(values=col1) 
    
    if(line_visible){
      p3 <- p3+ theme_genes()+
        scale_x_continuous(limits = c(start,end), expand = c(0, 0))+             #x-axis range is exactly as the lineplots
        theme(axis.line.x.bottom=element_line(size=1),         #cosmetics
              axis.text=element_text(size=si, face="bold"),
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",                          #no legend
              plot.background = element_rect(fill='transparent', color=NA))
    }else{
      p3 <- p3+ theme_void()+
        scale_x_continuous(limits = c(start,end), expand = c(0, 0))+             #x-axis range is exactly as the lineplots
        theme(axis.line.x.bottom=element_line(size=1),         #cosmetics
              axis.text=element_text(size=si, face="bold"),
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",                          #no legend
              plot.background = element_rect(fill='transparent', color=NA))
    }
    
    if(any(!is.na(G$from[gf]))){ if(subgenes){
      p3 <- p3 + geom_subgene_arrow(aes(xmin = as.integer(startall), xmax = as.integer(endall), 
                                        y = as.character(lane), xsubmin = as.integer(from),
                                        xsubmax = as.integer(to)), color="black",
                                    fill=G$subcolor[gf],
                                    arrowhead_height = unit(lb, "mm"), 
                                    arrowhead_width = unit(0, "mm"), 
                                    arrow_body_height = unit(lb, "mm"))
    }}
    p3 <- p3 +  geom_text(aes(x = (startall+endall)/2, label = gene, fontface = label_type),
                          size = msize)
    
    
  }
  
  if(length(gr)==0){                                  #in case there is no reverse Gff-element...
    p4 <- ggplot() +geom_gene_arrow()+  #write empty plot
      theme_void()+
      scale_x_continuous(limits = c(start,end))+
      theme(axis.title = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="none",
            plot.background = element_rect(fill='transparent', color=NA))+
      easy_remove_x_axis()
  }else{
    col1 <- unique(G$color[gr])
    names(col1)<- unique(G$color[gr])
    p4 <- ggplot(G[gr,], aes(y = as.character(10-as.integer(lane)))) + #write Gene-map (reverse)
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start),            
                          xmax = as.integer(end), fill = color),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(lw, "mm"), 
                      arrow_body_height = unit(lb, "mm")) +    #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start1), 
                          xmax = as.integer(end1), fill = color),
                      arrowhead_height = unit(lb*0.6, "mm"), 
                      arrowhead_width = unit(lb*0.3, "mm"), 
                      arrow_body_height = unit(lb, "mm"))+    #secound Gff-elements with ending outside area with irregular arrowheads
      scale_fill_manual(values=col1)     #colorization
    
    
    if(line_visible){
      p4 <- p4+ theme_genes()+ 
        scale_x_continuous(limits = c(start,end))+         #x-axis range is exactly as the lineplots
        theme(axis.text=element_text(size=12, face="bold"),    #text of axis
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",                          #no legend
              axis.ticks.x = element_blank(),                  
              axis.text.x = element_blank(),
              plot.background = element_rect(fill='transparent', color=NA))+
        easy_remove_x_axis()    #no x-axis
    }else{
      p4 <- p4+ theme_void()+ 
        scale_x_continuous(limits = c(start,end))+         #x-axis range is exactly as the lineplots
        theme(axis.text=element_text(size=12, face="bold"),    #text of axis
              axis.title = element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              legend.position="none",                          #no legend
              axis.ticks.x = element_blank(),                  
              axis.text.x = element_blank(),
              plot.background = element_rect(fill='transparent', color=NA))+
        easy_remove_x_axis()    #no x-axis
    }
    
    if(any(!is.na(G$from[gr]))){ if(subgenes){
      p4 <- p4 + geom_subgene_arrow(aes(xmin = as.integer(startall), xmax = as.integer(endall), 
                                        y = as.character(10-as.integer(lane)), xsubmin = as.integer(from),
                                        xsubmax = as.integer(to)), color="black",
                                    fill=G$subcolor[gr],
                                    arrowhead_height = unit(lb, "mm"), 
                                    arrowhead_width = unit(0, "mm"),
                                    arrow_body_height=unit(lb, "mm"))
    }}
    p4 <- p4+  geom_text(aes(x = (startall+endall)/2, label = gene, fontface = label_type),
                         size = msize)
  }
  
  
  p_lab <- ggplot() +               #write empty plot, with text
    annotate(geom = "text", x = 0, y = 0, label = "                             number of reads", 
             fontface = "bold", angle = 90, size=6) +
    coord_cartesian(clip = "off")+
    theme_void()                    #no background
  
  p_xlab <- ggplot() +               #write empty plot, with text
    annotate(geom = "text", x = mean(n), y = 0, label = "position", fontface = "bold", size=4) +
    coord_cartesian(clip = "off")+
    theme_void() 
  
  ap <- p2 %>% 
    insert_top(p4, height=1/graph_size) %>% 
    insert_top(p3, height = 1/graph_size) %>%
    insert_top(p1, height = 1)   #combine plots with each other, graph_size free variable
  ap <- ap %>% #insert_left(p_lab, width = 0.1)%>%  #combined complete plots with y-axis title
    insert_bottom(p_xlab, height = 0.15)
  
  
  
  return(ap) #return plot
  
}



load_Gff <- function(input){
  G <- ape::read.gff(input)
  G$ID <- str_extract(G$attributes,one_or_more(WRD)%R%DOT%R%DGT)
  G$Gen <- str_extract(G$attributes,"gene-"%R%one_or_more(WRD))
  G$Gen <- str_sub(G$Gen,6,-1)
  G$strand <- sapply(G$strand, function(i) c(1,-1)[which(c("+","-")==i)])
  b <- which(!duplicated(G[,c(4,5,7,11)]) & !is.na(G$Gen))
  G <- G[b,c(1,11,3:5,7)]
  colnames(G) <- c("molecule", "gene", "type", "start", "end", "orientation")
  G$color <- "#4682B4"
  G$color[G$type !="gene"] <- "#D6DA18"
  G$from <- NA
  G$to <- NA
  G$subcolor <- NA
  G$label_type <- "bold"
  G$label_type[G$type =="gene"] <- "bold.italic"
  G$lane <- 1
  return(G)
}




