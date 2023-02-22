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

###program
RNAplot <- function(Data, Gff, start, end, alpha= 0.8, graph_size = 3,color=c(), colorize_map=TRUE,
                    promlim=0, promsearch=c(), subgenes = FALSE, label = TRUE, 
                    incom_genes ="_partial", max_read=NA, filter=TRUE, rename_graphs=NULL, 
                    lane_manuel_selection =FALSE, line_visible=TRUE, arrow_body_height= 7, 
                    arrowhead_height=10, arrowhead_width=8){
  if(is.null(alpha)) alpha <- 0.8
  if(is.null(graph_size)) graph_size <- 3
  if(is.null(subgenes)) subgenes <- TRUE
  if(is.null(label)) label <- TRUE
  if(is.null(line_visible)) line_visible <- TRUE
  if(is.null(incom_genes)) incom_genes <- "_partial"
  if(is.null(arrow_body_height)) arrow_body_height <- 7
  if(is.null(arrowhead_height)) arrowhead_height <- 10
  if(is.null(arrowhead_width)) arrowhead_width <- 8
  if(is.null(max_read)) max_read <- NA
  
  
  
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
  
  
  z <- sort(sapply(unique(D$Name), function(i) sum(D$y[which(D$Name==i)])), decreasing =TRUE)
  D$Name <- factor(D$Name, levels = names(z))#sorts samples by total number of reads
  
  f <- which(str_detect(D$Name,"fwd"))              #which RNA-reads are forward
  
  Nf <- names(z)[which(str_detect(names(z),"fwd"))]
  Nf <- str_sub(Nf, 1, str_locate(Nf,"fwd")[,1]-1)  #get Names of RNA-reads, without "fwd"
  
  if(length(color)==0){
    col <- brewer.pal(max(length(unique(Nf)),3),"Paired")    #create colors for RNA-reads
    col <- as.character(sapply(Nf, function(i) col[which(N %in% i)]))
  }else{
    col <- as.character(sapply(Nf, function(i) color[which(N1 %in% i)]))
  }
  
  Nf1 <- rename_graphs[sapply(Nf, function(i) which(N==i))]
  
  if(!is.na(max_read)){
    D$y[D$y>max_read] <- max_read
    D$y[D$y<(-max_read)] <- -max_read
  }
  
  
  p1 <- ggplot(D[f,], aes(x=x, y=y, fill=Name)) +                        #higher lineplot
    geom_area(position="identity", alpha=alpha, size=0.4,color="#000000")+ #creates a filled lineplot
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
    scale_y_continuous(expand = c(0, 0),limits = c(0,max_read))+         #y-axis starts at 0
    labs(fill = "")+                                                     #no Name for legend
    scale_fill_manual(values=col, labels = Nf1)                          #colorization and naming of lines
  
  
  
  D$Name <- factor(D$Name, levels = rev(names(z)))                      #reorder Names by size
  r <- which(str_detect(D$Name,"rev"))                                  #which RNArun is reverse
  Nr <- rev(names(z))[which(str_detect(rev(names(z)),"rev"))]
  Nr <- str_sub(Nr, 1, str_locate(Nr,"rev")[,1]-1)
  col2 <- as.character(sapply(Nr, function(i) col[which(Nf %in% i)]))   #resort colors for plot 2
  
  
  
  p2 <- ggplot(D[r,], aes(x=x, y=y, fill=Name)) +     #lower lineplot
    geom_area(position="identity", alpha=alpha, size=0.4, show.legend = FALSE, color="#000000")+
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
    scale_y_continuous(expand = c(0, 0),limits = c(-max_read,0), labels=abs)+      #y-axis starts at zero and shows absolut numbers
    scale_fill_manual(values=col2, labels = Nr)           #reorder color
  
  gf <- which(Gff$orientation == 1  & ((Gff$start>start & Gff$start<end) |
                                         (Gff$end>start & Gff$end<end)))      #get all relevant Gff-elements
  gr <- which(Gff$orientation == -1  & ((Gff$start>start & Gff$start<end) |
                                          (Gff$end>start & Gff$end<end)))
  
  lb <- arrow_body_height #width of gene-arrow
  lh <-arrowhead_height
  lw <-arrowhead_width
  
  if(length(c(gf, gr))>0){
    G <- Gff[c(gf, gr),]
    G$start <- sapply(c(gf, gr), function(i) min(Gff$start[i],Gff$end[i]))  #gives the lower value to column "start"
    G$end <- sapply(c(gf, gr), function(i) max(Gff$start[i],Gff$end[i]))    #gives the higher value to column "end"
    
    s <- which((G$start<start & G$orientation==-1) |
                 G$end>end & G$orientation==1)         #which Gff-elements would have an arrow and over the edge
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
    
    if(!colorize_map){ G$color <- "#8dd3c7"}     #if colorize_map is FALSE, then give standard color
    
    gf <- which(G$orientation==1) #which Gff-element is forward
    gr <- which(G$orientation==-1)#which Gff-element is reverse
  }
  
  if(promlim>0){
    li <- promlim
    ln <- length(n)/40
    if(length(promsearch)==0){
      f <- which(str_detect(colnames(Data2),"fwd"))
      r <- which(str_detect(colnames(Data2),"rev"))
    }else{
      f <- which(colnames(Data2) %in% str_c(promsearch,"fwd"))
      r <- which(colnames(Data2) %in% str_c(promsearch,"rev"))
    }
    
    po <- sapply(n, function(i) any(Data2[i,f]>li))
    posf <- n[which(po==FALSE & po!=c(po[2:length(po)],FALSE))+1]
    
    po <- sapply(n, function(i) any(Data2[i,r]< -(li)))
    posr <- n[which(po==TRUE & po!=c(po[2:length(po)],TRUE))]
    Pf <- data.frame(x = posf, x2 = posf+ln, y1 = rep(1, length(posf)), y2 = rep(1.4, length(posf)))
    Pr <- data.frame(x = posr, x2 = posr-ln, y1 = rep(1, length(posr)), y2 = rep(0.6, length(posr)))
    
  }
  
  if(lane_manuel_selection & length(c(gf, gr))>0){
    G$lane <- 1
    for(i in 1:nrow(G)){
      G$lane[i] <- as.integer(select.list(as.character(1:9), title=str_c("Which lane is ", G$gene[i],"?"), 
                                          graphics = T,multiple = FALSE))
    }
  }
  if(length(unique(G$lane))>1){
    g <- unique(G$lane)
    l1 <- length(unique(G$lane))
    r1 <-1:(l1*2)+nrow(G)
    gf <- c(gf,1:(l1)+nrow(G))
    gr <- c(gr,1:l1+nrow(G)+l1)
    G[r1,] <-NA
    G$orientation[r1] <- c(rep(1,l1),rep(-1,l1))
    G$lane[r1] <- rep(g,2)
    G$molecule[r1] <- G$molecule[1]
    
  }
  
  
  
  if(length(gf)==0){     #in case there is no forward Gff-element...
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
  }else{
    col1 <- unique(G$color[gf])
    names(col1)<- unique(G$color[gf])
    p3 <- ggplot(G[gf,], aes(y = as.character(lane))) +
      geom_gene_arrow(aes(forward = as.integer(orientation), fill=color,             #write Gene-map
                          xmin = as.integer(start), xmax = as.integer(end)),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(max(0,lw), "mm"), 
                      arrow_body_height = unit(lb, "mm"))+     #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = as.integer(orientation), xmin = as.integer(start1), 
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
    if(label==TRUE){
      p3 <- p3 +  geom_gene_label(aes(xmin = as.integer(startall),       #labeling of arrows
                                      xmax = as.integer(endall),label = gene,
                                      fontface = label_type), grow = TRUE, 
                                  padding.x = grid::unit(0, "mm"),
                                  padding.y = grid::unit(0, "lines"))
      
    }
  }
  
  if(promlim > 0){
    if(nrow(Pf)>0){
      p3 <- p3 + geom_segment(aes(x = x, y = y1, xend = x, yend = y2, colour = "segment"), data = Pf, color="#000000",
                              lineend = "round", size=1)
      p3 <- p3 +geom_segment(aes(x = x, y = y2, xend = x2, yend = y2), data = Pf,
                             arrow = arrow(length = unit(0.2, "cm"), type = "open"),
                             lineend = "round", size=1)
    }}
  
  
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
    p4 <- ggplot(G[gr,], aes(y = as.character(10-lane))) + #write Gene-map (reverse)
      geom_gene_arrow(aes(forward = as.integer(orientation), xmin = as.integer(start),            
                          xmax = as.integer(end), fill = color),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(lw, "mm"), 
                      arrow_body_height = unit(lb, "mm")) +    #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = as.integer(orientation), xmin = as.integer(start1), 
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
                                        y = as.character(10-lane), xsubmin = as.integer(from),
                                        xsubmax = as.integer(to)), color="black",
                                    fill=G$subcolor[gr],
                                    arrowhead_height = unit(lb, "mm"), 
                                    arrowhead_width = unit(0, "mm"),
                                    arrow_body_height=unit(lb, "mm"))
    }}
    if(label==TRUE){
      p4 <- p4 +  geom_gene_label(aes(xmin = as.integer(startall),       #labeling of arrows
                                      xmax = as.integer(endall),label = gene,
                                      fontface = label_type), grow = TRUE, 
                                  padding.x = grid::unit(0, "mm"),
                                  padding.y = grid::unit(0, "lines"))
      
    }
    
  }
  
  if(promlim > 0){
    if(nrow(Pr)>0){
      p4 <- p4 + geom_segment(aes(x = x, y = y1, xend = x, yend = y2, colour = "segment"), data = Pr, color="#000000",
                              lineend = "round", size=1)
      p4 <- p4 +geom_segment(aes(x = x, y = y2, xend = x2, yend = y2), data = Pr,
                             arrow = arrow(length = unit(0.2, "cm"), type = "open"),
                             lineend = "round", size=1)
    }}
  
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







