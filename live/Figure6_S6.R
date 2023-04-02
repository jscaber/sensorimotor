library("tidyverse")
library("patchwork")
library("ggpubr")
library("cowplot")
library(Cairo)

setwd("set here")

file_list <- list.files(path = "set here", pattern = "_tracked.txt$", recursive = TRUE)

direction <- list()
max_speed<- list()
type <-  list()
direction <-  list()
max_speed <-  list()
average_speed <-  list()
pause_time <-  list()
pause_percentage <-  list()
total_frames <-  list()
distance_traveled <-  list()
plane <- list()
f_name <- list()
p_name <- list()

j <- 1

for(fn in file_list){
  filename <- paste0("set here ",fn)
  if(file.info(filename)$size == 0) { next } 
  lines <- readLines(filename)

  tables <- split(lines, cumsum( grepl("^#", lines)))
  names(tables) <- vapply(tables, function(table)sub("^#", "", table[1]), "")
  data <- lapply(tables, function(text)read.table(text=text, header=FALSE))
  data <- lapply(data, setNames, nm = c("date2","x","y"))
  
  # This code recognises particles that are not horizonal
  # If these are present, the imaging almost certainly
  euclidian <- FALSE
  for(i in 1:length(data)) {
    if(diff(range(data[[i]][,"y"])) > 100){
      euclidian <- TRUE
      break
    }
  }
  
  for(i in 1:length(data)) {
    
    # February was imaged the other way around
    if(grepl("Feb",fn)){
      data[[i]]$x <- rev(data[[i]]$x)
      data[[i]]$y <- rev(data[[i]]$y)
    }
    
    position <- data[[i]]$x
    if(euclidian == TRUE){
      diffs <- sqrt(diff(data[[i]]$x)^2 + diff(data[[i]]$y)^2)
      diffs <- sign(diff(data[[i]]$x)) * diffs
      position <- diffinv(diff(diffs), xi = 200)
      plane[[j]] <- "euclidian"
    } else {
      plane[[j]] <- "linear"
    }
    
    total_frames[[j]] <- length(position)
    distance_traveled[[j]] <- position[length(position)] - position[1]
    f_name[[j]] <- fn
    p_name[[j]] <- names(data[i])
    

    
    if(abs(position[1] - position[length(position)]) < length(position)*3){
      type[[j]] <- "Dynamic Stationary"
      direction[[j]] <- NA
      max_speed[[j]] <- NA
      average_speed[[j]] <- NA
      pause_time[[j]] <- length(position)-1
    } else {
      type[[j]] <- "Moving"
      if(position[1] - position[length(position)] <0){
        direction[[j]] <- "Right"
      } else {
       direction[[j]] <- "Left"
      }
      max_speed[[j]] <- diff(position,2)[which.max(abs(diff(position,2)/2))]/2
      speed <- diff(position)
      pause_time[[j]] <- sum(abs(speed) < 5)
      average_speed[[j]] <- mean(speed[abs(speed)>4])
    }
    j = j + 1
  }
}

results <- c(f_name, p_name, plane, type, average_speed,direction,distance_traveled, pause_time, total_frames, max_speed)
for (item in  results){
  item <- unlist(item, recursive = TRUE)
}

df <- cbind(f_name, p_name, plane, type, average_speed,direction,distance_traveled, pause_time, total_frames, max_speed)

df <- as.data.frame(df)
df$average_speed <- as.numeric(df$average_speed)
df$p_name <- as.character(df$p_name)
df$plane <- as.character(df$plane)
df$type <- as.character(df$type)
df$direction <- as.character(df$direction)
df$distance_traveled <- as.numeric(df$distance_traveled)
df$pause_time <- as.numeric(df$pause_time)
df$total_frames <- as.numeric(df$total_frames)
df$max_speed <- as.numeric(df$max_speed)
df <- df %>%
  separate(f_name, c("differentiation", "target", "line", "area"), "/")
df2 <- df[!grepl("Body|Axon",df$area),]


summarise_axonal_transport<- function(df){
  df$left = df$list_total_particles-df$list_stationary-df$list_right
  df$prop_right = df$list_right/(df$list_total_particles-df$list_stationary)
  df$prop_stationary = df$list_stationary/df$list_total_particles
  df$distance_left_pp = df$list_distance_left/df$left
  df$distance_right_pp = df$list_distance_right/df$list_right
  df$distance_total_pp = df$list_distance_tot/(df$list_total_particles-df$list_stationary)
  df$pauses_moving_particles = df$list_pause_time_moving/df$list_total_frames_moving
  df$pauses_all_particles = df$list_pause_time_total/df$list_total_frames
  return(df)
}

summarise_axonal_transport_raw<- function(df){
  df$pauses_moving_particles = df$pause_time_moving/df$total_frames_moving
  df$pauses_all_particles = df$pause_time_total/df$total_frames
  return(df)
}

differentiation <- list()
target <- list()
line <- list()
area <-list()
list_plane <-list()
list_stationary <- list()
list_right <- list()
list_distance_tot <- list()
list_distance_right <- list()
list_distance_left<- list()
list_distance_stat <- list()
list_speed_tot <- list()
list_speed_right <- list()
list_speed_left<- list()
list_maxspeed_tot <- list()
list_maxspeed_right <- list()
list_maxspeed_left<- list()
list_avmaxspeed_tot <- list()
list_avmaxspeed_right <-list()
list_avmaxspeed_left<- list()
list_pause_time_moving <- list()
list_pause_time_total <- list()
list_total_frames_moving <- list()
list_total_frames <-list()
list_total_particles <- list()

k=1

for(uniq_diff in unique(df2$differentiation)){
  tmp_diff <-  filter(df2, differentiation == uniq_diff)
  for(uniq_target in unique(tmp_diff$target)) {
    tmp_target <-  filter(tmp_diff, target == uniq_target)
    for(uniq_line in unique(tmp_target$line)) {
      tmp_line <- filter(tmp_target, line == uniq_line)
      for (uniq_area in unique(tmp_line$area)) {
        tmp_area <- filter(tmp_line, area == uniq_area)
        ds <- tmp_area[tmp_area$type == "Dynamic Stationary",]
        mov <- tmp_area[tmp_area$type == "Moving",]
        mov_r <- mov[mov$direction == "Right",]
        mov_l <- mov[mov$direction == "Left",]
        differentiation[[k]] <- uniq_diff
        target[[k]] <- uniq_target
        line[[k]] <- uniq_line
        area[[k]] <- uniq_area
        list_plane[[k]] <- tmp_area[[1,"plane"]]
        list_stationary[[k]] <- length(ds$type)
        list_right[[k]] <- length(mov_r$type)
        list_distance_tot[[k]] <- sum(abs(mov$distance_traveled))
        list_distance_right[[k]] <- sum(mov_r$distance_traveled)
        list_distance_left[[k]]<- sum(mov_l$distance_traveled)
        list_distance_stat[[k]] <- sum(abs(ds$distance_traveled))
        list_speed_tot[[k]] <- mean(abs(mov$average_speed))
        list_speed_right[[k]] <- mean(mov_r$average_speed)
        list_speed_left[[k]]<- mean(mov_l$average_speed)
        if(dim(mov)[1] != 0){
          list_maxspeed_tot[[k]] <- mov[which.max(abs(mov$max_speed)),]$max_speed
        }else { list_maxspeed_tot[[k]] <- NA}
        if(dim(mov_r)[1] != 0){
          list_maxspeed_right[[k]] <- mov_r[which.max(abs(mov_r$max_speed)),]$max_speed
        }else { list_maxspeed_right[[k]] <- NA}
        if(dim(mov_l)[1] != 0){
          list_maxspeed_left[[k]] <- mov_l[which.max(abs(mov_l$max_speed)),]$max_speed
        }else { list_maxspeed_left[[k]] <- NA}
        list_avmaxspeed_tot[[k]] <- mean(abs(mov$max_speed))
        list_avmaxspeed_right[[k]] <- mean(mov_r$max_speed)
        list_avmaxspeed_left[[k]]<- mean(mov_l$max_speed)  
        list_pause_time_moving[[k]]<- sum(mov$pause_time)
        list_pause_time_total[[k]] <- sum(tmp_area$pause_time)
        list_total_frames_moving[[k]] <- sum(mov$total_frames)
        list_total_frames[[k]] <- sum(tmp_area$total_frames)
        list_total_particles[[k]] <- length(tmp_area$p_name)
        k = k +1
      }
    }
  }
}

summary_areas <- cbind(differentiation,target,line,area,list_plane,list_stationary,list_right,list_distance_tot,list_distance_right,list_distance_left,list_distance_stat,list_speed_tot,list_speed_right,list_speed_left,list_maxspeed_tot,list_maxspeed_right,list_maxspeed_left,list_avmaxspeed_tot,list_avmaxspeed_right,list_avmaxspeed_left,list_pause_time_moving,list_pause_time_total,list_total_frames_moving,list_total_frames,list_total_particles)
summary_areas <- as.data.frame(summary_areas)
summary_areas[,6:25] = apply(summary_areas[,6:25], 2, function(x) as.numeric(x))
summary_areas[,1:5] = apply(summary_areas[,1:5], 2, function(x) as.character(x))

summary_areas <- summary_areas[summary_areas$list_plane == "linear", ]
# join the two, keeping all of df1's indices
df3 <- dplyr::semi_join(df2, summary_areas, by=c('differentiation', 'target', "line", "area"))

summary_areas <- summarise_axonal_transport(summary_areas)




differentiation <- list()
target <- list()
line <- list()
area <-list()
list_plane <-list()
list_stationary <- list()
list_right <- list()
list_distance_tot <- list()
list_distance_right <- list()
list_distance_left<- list()
list_distance_stat <- list()
list_speed_tot <- list()
list_speed_right <- list()
list_speed_left<- list()
list_maxspeed_tot <- list()
list_maxspeed_right <- list()
list_maxspeed_left<- list()
list_avmaxspeed_tot <- list()
list_avmaxspeed_right <-list()
list_avmaxspeed_left<- list()
list_pause_time_moving <- list()
list_pause_time_total <- list()
list_total_frames_moving <- list()
list_total_frames <-list()
list_total_particles <- list()

k=1

for(uniq_diff in unique(df3$differentiation)){
  tmp_diff <-  filter(df3, differentiation == uniq_diff)
  for(uniq_target in unique(tmp_diff$target)) {
    tmp_target <-  filter(tmp_diff, target == uniq_target)
    for(uniq_line in unique(tmp_target$line)) {
      tmp_line <- filter(tmp_target, line == uniq_line)
      ds <- tmp_line[tmp_line$type == "Dynamic Stationary",]
      mov <- tmp_line[tmp_line$type == "Moving",]
      mov_r <- mov[mov$direction == "Right",]
      mov_l <- mov[mov$direction == "Left",]
      differentiation[[k]] <- uniq_diff
      target[[k]] <- uniq_target
      line[[k]] <- uniq_line
      area[[k]] <- uniq_area
      list_plane[[k]] <- tmp_line[[1,"plane"]]
      list_stationary[[k]] <- length(ds$type)
      list_right[[k]] <- length(mov_r$type)
      list_distance_tot[[k]] <- sum(abs(mov$distance_traveled))
      list_distance_right[[k]] <- sum(mov_r$distance_traveled)
      list_distance_left[[k]]<- sum(mov_l$distance_traveled)
      list_distance_stat[[k]] <- sum(abs(ds$distance_traveled))
      list_speed_tot[[k]] <- mean(abs(mov$average_speed))
      list_speed_right[[k]] <- mean(mov_r$average_speed)
      list_speed_left[[k]]<- mean(mov_l$average_speed)
      if(dim(mov)[1] != 0){
        list_maxspeed_tot[[k]] <- mov[which.max(abs(mov$max_speed)),]$max_speed
      }else { list_maxspeed_tot[[k]] <- NA}
      if(dim(mov_r)[1] != 0){
        list_maxspeed_right[[k]] <- mov_r[which.max(abs(mov_r$max_speed)),]$max_speed
      }else { list_maxspeed_right[[k]] <- NA}
      if(dim(mov_l)[1] != 0){
        list_maxspeed_left[[k]] <- mov_l[which.max(abs(mov_l$max_speed)),]$max_speed
      }else { list_maxspeed_left[[k]] <- NA}
      list_avmaxspeed_tot[[k]] <- mean(abs(mov$max_speed))
      list_avmaxspeed_right[[k]] <- mean(mov_r$max_speed)
      list_avmaxspeed_left[[k]]<- mean(mov_l$max_speed)  
      list_pause_time_moving[[k]]<- sum(mov$pause_time)
      list_pause_time_total[[k]] <- sum(tmp_line$pause_time)
      list_total_frames_moving[[k]] <- sum(mov$total_frames)
      list_total_frames[[k]] <- sum(tmp_line$total_frames)
      list_total_particles[[k]] <- length(tmp_line$p_name)
      k = k +1
    }
  }
}

summary_lines <- cbind(differentiation,target,line,area,list_plane,list_stationary,list_right,list_distance_tot,list_distance_right,list_distance_left,list_distance_stat,list_speed_tot,list_speed_right,list_speed_left,list_maxspeed_tot,list_maxspeed_right,list_maxspeed_left,list_avmaxspeed_tot,list_avmaxspeed_right,list_avmaxspeed_left,list_pause_time_moving,list_pause_time_total,list_total_frames_moving,list_total_frames,list_total_particles)
summary_lines <- as.data.frame(summary_lines)
summary_lines[,6:25] = apply(summary_lines[,6:25], 2, function(x) as.numeric(x))
summary_lines[,1:5] = apply(summary_lines[,1:5], 2, function(x) as.character(x))

summary_lines <- summary_lines[summary_lines$list_total_particles >= 50, ]
# join the two, keeping all of df1's indices
df3 <- dplyr::semi_join(df3, summary_lines, by=c('differentiation', 'target', "line"))


summary_lines <- summarise_axonal_transport(summary_lines)



exploratory <- function(df, stationary_threshold = 1, title, string){
  df$group <- "PT-CRISPR (C22)"
  if(sum(grepl("C2|C92",df$line))!=0){
    df[grep("C2|C92",df$line),]$group <- "PT-C9-2"
  }
  if(sum(grepl("C1|C91",df$line))!=0){
    df[grep("C1|C91",df$line),]$group <- "PT-C9-1"
  }
  if(sum(grepl("C4|C94",df$line))!=0){
      df[grep("C4|C94",df$line),]$group <- "PT-C9-4"
  }
  if(sum(grepl("OX",df$line))!=0){
    df[grep("OX",df$line),]$group <- "CTR-OX3"
  }
  if(sum(grepl("180",df$line))!=0){
    df[grep("180",df$line),]$group <- "CTR-180"
  }
  if(sum(grepl("840",df$line))!=0){
    df[grep("840",df$line),]$group <- "CTR-840"
  }
  if(sum(grepl("841",df$line))!=0){
    df[grep("841",df$line),]$group <- "CTR-841"
  }
  if(sum(grepl("856",df$line))!=0){
    df[grep("856|853",df$line),]$group <- "CTR-856"
  }
  
  df$col <- "ISO"
  df[grep("C",df$line),]$col <- "C9"
  df[grep("OX|85|84|18",df$line),]$col <- "CTR"
  
  # Remove dead cells
  df2 <- df[df$prop_stationary < stationary_threshold,]
  df <- df[(df$col != "ISO"),]
  df2 <- df2[(df2$col != "ISO"),]
  
  df$col <- factor(df$col, levels=c("CTR", "C9"))
  df2$col <- factor(df2$col, levels=c("CTR", "C9"))
  
  
  df2$distance_total_pp = df2$distance_total_pp * 0.093
  df2$list_avmaxspeed_tot = df2$list_avmaxspeed_tot * 0.0465
  df2$list_maxspeed_tot = df2$list_maxspeed_tot * 0.0465
  df2$list_speed_tot = df2$list_speed_tot * 0.0465
  
  anterograde <- df2[df2$direction == "Left",]
  retrograde <- df2[df2$direction == "Right",]
  total <- df2
  my_list <- list(anterograde, retrograde,total)
  names(my_list) <- c("anterograde", "retrograde","total")
  
  p_retrograde <- list()
  p_anterograde <- list()
  p_total <- list()
  graphs <- list(p_anterograde, p_retrograde,p_total)  
  names(graphs)  <- c("anterograde", "retrograde","total")
  
  if(grepl("sn",string, ignore.case = TRUE)){
    colour = c("#8E98C0", "dodgerblue")
  }
  else{
    colour = c("#C09F98", "indianred")
  }

  for(j in names(my_list)){
    i = my_list[[j]]
    p_stationary <- ggplot(df, aes(x=group,y=prop_stationary, fill=col)) + geom_boxplot() +   theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Proportion", title =paste0("Proportion of stationary particles in\n", title))
    p_total <-ggplot(i, aes(x=group,y=list_total_particles, fill=col)) + geom_boxplot() +  theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Total particles", title =paste0("Total particles per experiment in\n", title))
    p_retro <-ggplot(i, aes(x=group,y=prop_right, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Proportion", title =paste0("Retrograde transport in\n", title))
    p_pauseall <- ggplot(i, aes(x=group,y=pauses_all_particles, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Proportion of time in pause/all particles", title =paste0("Pause time (all particles) in\n", title))
    p_pausemoving  <- ggplot(i, aes(x=group,y=pauses_moving_particles, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Proportion of time in pause/moving particles", title =paste0("Pause time (moving particles) in\n", title))
    p_meandist  <- ggplot(i, aes(x=group,y=distance_total_pp, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Mean distance travelled (\u03bCm)", title =paste0("Distance travelled in\n", title))
    p_meanmaxspeed  <- ggplot(i, aes(x=group,y=list_avmaxspeed_tot, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Mean maximum speed (\u03bCm/s)", title =paste0("Mean maximum speed in\n", title))
    p_maxspeed  <- ggplot(i, aes(x=group,y=abs(list_maxspeed_tot), fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Maximum speed (\u03bCm/s)", title =paste0("Fastest particle in\n", title))
    p_meanspeed  <- ggplot(i, aes(x=group,y=list_speed_tot, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            labs(x=NULL,  y ="Mean speed (\u03bCm/s)", title =paste0("Mean speed in\n", title))
    
    my_comparisons = list(c("CTR", "C9"))
    p_stationary2  <- ggplot(df, aes(x=col,y=prop_stationary, fill=col)) + geom_boxplot() +    theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Proportion", title =paste0("Proportion of stationary particles in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_total2  <- ggplot(i, aes(x=col,y=list_total_particles, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+ scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Total particles", title =paste0("Total particles per experiment in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_retro2  <- ggplot(i, aes(x=col,y=prop_right, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Proportion", title =paste0("Retrograde transport in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_pauseall2  <- ggplot(i, aes(x=col,y=pauses_all_particles, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Proportion of time in pause/all particles", title =paste0("Pause time (all particles) in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_pausemoving2  <- ggplot(i, aes(x=col,y=pauses_moving_particles, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Proportion of time in pause/moving particles", title =paste0("Pause time (moving particles) in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_meandist2  <- ggplot(i, aes(x=col,y=distance_total_pp, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Mean distance travelled (\u03bCm)", title =paste0("Distance travelled in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_meanmaxspeed2  <- ggplot(i, aes(x=col,y=list_avmaxspeed_tot, fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Mean maximum speed (\u03bCm/s)", title =paste0("Mean maximum speed in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_maxspeed2  <- ggplot(i, aes(x=col,y=abs(list_maxspeed_tot), fill=col)) + geom_boxplot() + theme_cowplot(rel_large = 1)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+scale_fill_manual(values=colour)+
            labs(x=NULL,  y ="Maximum speed (\u03bCm/s)", title =paste0("Fastest particle in\n", title))+ stat_compare_means(comparisons = my_comparisons)
    p_meanspeed2  <- ggplot(i, aes(x=col,y=list_speed_tot, fill=col)) + geom_boxplot() +scale_fill_manual(values=colour)+
            geom_point(position = position_jitter()) + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))+ theme_cowplot(rel_large = 1)+
            labs(x=NULL,  y ="Mean speed (\u03bCm/s)", title =paste0("Mean speed in\n", title))+ stat_compare_means(comparisons = my_comparisons) + 
      stat_summary(fun.y=mean, geom="text", show_guide = FALSE, hjust=-.5, vjust=-1, aes( label=round(..y.., digits=2)))
    
    graphs[[j]] <- list(p_stationary , p_total , p_retro , p_pauseall , p_pausemoving , p_meandist , p_meanmaxspeed ,p_maxspeed , p_meanspeed,
                        p_stationary2 , p_total2 , p_retro2 , p_pauseall2 , p_pausemoving2 , p_meandist2 , p_meanmaxspeed2 ,p_maxspeed2 , p_meanspeed2)
    names(graphs[[j]]) <- list("p_stationary" , "p_total" , "p_retro" , "p_pauseall" , "p_pausemoving" , "p_meandist" , "p_meanmaxspeed" ,"p_maxspeed" , "p_meanspeed",
                               "p_stationary2" , "p_total2" , "p_retro2" , "p_pauseall2" , "p_pausemoving2" , "p_meandist2" , "p_meanmaxspeed2" ,"p_maxspeed2" , "p_meanspeed2")
    
    svg(paste0(string,"_",j,".svg"), width = 40, height = 15)
    print(p_stationary + p_total + p_retro + p_pauseall + p_pausemoving + p_meandist + p_meanmaxspeed +p_maxspeed + p_meanspeed+
            p_stationary2 + p_total2 + p_retro2 + p_pauseall2 + p_pausemoving2 + p_meandist2 + p_meanmaxspeed2 +p_maxspeed2 + p_meanspeed2 +  plot_layout(ncol = 9)+
            plot_annotation(title = paste0(j),theme = theme(plot.title = element_text(hjust = 0.5,size = 28))))
    dev.off()
    
  }
  return(graphs)
}

exploratory_raw <- function(df, title, string){
  df$group <- "PT-CRISPR (C22)"
  df[grep("C2|C92",df$line),]$group <- "PT-C9-2"
  df[grep("C1|C91",df$line),]$group <- "PT-C9-1"
  if(sum(grepl("C4|C94",df$line))!=0){
    df[grep("C4|C94",df$line),]$group <- "PT-C9-4"
  }
  
  if(sum(grepl("OX",df$line))!=0){
    df[grep("OX",df$line),]$group <- "CTR-OX3"
  }
  if(sum(grepl("180",df$line))!=0){
    df[grep("180",df$line),]$group <- "CTR-180"
  }
  df[grep("840",df$line),]$group <- "CTR-840"
  if(sum(grepl("841",df$line))!=0){
  df[grep("841",df$line),]$group <- "CTR-841"
  }
  if(sum(grepl("856",df$line))!=0){
    df[grep("856|853",df$line),]$group <- "CTR-856"
  }
  
  df$col <- "ISO"
  df[grep("C",df$line),]$col <- "C9"
  df[grep("OX|85|84|18",df$line),]$col <- "CTR"
  
  df$abs_averagespeed <- abs(df$average_speed)
  df <- df[(df$col != "ISO"),]  
  
  df$col <- factor(df$col, levels=c("CTR", "C9"))
  
  df$distance_traveled = df$distance_traveled * 0.093
  df$max_speed = df$max_speed * 0.0465
  df$average_speed = df$average_speed * 0.0465
  
  
  
  anterograde <- df[df$direction == "Left",]
  retrograde <- df[df$direction == "Right",]
  total <- df
  my_list <- list(anterograde, retrograde,total)
  names(my_list) <- c("anterograde", "retrograde","total")
  

  p_retrograde <- list()
  p_anterograde <- list()
  p_total <- list()
  graphs <- list(p_anterograde, p_retrograde,p_total)  
  names(graphs)  <- c("anterograde", "retrograde","total")
  
  if(grepl("sn",string, ignore.case = TRUE)){
    colour = c("#8E98C0", "dodgerblue")
  }
  else{
    colour = c("#C09F98", "indianred")
  }
  
  for(j in names(my_list)){
    i = my_list[[j]]
    p_pause <- ggplot(i, aes(x=group,y=pause_time/total_frames, fill=col))  + geom_violin() + geom_boxplot(width = 0.1,) +theme_cowplot(rel_large = 1)+
              theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_fill_manual(values=colour)+
              labs(x=NULL,  y ="Proportion of time in pause/moving particles", title =paste0("Pause time (moving particles) in\n", title)) +
      stat_summary(aes(label =round(..y.., 2), y = stage(pause_time/total_frames, after_stat = 0)),size=3, fun = mean, geom = "label") 
    p_distance <- ggplot(i, aes(x=group,y=abs(distance_traveled), fill=col)) + geom_violin() + geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+
              theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=colour)+
              labs(x=NULL,  y ="Mean distance travelled (\u03bCm)", title =paste0("Distance travelled in\n", title))+
      stat_summary(aes(label =round(..y.., 1), y = stage(abs(distance_traveled), after_stat = 0)),size=3, fun = mean, geom = "label") 
    p_fastest <- ggplot(i, aes(x=group,y=abs(max_speed), fill=col))  + geom_violin() + geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+
      theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_fill_manual(values=colour)+
      labs(x=NULL,  y ="Maximum speed (\u03bCm/s)", title =paste0("Maximum speed in\n", title))+
      stat_summary(aes(label =round(..y.., 2), y = stage(abs(max_speed), after_stat = 0)),size=3, fun = mean, geom = "label") 
    p_average <- ggplot(i, aes(x=group,y=abs(average_speed), fill=col))  + geom_violin() + geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+
      theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_fill_manual(values=colour)+
      labs(x=NULL,  y ="Mean speed (\u03bCm/s)", title =paste0("Mean speed in\n", title))+ 
      stat_summary(aes(label =round(..y.., 2), y = stage(abs(average_speed), after_stat = 0)),size=3, fun = mean, geom = "label") 
    #my_comparisons = list(c("CTR", "C9"),c("ISO", "C9"), c("CTR","ISO"))
    my_comparisons = list(c("CTR", "C9"))
    p_pause2 <- ggplot(i, aes(x=col,y=pause_time/total_frames, fill=col))  + geom_violin() + 
      geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+ theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14))+
      labs(x=NULL,  y ="Proportion of time in pause/moving particles", title =paste0("Pause time (moving particles) in\n", title)) + 
      stat_compare_means(comparisons = my_comparisons)+  scale_fill_manual(values=colour)+
      stat_summary(aes(label =round(..y.., 2), y = stage(pause_time/total_frames, after_stat = 0)), fun = mean, geom = "label") 
    p_distance2 <- ggplot(i, aes(x=col,y=abs(distance_traveled), fill=col)) + geom_violin() +
      geom_boxplot(width = 0.1) + theme_cowplot(rel_large = 1)+ scale_fill_manual(values=colour)+
      theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14))+ labs(x=NULL,  y ="Mean distance travelled (\u03bCm)", title =paste0("Distance travelled in\n", title)) + 
      stat_compare_means(comparisons = my_comparisons)+ stat_compare_means(comparisons = my_comparisons) + 
      stat_summary(aes(label =round(..y.., 1), y = stage(abs(distance_traveled), after_stat = 0)), fun = mean, geom = "label") 
    p_fastest2 <- ggplot(i, aes(x=col,y=abs(max_speed), fill=col))  + geom_violin() + geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+ scale_fill_manual(values=colour)+
            theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14))+ labs(x=NULL,  y ="Maximum speed (\u03bCm/s)", title =paste0("Maximum speed in\n", title)) + stat_compare_means(comparisons = my_comparisons)+
      stat_summary(aes(label =round(..y.., 2), y = stage(abs(max_speed), after_stat = 0)), fun = mean, geom = "label") 
    p_average2 <- ggplot(i, aes(x=col,y=abs(average_speed), fill=col))  + geom_violin() + geom_boxplot(width = 0.1) +theme_cowplot(rel_large = 1)+ scale_fill_manual(values=colour)+
            theme(plot.title = element_text(hjust = 0.5), legend.position="none", text = element_text(size = 14))+ labs(x=NULL,  y ="Mean speed (\u03bCm/s)", title =paste0("Mean speed in\n", title)) + stat_compare_means(comparisons = my_comparisons)+
            stat_summary(aes(label =round(..y.., 2), y = stage(abs(average_speed), after_stat = 0)), fun = mean, geom = "label") 
    graphs[[j]] <- list(p_pause, p_distance,p_fastest,p_average,p_pause2, p_distance2,p_fastest2,p_average2)
    names(graphs[[j]]) <- c("p_pause", "p_distance","p_fastest","p_average","p_pause2", "p_distance2","p_fastest2","p_average2")

    svg(paste0(string,"_",j,".svg"), width = 20, height = 10)
      print(p_pause + p_distance +p_fastest + p_average +p_pause2 + p_distance2 +p_fastest2 + p_average2 + plot_layout(ncol = 4)+
        plot_annotation(title = paste0(j),theme = theme(plot.title = element_text(hjust = 0.5,size = 24))))
    dev.off()
    
    
  }
  ctr <- abs(df[df$col == "CTR",]$average_speed)
  c9 <- abs(df[df$col == "C9",]$average_speed)
  iso <- abs(df[df$col == "ISO",]$average_speed)
  sq <- seq(max(length(ctr), length(c9), length(iso)))
  write_tsv(data.frame(ctr[sq], c9[sq], iso[sq]),paste0(title,".tsv"))
  print(t.test(df[df$col == "C9",]$abs_averagespeed, df[df$col == "CTR",]$abs_averagespeed))
  #df2 <-df2[grep("Feb",  df2$differentiation, invert = TRUE), ]
  return(graphs)
}




hct_mn_raw <- df3[grepl("Motor",df3$differentiation),]
hct_mn_raw <- hct_mn_raw[grepl("HcT",hct_mn_raw$target),]
hct_mn_raw <- hct_mn_raw[grepl("Moving",hct_mn_raw$type),]
p_hct_mn_raw <- exploratory_raw(hct_mn_raw, "Motor Neurons (Cholera Toxin)", "raw_hct_mn")

hct_mn <- summary_lines[grepl("Motor",summary_lines$differentiation),]
hct_mn <- hct_mn[grepl("HcT",hct_mn$target),]
#two instances of C22 for one of the lines
hct_mn[grep("C22",hct_mn$line),]$line <- "C22"
hct_mn <-hct_mn[grep("Feb",  hct_mn$differentiation, invert = TRUE), ]
p_hct_mn <- exploratory(hct_mn, 1, "Motor Neurons (Cholera Toxin)", "group_hct_mn")

hct_sn_raw <- df3[grepl("Sensory",df3$differentiation),]
hct_sn_raw <- hct_sn_raw[grepl("HcT",hct_sn_raw$target),]
hct_sn_raw <- hct_sn_raw[grepl("Moving",hct_sn_raw$type),]
p_hct_sn_raw <- exploratory_raw(hct_sn_raw, "Sensory Neurons (Cholera Toxin)", "raw_hct_sn")

hct_SN <- summary_lines[grepl("Sensory",summary_lines$differentiation),]
hct_SN <- hct_SN[grepl("HcT",hct_SN$target),]
p_hct_SN <- exploratory(hct_SN, 1, "Sensory Neurons (Cholera Toxin)", "group_hct_sn")


lyso_mn_raw <- df3[grepl("Motor",df3$differentiation),]
lyso_mn_raw <- lyso_mn_raw[grepl("Lyso",lyso_mn_raw$target),]
lyso_mn_raw <- lyso_mn_raw[grepl("Moving",lyso_mn_raw$type),]
p_lyso_mn_raw <- exploratory_raw(lyso_mn_raw, "Motor Neurons (Lysotracker)", "raw_lyso_mn")

lyso_MN <- summary_lines[grepl("Motor",summary_lines$differentiation),]
lyso_MN <- lyso_MN[grepl("Lyso",lyso_MN$target),]
p_lyso_MN <- exploratory(lyso_MN, 1, "Motor Neurons (Lysotracker)", "group_lyso_mn")


lyso_sn_raw <- df3[grepl("Sensory",df3$differentiation),]
lyso_sn_raw <- lyso_sn_raw[grepl("Lyso",lyso_sn_raw$target),]
lyso_sn_raw <- lyso_sn_raw[grepl("Moving",lyso_sn_raw$type),]
p_lyso_sn_raw <- exploratory_raw(lyso_sn_raw, "Sensory Neurons (Lysotracker)", "raw_lyso_sn")

lyso_SN <- summary_lines[grepl("Sensory",summary_lines$differentiation),]
lyso_SN <- lyso_SN[grepl("Lyso",lyso_SN$target),]
lyso_SN <- exploratory(lyso_SN, 1, "Sensory Neurons (Lysotracker)","group_lyso_sn")



mito_mn_raw <- df3[grepl("Motor",df3$differentiation),]
mito_mn_raw <- mito_mn_raw[grepl("Mito",mito_mn_raw$target),]
mito_mn_raw <- mito_mn_raw[grepl("Moving",mito_mn_raw$type),]
p_mito_mn_raw <- exploratory_raw(mito_mn_raw, "Motor Neurons (Mitotracker)", "raw_mito_mn")

mito_MN <- summary_lines[grepl("Motor",summary_lines$differentiation),]
mito_MN <- mito_MN[grepl("Mito",mito_MN$target),]
p_mito_MN <- exploratory(mito_MN, 1, "Motor Neurons (Mitotracker)","group_mito_mn")


mito_sn_raw <- df3[grepl("Sensory",df3$differentiation),]
mito_sn_raw <- mito_sn_raw[grepl("Mito",mito_sn_raw$target),]
mito_sn_raw <- mito_sn_raw[grepl("Moving",mito_sn_raw$type),]
p_mito_sn_raw <- exploratory_raw(mito_sn_raw, "Sensory Neurons (Mitotracker)", "raw_mito_sn")

mito_SN <- summary_lines[grepl("Sensory",summary_lines$differentiation),]
mito_SN <- mito_SN[grepl("Mito",mito_SN$target),]
p_mito_SN <- exploratory(mito_SN, 1, "Sensory Neurons (Mitotracker)", "group_mito_sn")


syto_mn_raw <- df3[grepl("Motor",df3$differentiation),]
syto_mn_raw <- syto_mn_raw[grepl("Syto",syto_mn_raw$target),]
syto_mn_raw <- syto_mn_raw[grepl("Moving",syto_mn_raw$type),]
p_syto_mn_raw <- exploratory_raw(syto_mn_raw, "Motor Neurons (Syto RNA Select)" , "raw_syto_mn")

syto_MN <- summary_lines[grepl("Motor",summary_lines$differentiation),]
syto_MN <- syto_MN[grepl("Syto",syto_MN$target),]
p_syto_MN <- exploratory(syto_MN, 1, "Motor Neurons (Syto RNA)", "group_syto_mn")

syto_sn_raw <- df3[grepl("Sensory",df3$differentiation),]
syto_sn_raw <- syto_sn_raw[grepl("Syto",syto_sn_raw$target),]
syto_sn_raw <- syto_sn_raw[grepl("Moving",syto_sn_raw$type),]
p_syto_sn_raw <- exploratory_raw(syto_sn_raw, "Sensory Neurons (Syto RNA Select)", "raw_syto_sn")

syto_SN <- summary_lines[grepl("Sensory",summary_lines$differentiation),]
syto_SN <- syto_SN[grepl("Syto",syto_SN$target),]
p_syto_SN <- exploratory(syto_SN, 1, "Sensory Neurons (Syto RNA)", "group_syto_sn")

mn_vs_sens <- df2[grepl("Syto",df2$target),]
mn_vs_sens <- mn_vs_sens[grepl("Moving",mn_vs_sens$type),]
mn_vs_sens$group <- "Motor"
mn_vs_sens[grep("Sens",mn_vs_sens$differentiation),]$group <- "Sensory"
print(ggplot(mn_vs_sens, aes(x=group,y=pause_time/total_frames, fill =differentiation)) + geom_boxplot() + theme(text = element_text(size = 14))+
        labs(x="Celltype",  y ="Proportion of time in pause/moving particles", title ="Pause time (moving particles)"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(distance_traveled), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x=NULL,  y ="Mean distance travelled", title ="Distance travelled"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(max_speed), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x=NULL,  y ="Maximum speed", title ="Fastest particle"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(average_speed), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x=NULL,  y ="Mean speed (\u03bCm/s)", title ="Mean speed"))

save(file = paste0("maindataframe.RData"), df, df2)


(p_lyso_mn_raw$total$p_average2+labs(title="total")) + (((p_lyso_mn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_lyso_mn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_lyso_mn_raw$total$p_average+labs(title = "breakdown by subject")))+ plot_spacer()+
  (p_lyso_sn_raw$total$p_average2+labs(title="total")) + (((p_lyso_sn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_lyso_sn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_lyso_sn_raw$total$p_average+labs(title = "breakdown by subject"))) + plot_spacer()+(p_hct_sn_raw$total$p_average2/p_hct_sn_raw$total$p_average2)+
  plot_layout(widths = c(0.5,1,0.3,0.5,.3,.5))
  

group1 <- p_lyso_mn_raw$total$p_average2+labs(title="total") + (((p_lyso_mn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_lyso_mn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_lyso_mn_raw$total$p_average+labs(title = "breakdown by subject"))) + plot_annotation(title = "Lysotracker (motor neurons)", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))) + plot_layout(widths = c(0.5,1))
group2 <-   (p_lyso_sn_raw$total$p_average2+labs(title="total")) + (((p_lyso_sn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_lyso_sn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_lyso_sn_raw$total$p_average+labs(title = "breakdown by subject"))) + plot_annotation(title = "Lysotracker (sensory neurons)", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))) +plot_layout(widths = c(0.5,1))
CairoSVG("Figure6A.svg", width = 16, height = 9,pointsize = 16)
plot_grid(group1, group2,ncol = 2, rel_widths = c(2,2))
dev.off()

group1 <- p_mito_mn_raw$total$p_average2+labs(title="total") + (((p_mito_mn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_mito_mn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_mito_mn_raw$total$p_average+labs(title = "breakdown by subject"))) + plot_annotation(title = "Mitotracker (motor neurons)", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))) + plot_layout(widths = c(0.5,1))
group2 <-   (p_mito_sn_raw$total$p_average2+labs(title="total")) + (((p_mito_sn_raw$anterograde$p_average2+labs(title="anterograde"))+(p_mito_sn_raw$retrograde$p_average2+labs(title="retrograde")))/(p_mito_sn_raw$total$p_average+labs(title = "breakdown by subject"))) + plot_annotation(title = "Mitotracker (sensory neurons)", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))) +plot_layout(widths = c(0.5,1))
CairoSVG("Figure6B.svg", width = 16, height = 9,pointsize = 16)
plot_grid(group1, group2, ncol = 2)
dev.off()


group3a <- (p_hct_mn_raw$total$p_average2+labs(title="motor neurons"))/(p_hct_sn_raw$total$p_average2+labs(title="sensory neurons")) + plot_annotation(title = "Tetanus Toxin", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
group3b <- ((p_lyso_mn_raw$total$p_pause2+labs(title="Lysotracker MNs"))+(p_mito_mn_raw$total$p_pause2+labs(title="Mitotracker MNs")))/((p_lyso_sn_raw$total$p_pause2+labs(title="Lysotracker SNs"))+(p_mito_sn_raw$total$p_pause2+labs(title="Mitotracker SNs"))) + plot_annotation(title = "Pause time", theme = theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold")))
CairoSVG("FigureS6.svg", width = 16, height = 12,pointsize = 16)
plot_grid(group3a, group3b, ncol = 2,rel_widths = c(2,4))
dev.off()
