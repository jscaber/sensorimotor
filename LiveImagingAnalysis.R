library("tidyverse")


file_list <- list.files(path = "/Volumes/Elements/Microscopy/Live/", pattern = "_tracked.txt$", recursive = TRUE)

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
  filename <- paste0("/Volumes/Elements/Microscopy/Live/",fn)
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

list_diff <- list()
list_target <- list()
list_line <- list()
list_area <-list()
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
        list_diff[[k]] <- uniq_diff
        list_target[[k]] <- uniq_target
        list_line[[k]] <- uniq_line
        list_area[[k]] <- uniq_area
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

summary_areas <- cbind(list_diff,list_target,list_line,list_area,list_plane,list_stationary,list_right,list_distance_tot,list_distance_right,list_distance_left,list_distance_stat,list_speed_tot,list_speed_right,list_speed_left,list_maxspeed_tot,list_maxspeed_right,list_maxspeed_left,list_avmaxspeed_tot,list_avmaxspeed_right,list_avmaxspeed_left,list_pause_time_moving,list_pause_time_total,list_total_frames_moving,list_total_frames,list_total_particles)
summary_areas <- as.data.frame(summary_areas)
summary_areas[,6:25] = apply(summary_areas[,6:25], 2, function(x) as.numeric(x))
summary_areas[,1:5] = apply(summary_areas[,1:5], 2, function(x) as.character(x))
summary_areas <- summarise_axonal_transport(summary_areas)




list_diff <- list()
list_target <- list()
list_line <- list()
list_area <-list()
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
      ds <- tmp_line[tmp_line$type == "Dynamic Stationary",]
      mov <- tmp_line[tmp_line$type == "Moving",]
      mov_r <- mov[mov$direction == "Right",]
      mov_l <- mov[mov$direction == "Left",]
      list_diff[[k]] <- uniq_diff
      list_target[[k]] <- uniq_target
      list_line[[k]] <- uniq_line
      list_area[[k]] <- uniq_area
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

summary_lines <- cbind(list_diff,list_target,list_line,list_area,list_plane,list_stationary,list_right,list_distance_tot,list_distance_right,list_distance_left,list_distance_stat,list_speed_tot,list_speed_right,list_speed_left,list_maxspeed_tot,list_maxspeed_right,list_maxspeed_left,list_avmaxspeed_tot,list_avmaxspeed_right,list_avmaxspeed_left,list_pause_time_moving,list_pause_time_total,list_total_frames_moving,list_total_frames,list_total_particles)
summary_lines <- as.data.frame(summary_lines)
summary_lines[,6:25] = apply(summary_lines[,6:25], 2, function(x) as.numeric(x))
summary_lines[,1:5] = apply(summary_lines[,1:5], 2, function(x) as.character(x))



summary_lines <- summarise_axonal_transport(summary_lines)


exploratory <- function(df, stationary_threshold = 1, title){
  df$group <- "PT-CRISPR (C22)"
  df[grep("C2|C92",df$list_line),]$group <- "PT-C9-2"
  df[grep("C1|C91",df$list_line),]$group <- "PT-C9-1"
  df[grep("C4|C94",df$list_line),]$group <- "PT-C9-4"
  df[grep("OX",df$list_line),]$group <- "CTR-OX3"
  if(sum(grepl("180",df$list_line))!=0){
    df[grep("180",df$list_line),]$group <- "CTR-180"
  }
  df[grep("840",df$list_line),]$group <- "CTR-840"
  df[grep("841",df$list_line),]$group <- "CTR-841"
  df[grep("856|853",df$list_line),]$group <- "CTR-856"
  
  df$col <- "CRISPR"
  df[grep("C",df$list_line),]$col <- "C9"
  df[grep("OX|85|84|18",df$list_line),]$col <- "CTR"
  
  # Remove dead cells
  df2 <- df[df$prop_stationary < stationary_threshold,]
  
  print(ggplot(df, aes(x=group,y=prop_stationary, fill=col)) + geom_boxplot() +   
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Proportion", title =paste("Proportion of stationary particles in", title)))
  print(ggplot(df2, aes(x=group,y=list_total_particles, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Total particles", title =paste("Total particles per experiment in", title)))
  print(ggplot(df2, aes(x=group,y=prop_right, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Proportion", title =paste("Retrograde transport in", title)))
  print(ggplot(df2, aes(x=group,y=pauses_all_particles, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Proportion of time in pause/all particles", title =paste("Pause time (all particles) in", title)))
  print(ggplot(df2, aes(x=group,y=pauses_moving_particles, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Proportion of time in pause/moving particles", title =paste("Pause time (moving particles) in", title)))
  print(ggplot(df2, aes(x=group,y=distance_total_pp, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Mean distance travelled", title =paste("Distance travelled in", title)))
  print(ggplot(df2, aes(x=group,y=list_avmaxspeed_tot, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Mean maximum speed", title =paste("Mean maximum speed in", title)))
  print(ggplot(df2, aes(x=group,y=abs(list_maxspeed_tot), fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Maximum speed", title =paste("Fastest particle in", title)))
  print(ggplot(df2, aes(x=group,y=list_speed_tot, fill=col)) + geom_boxplot() +
          geom_point(position = position_jitter()) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Mean speed", title =paste("Mean speed in", title)))
}

exploratory_raw <- function(df, title){
  df$group <- "PT-CRISPR (C22)"
  df[grep("C2|C92",df$line),]$group <- "PT-C9-2"
  df[grep("C1|C91",df$line),]$group <- "PT-C9-1"
  df[grep("C4|C94",df$line),]$group <- "PT-C9-4"
  df[grep("OX",df$line),]$group <- "CTR-OX3"
  if(sum(grepl("180",df$line))!=0){
    df[grep("180",df$line),]$group <- "CTR-180"
  }
  df[grep("840",df$line),]$group <- "CTR-840"
  df[grep("841",df$line),]$group <- "CTR-841"
  df[grep("856|853",df$line),]$group <- "CTR-856"
  
  df$col <- "CRISPR"
  df[grep("C",df$line),]$col <- "C9"
  df[grep("OX|85|84|18",df$line),]$col <- "CTR"

  print(ggplot(df, aes(x=group,y=pause_time/total_frames, fill=col)) + geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(), size = .5) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Proportion of time in pause/moving particles", title =paste("Pause time (moving particles) in", title)))
  print(ggplot(df, aes(x=group,y=abs(distance_traveled), fill=col)) + geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(), size = .5) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Mean distance travelled", title =paste("Distance travelled in", title)))
  print(ggplot(df, aes(x=group,y=abs(max_speed), fill=col)) + geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(), size = .5) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Maximum speed", title =paste("Fastest particle in", title)))
  print(ggplot(df, aes(x=group,y=abs(average_speed), fill=col)) + geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(), size = .5) + theme(text = element_text(size = 14))+
          labs(x="Patient",  y ="Mean speed", title =paste("Mean speed in", title)))
}



hct_mn_raw <- df2[grepl("Motor",df2$differentiation),]
hct_mn_raw <- hct_mn_raw[grepl("HcT",hct_mn_raw$target),]
hct_mn_raw <- hct_mn_raw[grepl("Moving",hct_mn_raw$type),]
exploratory_raw(hct_mn_raw, "Motor Neurons (Cholera Toxin)")

hct_mn <- summary_lines[grepl("Motor",summary_lines$list_diff),]
hct_mn <- hct_mn[grepl("HcT",hct_mn$list_target),]
#two instances of C22 for one of the lines
hct_mn[grep("C22",hct_mn$list_line),]$list_line <- "C22"
exploratory(hct_mn, 1, "Motor Neurons (Cholera Toxin)")

hct_sn_raw <- df2[grepl("Motor",df2$differentiation),]
hct_sn_raw <- hct_sn_raw[grepl("HcT",hct_sn_raw$target),]
hct_sn_raw <- hct_sn_raw[grepl("Moving",hct_sn_raw$type),]
exploratory_raw(hct_sn_raw, "Sensory Neurons (Cholera Toxin)")

hct_SN <- summary_lines[grepl("Sensory",summary_lines$list_diff),]
hct_SN <- hct_SN[grepl("HcT",hct_SN$list_target),]
exploratory(hct_SN, 1, "Sensory Neurons (Cholera Toxin)")

lyso_mn_raw <- df2[grepl("Motor",df2$differentiation),]
lyso_mn_raw <- lyso_mn_raw[grepl("Lyso",lyso_mn_raw$target),]
lyso_mn_raw <- lyso_mn_raw[grepl("Moving",lyso_mn_raw$type),]
exploratory_raw(lyso_mn_raw, "Motor Neurons (Lysotracker)")

lyso_MN <- summary_lines[grepl("Motor",summary_lines$list_diff),]
lyso_MN <- lyso_MN[grepl("Lyso",lyso_MN$list_target),]
exploratory(lyso_MN, 1, "Motor Neurons (Lysotracker)")

lyso_sn_raw <- df2[grepl("Sensory",df2$differentiation),]
lyso_sn_raw <- lyso_sn_raw[grepl("Lyso",lyso_sn_raw$target),]
lyso_sn_raw <- lyso_sn_raw[grepl("Moving",lyso_sn_raw$type),]
exploratory_raw(lyso_sn_raw, "Sensory Neurons (Lysotracker)")

lyso_SN <- summary_lines[grepl("Sensory",summary_lines$list_diff),]
lyso_SN <- lyso_SN[grepl("Lyso",lyso_SN$list_target),]
exploratory(lyso_SN, 1, "Sensory Neurons (Lysotracker)")

mito_mn_raw <- df2[grepl("Motor",df2$differentiation),]
mito_mn_raw <- mito_mn_raw[grepl("Mito",mito_mn_raw$target),]
mito_mn_raw <- mito_mn_raw[grepl("Moving",mito_mn_raw$type),]
exploratory_raw(mito_mn_raw, "Motor Neurons (Mitotracker)")

mito_MN <- summary_lines[grepl("Motor",summary_lines$list_diff),]
mito_MN <- mito_MN[grepl("Mito",mito_MN$list_target),]
exploratory(mito_MN, 1, "Motor Neurons (Mitotracker)")

mito_sn_raw <- df2[grepl("Sensory",df2$differentiation),]
mito_sn_raw <- mito_sn_raw[grepl("Mito",mito_sn_raw$target),]
mito_sn_raw <- mito_sn_raw[grepl("Moving",mito_sn_raw$type),]
exploratory_raw(mito_sn_raw, "Sensory Neurons (Mitotracker)")

mito_SN <- summary_lines[grepl("Sensory",summary_lines$list_diff),]
mito_SN <- mito_SN[grepl("Mito",mito_SN$list_target),]
exploratory(mito_SN, 1, "Sensory Neurons (Mitotracker)")

syto_mn_raw <- df2[grepl("Motor",df2$differentiation),]
syto_mn_raw <- syto_mn_raw[grepl("Syto",syto_mn_raw$target),]
syto_mn_raw <- syto_mn_raw[grepl("Moving",syto_mn_raw$type),]
exploratory_raw(syto_mn_raw, "Motor Neurons (Syto RNA Select)")

syto_MN <- summary_lines[grepl("Motor",summary_lines$list_diff),]
syto_MN <- syto_MN[grepl("Syto",syto_MN$list_target),]
exploratory(syto_MN, 1, "Motor Neurons (Syto RNA)")

syto_sn_raw <- df2[grepl("Sensory",df2$differentiation),]
syto_sn_raw <- syto_sn_raw[grepl("Syto",syto_sn_raw$target),]
syto_sn_raw <- syto_sn_raw[grepl("Moving",syto_sn_raw$type),]
exploratory_raw(syto_sn_raw, "Sensory Neurons (Syto RNA Select)")

syto_SN <- summary_lines[grepl("Sensory",summary_lines$list_diff),]
syto_SN <- syto_SN[grepl("Syto",syto_SN$list_target),]
exploratory(syto_SN, 1, "Sensory Neurons (Syto RNA)")

mn_vs_sens <- df2[grepl("Syto",df2$target),]
mn_vs_sens <- mn_vs_sens[grepl("Moving",mn_vs_sens$type),]
mn_vs_sens$group <- "Motor"
mn_vs_sens[grep("Sens",mn_vs_sens$differentiation),]$group <- "Sensory"
print(ggplot(mn_vs_sens, aes(x=group,y=pause_time/total_frames, fill =differentiation)) + geom_boxplot() + theme(text = element_text(size = 14))+
        labs(x="Celltype",  y ="Proportion of time in pause/moving particles", title ="Pause time (moving particles)"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(distance_traveled), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x="Patient",  y ="Mean distance travelled", title ="Distance travelled"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(max_speed), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x="Patient",  y ="Maximum speed", title ="Fastest particle"))
print(ggplot(mn_vs_sens, aes(x=group,y=abs(average_speed), fill=differentiation)) + geom_boxplot(outlier.shape = NA) +
        theme(text = element_text(size = 14))+
        labs(x="Patient",  y ="Mean speed", title ="Mean speed"))
