library(tidyverse)
library(MuMIn)
library(ggplot2)
library(dplyr)


ang <- seq(45,157.5,length.out = 11)
k <- list.files("./Project_8_data", full.names = T,pattern = ".csv")
print(k)
k.l <- list()
for(i in k){
  met.dat <- unlist(strsplit(i,"_"))
  sub <- met.dat[2]
  ang <- as.numeric(met.dat[3])
  exp <- gsub("\\..+","",met.dat[4])
  k.l[[i]] <- read_delim(i,delim = " ", col_names = c("Reading","Force","Unit"), id="Experiment",progress = FALSE) %>%select(Force)%>%
    mutate(sub=sub,ang=ang,exp=exp)
}

data <- do.call(rbind, k.l)
data <- data%>%
  group_by(sub,exp,ang)%>%
  summarise(max.force=max(abs(Force), na.rm = TRUE), n=n())
data <- na.omit(data)



data.joined <- data%>%
  group_by(sub,exp)%>%
  summarize(max.f2 = max(max.force))%>%
  left_join(data)%>%
  mutate(norm.force=max.force/max.f2)%>%
  ggplot(aes(ang, norm.force))+geom_point()+geom_point(aes(x=ang[which.max(norm.force)], y=norm.force[which.max(norm.force)]), col="red", size=4)+facet_wrap(~exp, ncol=5)


angle <- seq(45, 157.5, length.out=11)
data.joined %>%
  ggplot(aes(ang, norm.force))+geom_point()+geom_point(aes(x=ang[which.max(norm.force)], y=norm.force[which.max(norm.force)]), col="red", size=4)+facet_wrap(~exp, ncol=5)
data.joined$ang[which.max(data.joined$norm.force)]


