library(tidyverse)
library(MuMIn)
library(ggplot2)
library(dplyr)

dat.f <- list.files(path="./Project_8_data")
dat.l <- list()
for(i in dat.f){
  metadata <- unlist(strsplit(i,'_'))
  subj <- metadata[2]
  ang <- metadata[3]
  exp <- gsub('.csv','',metadata[4])
  trial.force <- read_csv(paste0('Project_8_data/',i),col_names=FALSE)
  trial.force$X1 <- gsub("[a-zA-Z :]","",trial.force$X1)
  mvcf <- trial.force %>% 
    pull(X1) %>% 
    as.numeric() %>% 
    mean()
  dat.l[[i]] <- tibble(
    subject=subj,
    ang=as.numeric(ang),
    exp=exp,
    force=mvcf
  )
}
dat <- do.call(rbind,dat.l)


dat.p <- dat %>% 
  group_by(exp,subject) %>% 
  mutate(force=force/max(force)) %>% 
  print() 

dat.p2 <- dat.p %>% 
  na.omit()

AICs <- dat.p2%>%
  group_by(subject,exp)%>%
  summarize(
    m2=AICc(lm(force~poly(ang,2))), #second order
    m3=AICc(lm(force~poly(ang,3))), #third order
    m4=AICc(lm(force~poly(ang,4))) #fourth order
  )%>%
  pivot_longer(m2:m4,names_to="model",values_to="AICc")

x.pred <- seq(45,157.5,length.out = 1000)
fits <- dat.p2%>%
  group_by(subject,exp)%>%
  summarize(
    m2=predict(lm(force~poly(ang,2)),newdata=data.frame(ang=x.pred)), #second order
    m3=predict(lm(force~poly(ang,3)),newdata=data.frame(ang=x.pred)), #third order
    m4=predict(lm(force~poly(ang,4)),newdata=data.frame(ang=x.pred)) #fourth order
  )%>%
  pivot_longer(m2:m4,names_to="model")%>%
  group_by(subject,exp,model)%>%
  summarize(theta_max=x.pred[which.max(value)])

dat.p2 %>%
  ggplot(aes(ang,force,col=exp))+geom_point()


best.models <- fits%>%
  left_join(AICs)%>%
  group_by(subject,exp)%>%
  mutate(best=AICc==min(AICc))%>%
  filter(best==TRUE)%>%
  dplyr::select(-best)%>%
  print()

best.models%>%
  pivot_wider(id_cols=subject,names_from = exp,values_from=theta_max)%>%
  mutate(shift=fatigue-control) %>% 
  ungroup() %>%
  na.omit() %>% 
  summarise(mean.shift=mean(shift),se.shift=sd(shift)/sqrt(length(shift)))


