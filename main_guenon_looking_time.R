# Statistical analysis of guenon (Cercopithecus nictitans and C. mona) looking time experiments

# Sandra Winters <sandra.winters@bristol.ac.uk>

# Citation: 
# Winters S, Allen WL, Higham JP. 2019. The structure of species discrimination signals across a primate radiation. eLife. https://doi.org/10.7554/eLife.47428

library(tidyr)
library(lme4)
library(car)

source('formatLTdata.R')

#import & format data
dat.raw <- read.csv('looking time data.csv',header=T)

dat <- formatLTdata(dat.raw)
dat <- dat[[2]]

#separate by species
dat.putty <- droplevels(dat[dat$species=='nictitans',])
dat.mona <- droplevels(dat[dat$species=='mona',])

#PUTTY MODELS
y.putty <- cbind(dat.putty$looks_5s,dat.putty$total_looks_5s-dat.putty$looks_5s)
size.putty <- round(mean(dat.putty$total_looks_5s))

#null model
putty.null <- glmer(y.putty~1+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.null)
qqp(residuals(putty.null),"binom",size=size.putty,prob=0.5)

#LRT - image type *(trend)
putty.type <- glmer(y.putty~im_type+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.type)
qqp(residuals(putty.type),"binom",size=size.putty,prob=0.5)
anova(putty.type,putty.null)
summary(putty.type) #trend for conspecific bias

#LRT - image trait *
putty.trait <- glmer(y.putty~im_trait+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.trait)
qqp(residuals(putty.trait),"binom",size=size.putty,prob=0.5)
anova(putty.trait,putty.null)
summary(putty.trait) #significant bias for shared trait

#LRT- sex 
putty.sex <- glmer(y.putty~subject_sex+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.sex)
qqp(residuals(putty.sex),"binom",size=size.putty,prob=0.5)
anova(putty.sex,putty.null)

#LRT - age 
putty.age <- glmer(y.putty~log(age.years)+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.age)
qqp(residuals(putty.age),"binom",size=size.putty,prob=0.5)
anova(putty.age,putty.null)

#LRT - origin 
putty.origin <- glmer(y.putty~origin+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.origin)
qqp(residuals(putty.origin),"binom",size=size.putty,prob=0.5)
anova(putty.origin,putty.null)

#LRT - presentation spot (R/1 v. L/2) *
putty.pres <- glmer(y.putty~pres_spot+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.pres)
qqp(residuals(putty.pres),"binom",size=size.putty,prob=0.5)
anova(putty.pres,putty.null)
summary(putty.pres) #significant bias for image on right

#LRT - eye contact 
putty.eyecont <- glmer(y.putty~eye_contact+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.eyecont)
qqp(residuals(putty.eyecont),"binom",size=size.putty,prob=0.5)
anova(putty.eyecont,putty.null)

#LRT - familiarity 
putty.famil <- glmer(y.putty~im_familiarity+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.famil)
qqp(residuals(putty.famil),"binom",size=size.putty,prob=0.5)
anova(putty.famil,putty.null)

#LRT - sex of stimulus 
putty.stimSex <- glmer(y.putty~im_sex+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.stimSex)
qqp(residuals(putty.stimSex),"binom",size=size.putty,prob=0.5)
anova(putty.stimSex,putty.null)

#LRT - trial order  
putty.trialOrder <- glmer(y.putty~trial_order+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.trialOrder)
qqp(residuals(putty.trialOrder),"binom",size=size.putty,prob=0.5)
anova(putty.trialOrder,putty.null)

#LRT - pattern 
putty.pattern <- glmer(y.putty~pattern+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.pattern)
qqp(residuals(putty.pattern),"binom",size=size.putty,prob=0.5)
anova(putty.pattern,putty.null)

#LRT - icc 
putty.icc <- glmer(y.putty~icc+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.icc)
qqp(residuals(putty.icc),"binom",size=size.putty,prob=0.5)
anova(putty.icc,putty.null)

#model with all significant factors (incl trend for im_type)
putty.signif <- glmer(y.putty~im_type+im_trait+pres_spot+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.signif)
qqp(residuals(putty.signif),"binom",size=size.putty,prob=0.5)
anova(putty.signif,putty.null)
summary(putty.signif)

#drop each factor from significant model & run LRTs
#drop type
putty.signif.typeRm <- glmer(y.putty~im_trait+pres_spot+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.signif.typeRm)
qqp(residuals(putty.signif.typeRm),"binom",size=size.putty,prob=0.5)
anova(putty.signif.typeRm,putty.signif)

#drop trait *
putty.signif.traitRm <- glmer(y.putty~im_type+pres_spot+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.signif.traitRm)
qqp(residuals(putty.signif.traitRm),"binom",size=size.putty,prob=0.5)
anova(putty.signif.traitRm,putty.signif)

#drop presentation spot *
putty.signif.presSpotRm <- glmer(y.putty~im_type+im_trait+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.signif.presSpotRm)
qqp(residuals(putty.signif.presSpotRm),"binom",size=size.putty,prob=0.5)
anova(putty.signif.presSpotRm,putty.signif)

#final model
putty.final <- glmer(y.putty~im_trait+pres_spot+(1|group/subject/trial_subject),data=dat.putty,family=binomial)
plot(putty.final)
qqp(residuals(putty.final),"binom",size=size.putty,prob=0.5)
anova(putty.final,putty.null)
summary(putty.final) #significant bias for shared trait, image on right

exp(coef(summary(putty.final))[ , "Estimate"]) #odds ratio
1/exp(coef(summary(putty.final))[ , "Estimate"]) 

#look at only con v. hetero w/o nose spot (trial 1)
dat.putty.t1 <- dat.putty[dat.putty$condition=='not_shared',]
y.putty.t1 <- cbind(dat.putty.t1$looks_5s,dat.putty.t1$total_looks_5s-dat.putty.t1$looks_5s)

putty.null.t1 <- glmer(y.putty.t1~pres_spot+(1|group/subject/trial_subject),data=dat.putty.t1,family=binomial)
plot(putty.null.t1)
qqp(residuals(putty.null.t1),"binom",size=size.putty,prob=0.5)
summary(putty.null.t1) #significant right side bias

putty.type.t1 <- glmer(y.putty.t1~im_type+pres_spot+(1|group/subject/trial_subject),data=dat.putty.t1,family=binomial) #note: type & trait are the same here
plot(putty.type.t1)
qqp(residuals(putty.type.t1),"binom",size=size.putty,prob=0.5)
anova(putty.type.t1,putty.null.t1)
summary(putty.type.t1) #significant conspecific bias, right side bias

exp(coef(summary(putty.type.t1))[ , "Estimate"]) #odds ratio
1/exp(coef(summary(putty.type.t1))[ , "Estimate"])


#MONA MODELS
y.mona <- cbind(dat.mona$looks_5s,dat.mona$total_looks_5s-dat.mona$looks_5s)
size.mona <- round(mean(dat.mona$total_looks_5s))

#null model
mona.null <- glmer(y.mona~1+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.null)
qqp(residuals(mona.null),"binom",size=size.mona,prob=0.5)

#LRT - image type *
mona.type <- glmer(y.mona~im_type+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.type)
qqp(residuals(mona.type),"binom",size=size.mona,prob=0.5)
anova(mona.type,mona.null)
summary(mona.type) #significant conspecific bias

#LRT - image trait *
mona.trait <- glmer(y.mona~im_trait+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.trait)
qqp(residuals(mona.trait),"binom",size=size.mona,prob=0.5)
anova(mona.trait,mona.null)
summary(mona.trait) #significant bias for different trait

#LRT - type*trait interaction
mona.interact <- glmer(y.mona~im_type*im_trait+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.interact)
qqp(residuals(mona.interact),"binom",size=size.mona,prob=0.5)
anova(mona.interact,mona.null)
anova(mona.interact,mona.type)
anova(mona.interact,mona.trait)
summary(mona.interact)

#LRT - sex 
mona.sex <- glmer(y.mona~subject_sex+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.sex)
qqp(residuals(mona.sex),"binom",size=size.mona,prob=0.5)
anova(mona.sex,mona.null)

#LRT - age 
mona.age <- glmer(y.mona~log(age.years)+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.age)
qqp(residuals(mona.age),"binom",size=size.mona,prob=0.5)
anova(mona.age,mona.null)

#LRT - origin 
mona.origin <- glmer(y.mona~origin+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.origin)
qqp(residuals(mona.origin),"binom",size=size.mona,prob=0.5)
anova(mona.origin,mona.null)

#LRT - presentation spot (R/1 v. L/2)
mona.pres <- glmer(y.mona~pres_spot+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.pres)
qqp(residuals(mona.pres),"binom",size=size.mona,prob=0.5)
anova(mona.pres,mona.null)

#LRT - eye contact
mona.eyecont <- glmer(y.mona~eye_contact+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.eyecont)
qqp(residuals(mona.eyecont),"binom",size=size.mona,prob=0.5)
anova(mona.eyecont,mona.null)

#LRT - sex of stimulus 
mona.stimSex <- glmer(y.mona~im_sex+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.stimSex)
qqp(residuals(mona.stimSex),"binom",size=size.mona,prob=0.5)
anova(mona.stimSex,mona.null)

#LRT - trial order  
mona.trialOrder <- glmer(y.mona~trial_order+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.trialOrder)
qqp(residuals(mona.trialOrder),"binom",size=size.mona,prob=0.5)
anova(mona.trialOrder,mona.null)

#LRT - pattern 
mona.pattern <- glmer(y.mona~pattern+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.pattern)
qqp(residuals(mona.pattern),"binom",size=size.mona,prob=0.5)
anova(mona.pattern,mona.null)

#LRT - icc 
mona.icc <- glmer(y.mona~icc+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.icc)
qqp(residuals(mona.icc),"binom",size=size.mona,prob=0.5)
anova(mona.icc,mona.null)

#model with all significant factors 
mona.signif <- glmer(y.mona~im_type+im_trait+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.signif)
qqp(residuals(mona.signif),"binom",size=size.mona,prob=0.5)
anova(mona.signif,mona.null)
summary(mona.signif)

#drop each factor from significant model & run LRTs
#drop type *
mona.signif.typeRm <- glmer(y.mona~im_trait+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.signif.typeRm)
qqp(residuals(mona.signif.typeRm),"binom",size=size.mona,prob=0.5)
anova(mona.signif.typeRm,mona.signif)

#drop trait *
mona.signif.traitRm <- glmer(y.mona~im_type+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.signif.traitRm)
qqp(residuals(mona.signif.traitRm),"binom",size=size.mona,prob=0.5)
anova(mona.signif.traitRm,mona.signif)

#add type*trait interaction
mona.signif.interact <- glmer(y.mona~im_type*im_trait+(1|group/subject/trial_subject),data=dat.mona,family=binomial)
plot(mona.signif.interact)
qqp(residuals(mona.signif.interact),"binom",size=size.mona,prob=0.5)
anova(mona.signif.interact,mona.null)
summary(mona.signif.interact)

#final model
mona.final <- mona.signif.interact
anova(mona.final,mona.null)
summary(mona.final) #significant bias for conspecifics, non-shared traits, interaction

exp(coef(summary(mona.final))[ , "Estimate"]) #odds ratio
1/exp(coef(summary(mona.final))[ , "Estimate"])

#look at only con v. hetero w/o eyebrow patches
dat.mona.t1 <- dat.mona[dat.mona$condition=='not_shared',]
y.mona.t1 <- cbind(dat.mona.t1$looks_5s,dat.mona.t1$total_looks_5s-dat.mona.t1$looks_5s)

mona.null.t1 <- glmer(y.mona.t1~1+(1|group/subject/trial_subject),data=dat.mona.t1,family=binomial)
plot(mona.null.t1)
qqp(residuals(mona.null.t1),"binom",size=size.mona,prob=0.5)

mona.type.t1 <- glmer(y.mona.t1~im_type+(1|group/subject/trial_subject),data=dat.mona.t1,family=binomial) #note: type & trait are the same here
plot(mona.type.t1)
qqp(residuals(mona.type.t1),"binom",size=size.mona,prob=0.5)
anova(mona.type.t1,mona.null.t1)
summary(mona.type.t1) #significant bias for shared trait

exp(coef(summary(mona.type.t1))[ , "Estimate"]) #odds ratio
1/exp(coef(summary(mona.type.t1))[ , "Estimate"]) 

