#Formats looking time data, generating summary data from raw coded looks
#(Currently specific to CERCOPAN guenon project)

#INPUT: 
# dat: data frame of raw looking time data (required)
# secs: seconds of looking to include for each subject; 'all' does not restrict (default = 5)
# exclude: list of subjects to exclude from returned analyses (default = Mimi's Baby, Nacera, Happiness, Maya, Wizkid, Bobierto)
    # Mimi's Baby, Nacera, Wizkid, Bobierto - juveniles
    # Happiness, Maya - didn't participate in all trials

#OUTPUT: 
# data frame of formatted looking time data

#NOTE: left and right are relative, i.e. from perspective of stimuli, im1 is on left & im2 is on right, from perspective of 
#monkey, left is bias toward im2 and right is bias toward im1


formatLTdata <- function(dat,secs=5,exclude=c("mimisbaby","nacera","happiness","maya","wizkid","bobierto","?")) {
  library(tidyr)
  library(dplyr)
  
  #import trial & subject data
  trial_info<- read.csv('trial info.csv',header=T)
  subject_info <- read.csv('subject info.csv',header=T)
  stimulus_info <- read.csv('stimulus info.csv',header=T)
  
  subject_info$subject <- tolower(subject_info$subject)
  subject_info$subject <- gsub('\'','',subject_info$subject)
  subject_info$subject <- gsub('-','',subject_info$subject)
  subject_info$subject <- gsub(' ','',subject_info$subject)
  
  trial_info$familiarity_im1 <- ordered(trial_info$familiarity_im1,levels=c('not_present','not_visible','visible'))
  trial_info$familiarity_im2 <- ordered(trial_info$familiarity_im2,levels=c('not_present','not_visible','visible'))
  trial_info$trial_order <- as.ordered(trial_info$trial_order)

  trial_info$im1_eyeContact <- with(stimulus_info,eyeContact[match(trial_info$im1,image)])
  trial_info$im2_eyeContact <- with(stimulus_info,eyeContact[match(trial_info$im2,image)])
  
  #create data frame (wide style)
  dat$subject <- tolower(dat$subject)
  dat$subject <- gsub('\'','',dat$subject)
  dat$subject <- gsub('-','',dat$subject)
  dat$subject <- gsub(' ','',dat$subject)
  
  dat.wide <- merge(subject_info,trial_info,by='group')

  dat.wide$trial_subject <- mapply(paste,dat.wide$species,sprintf("%02d",dat.wide$trial),dat.wide$subject,MoreArgs=list(sep='.'))
  dat$trial_subject <- mapply(paste,dat$species,sprintf("%02d",dat$trial),dat$subject,MoreArgs=list(sep='.'))
  
  
  #calculate R & L looks per subject per trial for all looks
  for (i in 1:dim(dat.wide)[1]) {
    dat.wide$L_all[i] = sum(dat$duration[(dat$trial_subject==dat.wide$trial_subject[i] & dat$look=='L')])
    dat.wide$R_all[i] = sum(dat$duration[(dat$trial_subject==dat.wide$trial_subject[i] & dat$look=='R')])
  }
  
  #calculate total looks
  dat.wide$total_looks <- dat.wide$L_all+dat.wide$R_all
  
  #calculate L & R looks per subject per trial for given secs
  if (secs!='all') {
    dat.wide <- cbind(dat.wide,dat.wide$L_all)
    colnames(dat.wide)[length(colnames(dat.wide))] <- paste('L_',secs,'s',sep='')
    dat.wide <- cbind(dat.wide,dat.wide$R_all)
    colnames(dat.wide)[length(colnames(dat.wide))] <- paste('R_',secs,'s',sep='')

    nframes <- secs*30 #30 fps video
    buffer <- 30 #keep going until 1 sec of non-looking
    for (i in 1:dim(dat.wide)[1]) {
      if (dat.wide$total_looks[i]>nframes) {
        tmp <- dat[dat$trial==dat.wide[i,]$trial & dat$subject==as.character(dat.wide[i,]$subject),]
        tmp <- tmp[tmp$look=='R' | tmp$look=='L',]
        line <- 1
        final <- 0
        R <- 0
        L <- 0
        while (final==0) { # (R+L<nframes)
          #sum R & L durations
          if (tmp$look[line]=='R') {
            R <- R + tmp$duration[line]
          } else if (tmp$look[line]=='L') {
            L <- L + tmp$duration[line]
          }
          
          #check to see if nframes is over threshold & current look is followed by 1 sec of non-looking
          # if (R+L>=nframes) { #( line==dim(tmp)[1] || ((R+L>=nframes) && ((tmp$first[line+1]-tmp$last[line])>=buffer)) )
          if ( line+1>dim(tmp)[1] || ( R+L>=nframes && ((tmp$first[line+1]-tmp$last[line])>=buffer) ) ) {
            final <- 1
          }
          
          line <- line + 1
        }
        dat.wide[i,paste('L_',secs,'s',sep='')] <- L
        dat.wide[i,paste('R_',secs,'s',sep='')] <- R
      } #end if
    } #end for
  
  dat.wide <- cbind(dat.wide,dat.wide[,paste('L_',secs,'s',sep='')]+dat.wide[,paste('R_',secs,'s',sep='')])
  colnames(dat.wide)[length(colnames(dat.wide))] <- paste('total_looks_',secs,'s',sep='')
  } #end if
  
  #remove excluded subjects
  dat.wide <- dat.wide[!dat.wide$subject %in% exclude,]
  
  #convert to long format, log transform looks, & get image variables
  if (secs=='all') {
    dat.long <- gather(dat.wide,pres_spot,looks,L_all,R_all)
  } else {
    dat.long <- gather_(dat.wide,'pres_spot','looks',c(paste('L_',secs,'s',sep=''),paste('R_',secs,'s',sep='')))
  }
  
  dat.long$pres_spot <- ifelse(substr(dat.long$pres_spot,1,1)=='L',2,1)
  dat.long$pres_spot <- as.factor(dat.long$pres_spot)
  
  colnames(dat.long)[colnames(dat.long)=='sex'] <- 'subject_sex'
  
  dat.long$looks_prop <- dat.long$looks/dat.long$total_looks
  dat.long$im_spp <- ifelse(dat.long$pres_spot==1,as.character(dat.long$spp1),as.character(dat.long$spp2)) #the image for that line (right or left in trial)
  dat.long$im_type <- ifelse(dat.long$im_spp==dat.long$species,'conspecific','heterospecific') #whether the image for that line is a con or heterspp
  dat.long$im_type <- as.factor(dat.long$im_type)
  dat.long$stim_im <- ifelse(dat.long$pres_spot==1,as.character(dat.long$im1),as.character(dat.long$im2))
  dat.long$im_sex <- substr(as.character(sapply(dat.long$stim_im,function(x){strsplit(as.character(x),'_')}[[1]][2])),1,1)
  dat.long$im_familiarity <- ifelse(dat.long$pres_spot==1,as.character(dat.long$familiarity_im1),as.character(dat.long$familiarity_im2))
  dat.long$im_familiarity <- ordered(dat.long$im_familiarity,levels=c('not_present','not_visible','visible'))
  dat.long$eye_contact <- ifelse(dat.long$pres_spot==1,as.character(dat.long$im1_eyeContact),as.character(dat.long$im2_eyeContact)) #whether the stimulus image is making eye contact

  #specify whether image & paired image shares trait w/ subject
  dat.long$im_trait <- ''
  dat.long$pair_trait <- ''
  for (i in 1:dim(dat.long)[1]) {
    if (dat.long$im_type[i]=="conspecific") {
      dat.long$im_trait[i] <- ifelse(dat.long$conspp_mod[i]=="no","shared","not_shared")
      dat.long$pair_trait[i] <- ifelse(dat.long$heterospp_trait[i]=="shared","shared","not_shared")
    }
    else {
      dat.long$im_trait[i] <- ifelse(dat.long$heterospp_trait[i]=="shared","shared","not_shared")
      dat.long$pair_trait[i] <- ifelse(dat.long$conspp_mod[i]=="no","shared","not_shared")
    }
  }
  dat.long$im_trait <- as.factor(dat.long$im_trait)
  
  # dat.long$pres_spot <- gsub(paste('_',secs,'s',sep=''),'',dat.long$pres_spot)
  colnames(dat.long)[which(colnames(dat.long)=='looks')] <- paste('looks_',secs,'s',sep='')
  
  #only keep relevant variables in long format
  dat.long = dat.long[,colnames(dat.long) %in% c('group','species','subject','dob','age.years','subject_sex','origin','trial',
                                                 'trial_subject','n_subjects','condition','trial_order','pres_spot','pattern','icc',
                                                 'total_looks',paste('total_looks_',secs,'s',sep=''),paste('looks_',secs,'s',sep=''),
                                                 'looks_prop','im_spp','im_type','stim_im','im_sex','im_familiarity','eye_contact',
                                                 'im_trait','pair_trait')]
  
  return(list(dat.wide,dat.long))
}
