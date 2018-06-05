library(ggplot2)
library(cowplot)
library(plyr)

# Set working directory
setwd("C:/Users/Clif/Dropbox/Completed manuscripts/Flea feeding")

# Set dodge width for jittering
dodge = position_dodge(width=1.5)

# Import flea testing qPCR data
qPCR.dat = read.csv("2016-06-08_Flea_feeding_summary_qPCR.csv",header=T)

# Calculate mean and median Ct values
inf.Ct = data.frame(ddply(qPCR.dat, ~Day+logFeed, summarise,
                              mean=mean(Ct), median=median(Ct)))

pdf('Flea_testing_Ct.pdf', width=11, height=8.5)
# Plot flea testing Ct values
ggplot(data=qPCR.dat, aes(x=Day, y=Ct, 
                          group=factor(logFeed), colour=factor(logFeed))) +
  geom_line(data=inf.Ct, aes(x=Day, y=mean), size=1.5) +
  geom_jitter(size=4, alpha=0.8, width=0.8) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  scale_x_continuous(name="DPI", breaks=c(0, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Cycle threshold", breaks=seq(30, 45, 5), limits=c(30, 45))
dev.off()

# Calculate mean and median Ct values
inf.copies = data.frame(ddply(qPCR.dat, ~Day+logFeed, summarise,
                          mean=mean(logCopies), median=median(logCopies)))

# Plot flea testing log copy numbers
logCopies = ggplot(data=qPCR.dat, aes(x=Day, y=logCopies,
                                      group=factor(logFeed), colour=factor(logFeed))) +
  geom_line(data=inf.copies, aes(x=Day, y=mean), size=1.5) +
  geom_jitter(size=4, alpha=0.8, width=0.8) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  scale_x_continuous(name="DPI", breaks=c(0, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Log copy number", breaks=seq(0, 5, 1), limits=c(0, 5))

# Differences in infection
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==0),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==0),]$logCopies) # 0 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==0),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==0),]$logCopies) # 0 DPI
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==5),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==5),]$logCopies) # 5 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==5),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==5),]$logCopies, var.equal=F) # 5 DPI
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==7),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==7),]$logCopies) # 7 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==7),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==7),]$logCopies) # 7 DPI
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==9),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==9),]$logCopies) # 9 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==9),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==9),]$logCopies) # 9 DPI
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==11),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==11),]$logCopies) # 11 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==11),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==11),]$logCopies) # 11 DPI
var.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==13),]$logCopies,
         qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==13),]$logCopies) # 13 DPI
t.test(qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day==13),]$logCopies,
       qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day==13),]$logCopies) # 13 DPI

# Calculate decline in infection after 5 DPI
med.fleas = qPCR.dat[which(qPCR.dat$logFeed==7.77&qPCR.dat$Day>0),]
lm.med.fleas = lm(med.fleas$logCopies~med.fleas$Day)
summary(lm.med.fleas)

hi.fleas = qPCR.dat[which(qPCR.dat$logFeed==9.09&qPCR.dat$Day>0),]
lm.hi.fleas = lm(hi.fleas$logCopies~hi.fleas$Day)
summary(lm.hi.fleas)

# Import flea testing prevalence data, Clopper-Pearson exact confidence intervals
PCR.CI = read.csv("2016-12-28_Flea_feeding_PCR_CI.csv",header=T)

# Plot flea testing prevalence data
flea.prev = ggplot(data=PCR.CI, aes(x=Day, y=p, ymin=lower, ymax=upper, 
                                    group=factor(logFeed), colour=factor(logFeed),
                                    label=Tested)) +
  geom_linerange(position=dodge, size=1.5, alpha=0.2) +
  geom_line(position=dodge, size=1.5) +
  geom_point(stat="identity", position=dodge, size=4, alpha=0.8) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  guides(colour=guide_legend(title="Log feed")) +
  scale_x_continuous(name="DPI", breaks=c(0, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Proportion positive", breaks=seq(0, 1, 0.2)) +
  geom_text(aes(x=Day, y=-0.05), position=dodge)

pdf('flea_testing_prev_and_copies.pdf', height=8.5, width=11)
plot_grid(flea.prev, logCopies, labels=c('a','b'), ncol=1, nrow=2)
dev.off()

# Calculate the decline in prevalence after 5 DPI
med.prev = PCR.CI[which(PCR.CI$logFeed==7.77&PCR.CI$Day>0),]
lm.med.prev = lm(med.prev$p~med.prev$Day)
summary(lm.med.prev)

hi.prev = PCR.CI[which(PCR.CI$logFeed==9.09&PCR.CI$Day>0),]
lm.hi.prev = lm(hi.prev$p~hi.prev$Day)
summary(lm.hi.prev)

# Differences in prevalence between 0 and 5 DPI
prop.test(PCR.CI$x[c(2, 5)], PCR.CI$n[c(2, 5)]) # medium group
prop.test(PCR.CI$x[c(3, 6)], PCR.CI$n[c(3, 6)]) # high group

# Differences in prevalence between groups
prop.test(PCR.CI$x[2:3], PCR.CI$n[2:3]) # 0 DPI
prop.test(PCR.CI$x[5:6], PCR.CI$n[5:6]) # 5 DPI
prop.test(PCR.CI$x[8:9], PCR.CI$n[8:9]) # 7 DPI
prop.test(PCR.CI$x[11:12], PCR.CI$n[11:12]) # 9 DPI
prop.test(PCR.CI$x[14:15], PCR.CI$n[14:15]) # 11 DPI
prop.test(PCR.CI$x[17:18], PCR.CI$n[17:18]) # 13 DPI

# Import flea feces testing qPCR data
FqPCR.dat = read.csv("2016-12-28_Flea_feces_qPCR.csv",header=T)

# Calculate mean and median Ct values
flea.Ct = data.frame(ddply(FqPCR.dat, ~Day+logFeed, summarise,
                           mean=mean(Ct), median=median(Ct)))

pdf("Flea_feces_Ct.pdf",width=11,height=8.5)
# Plot flea feces testing Ct values
ggplot(data=FqPCR.dat, aes(x=Day, y=Ct, 
                          group=factor(logFeed), colour=factor(logFeed))) +
  geom_line(data=flea.Ct, aes(x=Day, y=mean), size=1.5) +
  geom_jitter(size=4, alpha=0.8, width=0.5) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  scale_x_continuous(name="DPI", breaks=c(0, 3, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Cycle threshold", breaks=seq(30, 45, 5), limits=c(30, 45))
dev.off()

# Calculate mean and median Ct values
flea.copies = data.frame(ddply(FqPCR.dat, ~Day+logFeed, summarise,
                           mean=mean(logCopies), median=median(logCopies)))
                         
pdf("Flea_feces_logCopies.pdf",width=11,height=8.5)
# Plot flea feces testing log copy numbers
ggplot(data=FqPCR.dat, aes(x=Day, y=logCopies,
                          group=factor(logFeed), colour=factor(logFeed))) +
  geom_line(data=flea.copies, aes(x=Day, y=mean), size=1.5) +
  geom_jitter(size=4, alpha=0.8, width=0.8) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  scale_x_continuous(name="DPI", breaks=c(0, 3, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Log copy number", breaks=seq(0, 5, 1), limits=c(0, 5))
dev.off()

# Calculate decline in infection after 5 DPI
med.feces = FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day>3),]
lm.med.feces = lm(med.feces$logCopies~med.feces$Day)
summary(lm.med.feces)

hi.feces = FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day>0),]
lm.hi.feces = lm(hi.feces$logCopies~hi.feces$Day)
summary(lm.hi.feces)

# Differences in infection
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==0),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==0),]$logCopies) # 0 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==3),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==3),]$logCopies) # 3 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==5),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==5),]$logCopies) # 5 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==7),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==7),]$logCopies) # 7 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==9),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==9),]$logCopies) # 9 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==11),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==11),]$logCopies) # 11 DPI
t.test(FqPCR.dat[which(FqPCR.dat$logFeed==7.77&FqPCR.dat$Day==13),]$logCopies,
       FqPCR.dat[which(FqPCR.dat$logFeed==9.09&FqPCR.dat$Day==13),]$logCopies) # 13 DPI

# Import flea feces testing prevalence data, Clopper-Pearson exact confidence intervals
FPCR.CI = read.csv("2016-12-28_Flea_feces_PCR_CI.csv",header=T)

pdf("Flea_feces_prevalence.pdf",width=11,height=8.5)
# Plot flea feces testing prevalence data
ggplot(data=FPCR.CI, aes(x=Day, y=p, ymin=lower, ymax=upper, 
                        group=factor(logFeed), colour=factor(logFeed),
                        label=Tested)) +
  geom_linerange(position=dodge, size=1.5, alpha=0.2) +
  geom_line(position=dodge, size=1.5) +
  geom_point(stat="identity", position=dodge, size=4, alpha=0.8) +
  theme_classic(base_size=20) +
  scale_colour_manual(values=c("#000000", "#E69F00", "#009E73"), name="Log feed") +
  guides(colour=guide_legend(title="Log feed")) +
  scale_x_continuous(name="DPI", breaks=c(0, 3, 5, 7, 9, 11, 13)) +
  scale_y_continuous(name="Proportion positive", breaks=seq(0, 1, 0.2)) +
  geom_text(aes(x=Day, y=-0.05), position=dodge)
dev.off()
