##### Script for handling of post-module information

Exp_mods <- as.data.frame(MEs)
rownames(Exp_mods) <- rownames(datExpr)

### Do tertiles for RED and Red Modules ############################

### for 4 quartiles

brks <- with(Exp_mods, quantile(MEbrown, probs = c(0, 0.25, 0.5, 0.75, 1)))
BrownV <- within(Exp_mods, quartile <- cut(MEbrown, breaks = brks, labels = 1:4, 
                                           include.lowest = TRUE))

brks <- with(Exp_mods, quantile(MEbrown, probs = c(0, 0.33, 0.66, 1)))
BrownV <- within(Exp_mods, tertile <- cut(MEbrown, breaks = brks, labels = 1:3, 
                                           include.lowest = TRUE))

brks <- with(Exp_mods, quantile(MEred, probs = c(0, 0.25, 0.5, 0.75, 1)))
RedV <- within(Exp_mods, quartile <- cut(MEred, breaks = brks, labels = 1:4, 
                                                   include.lowest = TRUE))
brks <- with(Exp_mods, quantile(MEred, probs = c(0, 0.33, 0.66, 1)))
RedV <- within(Exp_mods, tertile <- cut(MEred, breaks = brks, labels = 1:3, 
                                          include.lowest = TRUE))

### add the quartiles into the Exp_mods

Exp_mods$Brown.quar <- BrownV$quartile
Exp_mods$Red.quar <- RedV$quartile

Exp_mods$Brown.tert <- BrownV$tertile
Exp_mods$Red.tert <- RedV$tertile

### Add in the clinical data for use in PRISM

Exp_mods$Pam50 <- datClin$PAM50Call_RNAseq
Exp_mods$Survival_time <- datClin$OS_Time_nature2012
Exp_mods$Survival_event <- datClin$OS_event_nature2012

## Save out the file

write.table(Exp_mods, "Module_expression_set.txt", sep = "\t")

##########
library(plyr) ## for the "count" function, which is sweeeet!
count(Values$quartile)
## you can see that the distribution is very even / same kinda numbers in each quartile!



#### and you're off and running!!

#### Set the MethScore quartiles as a factor

Brown <- as.factor(Exp_mods$Brown.quar)
Red <- as.factor(Exp_mods$Red.quar)

Brown3 <- as.factor(Exp_mods$Brown.tert)
Red3 <- as.factor(Exp_mods$Red.tert)

Pam50 <- as.factor(datClin$PAM50Call_RNAseq)


boxplot(Exp_mods$MEbrown ~ Pam50)
boxplot(Exp_mods$MEred ~ Pam50)
boxplot(Exp_mods$MEcyan~ Pam50)

### need to do AOV and Tukeys on the brown and red


### Survivial curves
library(survival)

fit.diff = survdiff(Surv(OS_Time_nature2012,OS_event_nature2012 == 1) ~ factor(Brown3),
                    data=datClin) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Brown3)))-1),3) 
fit1 = survfit(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~Brown3,
               data=datClin,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n Brown Module Expression"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Brown3)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(Brown), data = datClin))

######## RED Mod

fit.diff = survdiff(Surv(OS_Time_nature2012,OS_event_nature2012 == 1) ~ factor(Red3),
                    data=datClin) 
chisq2 = signif(1-pchisq(fit.diff$chisq,length(levels(factor(Red3)))-1),3) 
fit1 = survfit(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~Red3,
               data=datClin,conf.type="log-log") 
kmcolours <- c("black", "red", "green", "blue")
plot(fit1, conf.int=F,col=kmcolours,xlab="Time to death (days)", 
     ylab="Survival",main=c("All subtype survival by\n Red Module Expression"), 
     lwd=4,mark.time=TRUE)
legend("bottomleft",legend=levels(factor(Red3)),
       fill = kmcolours, cex = 1)

###do another anova for signifigance
anova(coxph(Surv(OS_Time_nature2012,OS_event_nature2012 == 1)~factor(Red), data = datClin))