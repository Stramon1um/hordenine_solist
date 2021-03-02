library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)
library(kableExtra)
library(agricolae)

#Mac
setwd("~/Dropbox/Dottorato Mauro Maver/papers/Hordenina/dati_csv")

#
##
###
#### TIME COURSE 24H

TC_24<-read.csv('TC_24h.csv',header = TRUE, sep=";")

sapply(TC_24, class)

TC_24$time <- as.factor(TC_24$time) #o evito questa conversione oppure lascio GROUP=1 in ggplot

TC_24$molecule_order = factor(TC_24$molecule, levels = c("tyramine", "nmt", "hordenine"))

ggplot(TC_24, aes(time, avg, group=1))+
  geom_point(aes(fill=molecule, color=molecule),shape=21, size=3)+
  geom_line()+
  geom_errorbar(aes(x=time, ymin=avg-se, ymax=avg+se),
                width=.2,             # Width of the error bars
                position=position_dodge(0.9))+
  geom_smooth(method = lm,se=TRUE, aes(group=1),color="black")+
  facet_wrap(~molecule_order, nrow = 3, scales = "free_y")+
  ggtitle("Time course 24h")+
  ylab(bquote(''*mu~ 'mol' ~g^-1*'FW'))+
  xlab("hour")+
  theme_bw()+
  theme(legend.position = "none")


#
##
###
#### TIME COURSE 8

TC_8<-read.csv('TC_8d.csv',header = TRUE, sep=";")

sapply(TC_8, class)

TC_8$time <- as.factor(TC_8$time) #o evito questa conversione oppure lascio GROUP=1 in ggplot

TC_8$molecule_order = factor(TC_8$molecule, levels = c("tyramine", "nmt", "hordenine"))

ggplot(TC_8, aes(time, avg, group=1))+
  geom_point(aes(fill=molecule, color=molecule),shape=21, size=3)+
  geom_line()+
  geom_errorbar(aes(x=time, ymin=avg-se, ymax=avg+se),
                width=.2,             # Width of the error bars
                position=position_dodge(0.9))+
  geom_smooth(method = lm,se=TRUE, aes(group=1),color="black")+
  facet_wrap(~molecule_order, nrow = 3, scales = "free_y")+
  ggtitle("Time course 8 days")+
  ylab(bquote(''*mu~ 'mol' ~g^-1*'FW'))+
  xlab("hours")+
  theme_bw()+
  theme(legend.position = "none")

#
##
###
#### TIME COURSE CARENZE

TC_carenze<-read.csv('TC_carenze.csv',header = TRUE, sep=";")

sapply(TC_carenze, class)

#TC_carenze$time <- as.factor(TC_carenze$time) #o evito questa conversione oppure lascio GROUP=1 in ggplot

TC_carenze$molecule_order = factor(TC_carenze$molecule, levels = c("tyramine", "NMT", "hordenine"))

ggplot(TC_carenze, aes(time, avg))+
  geom_point(aes(fill=treatment, color=treatment),shape=21, size=3)+
  geom_line(aes(color=treatment))+
  geom_errorbar(aes(x=time, ymin=avg-se, ymax=avg+se),
                width=.2,             # Width of the error bars
                position=position_dodge(0.9))+
  #geom_smooth(method = lm,se=TRUE, aes(group=1),color="black")+
  facet_grid(col=vars(molecule_order), space = "free")+
  theme_bw()+
  theme(legend.position = "none")


ggplot(TC_carenze, aes(time, avg))+
  geom_point(aes(fill=treatment, color=treatment),shape=21, size=3)+
  #guides(fill=guide_legend(override.aes = list(size=2)))+
  geom_line(aes(color=treatment))+
  geom_errorbar(aes(x=time, ymin=avg-se, ymax=avg+se),
                width=.2,             # Width of the error bars
                position=position_dodge(0.9))+
  facet_grid(rows = vars(molecule_order), scales = "free_y")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'))+
  ggtitle("Time course Nutrient deficiencies")+
  ylab(bquote(''*mu~ 'mol' ~g^-1*'FW'))+
  xlab("hours")+
  theme_bw()+
  theme(legend.position = "top",
        legend.key.size = unit(3,"line"),
        legend.title = element_blank(),
        strip.text = element_text(size = 13, face = "bold"))


#
##
###
#### ADSORBIMENTO

adso<-read.csv('adsorbimento_quarryfield.csv',header = TRUE, sep=";")

sapply(adso, class)

names(adso)[1] <- "X" #Ce mg/L
names(adso)[2] <- "Y" #Qe mg/g 

Lang <- nls(formula = Y ~ Q*b*X/(1+b*X),  data = adso, start = list(Q = 300, b = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(Lang)

Qmax <- summary(Lang)$coefficients[1,1]
b    <- summary(Lang)$coefficients[2,1]

sigLang <- nls(formula = Y ~ Q*b*X/(1+b*X+s/X),  data = adso, start = list(Q = 10, b = 0.1, s = 10), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sigLang)

Qmax_sig <- summary(sigLang)$coefficients[1,1]
b_sig    <- summary(sigLang)$coefficients[2,1]
s_sig    <- summary(sigLang)$coefficients[3,1]


freu <- nls(formula = Y ~ K*(X)^(1/n),  data = adso, start = list(K = 300, n = 1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(freu)

Kf <- summary(freu)$coefficients[1,1]
n    <- summary(freu)$coefficients[2,1]

sips <- nls(formula = Y ~ (Ks*(X)^Bs)/(1+as*(X)^Bs),  data = adso, start = list(Ks = 29, Bs = 1, as = 0.1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sips)

Ks <- summary(sips)$coefficients[1,1]
Bs    <- summary(sips)$coefficients[2,1]
as    <- summary(sips)$coefficients[3,1]

sips2 <- nls(formula = Y ~ (qs*Ks*(X)^Bs)/(1+Ks*(X)^Bs),  data = adso, start = list(Ks = 0.1, Bs = 0.1, qs = 0.1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(sips2)

Ks_2 <- summary(sips2)$coefficients[1,1]
Bs_2    <- summary(sips2)$coefficients[2,1]
qs_2    <- summary(sips2)$coefficients[3,1]

langfreu <- nls(formula = Y ~ (qmax*b*(X)^n)/(1+(b*X)^n),  data = adso, start = list(qmax = 50, b = 0.01, n = -30), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(langfreu)

Qmax_lf <- summary(langfreu)$coefficients[1,1]
Klf    <- summary(langfreu)$coefficients[2,1]
nlf    <- summary(langfreu)$coefficients[3,1]

redlich <- nls(formula = Y ~ Kr*X/(1+ar*(X)^g),  data = adso, start = list(Kr = 1, g = -10, ar = 0.1), control = nls.control(maxiter = 1000,tol = 1e-05), algorith = "port")
summary(redlich)

Kr <- summary(redlich)$coefficients[1,1]
g    <- summary(redlich)$coefficients[2,1]
ar    <- summary(redlich)$coefficients[3,1]

f_freu <- function(x, Kf, n) {
  (Kf*(x)^(1/n))
}

f_lang <- function(x, Qmax, b) {
  (Qmax*b*x/(1+b*x))
}

f_siglang <- function(x, Qmax_sig, b_sig, s_sig) {
  (Qmax_sig*b_sig*x/(1+b_sig*x+s_sig/x))
}

f_sips <- function(x, Ks, Bs, as) {
  (Ks*(x)^Bs/(1+as*(x)^Bs))
}

f_redlich <- function(x, Kr, g, ar) {
  (Kr*x/(1+ar*(x)^g))
}

ggplot(data = adso, aes(x = X, y = Y))+
  geom_point(size = 4, shape = 18)+
  ggtitle("Non-linearized isotherm adsorption models - Hordenine in Quarryfield soil")+
  xlab("Ce (mg/L)")+
  ylab("qe (mg/g)")+
  geom_line(aes(color = "Freundlich"),size = 1, linetype = 1,stat="function", fun = function(x) f_freu(x, Kf=Kf, n=n))+
  geom_line(aes(color = "Sigmoidal Langmuir"),size = 1, linetype = 1,stat="function", fun = function(x) f_siglang(x, Qmax_sig=Qmax_sig, b_sig=b_sig, s_sig=s_sig))+
  geom_line(aes(color = "Langmuir"), size = 1, linetype = 1,stat="function", fun = function(x) f_lang(x, Qmax=Qmax, b=b))+
  geom_line(aes(color = "Sips"), size = 1, linetype = 1,stat="function", fun = function(x) f_sips(x, Ks=Ks, Bs=Bs, as=as))+
  geom_line(aes(color = "Redlich-Peterson"), size = 1, linetype = 1,stat="function", fun = function(x) f_redlich(x, Kr=Kr, g=g, ar=ar))+
  scale_color_manual(name="Model", values = c("red", "blue", "green", "purple", "black"))+
  annotate("label", size = 5, hjust = 0, fontface = 0, x = 2, y = max(adso$Y)-0.20, label = "Residual Sum of Squares (RSSs):\nFreundlich: 0.0581\nLangmuir: 0.0678\nSips: 0.0373\nRedlich-Peterson: 0.0561\nSigmoidal Langmuir: 0.0663")+ 
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        title = element_text(size = 15))


#
##
###
#### TC METABOLOMICA

TC_metabolomica<-read.csv('TC_metabolomica.csv',header = TRUE, sep=";")

#statistica per compound

plant.lm <- lm(area ~ treatment, data = subset(TC_metabolomica, compound=="tyramine"))
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test
tukey.test2 <- HSD.test(plant.av, trt = 'treatment')
tukey.test2

plant.lm <- lm(area ~ treatment, data = subset(TC_metabolomica, compound=="N-Methyltyramine"))
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test
tukey.test2 <- HSD.test(plant.av, trt = 'treatment')
tukey.test2

plant.lm <- lm(area ~ treatment, data = subset(TC_metabolomica, compound=="hordenine"))
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test
tukey.test2 <- HSD.test(plant.av, trt = 'treatment')
tukey.test2

plant.lm <- lm(area ~ treatment, data = subset(TC_metabolomica, compound=="hordenine"))
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test
tukey.test2 <- HSD.test(plant.av, trt = 'treatment')
tukey.test2




cdata_metab <- ddply(TC_metabolomica, c("treatment", "compound","t_test"), summarise,
                     N    = length(area),
                     mean = mean(area),
                     sd   = sd(area),
                     se   = sd / sqrt(N))

cdata_metab

cdata_metab$treatment_order=factor(cdata_metab$treatment, levels=c('C','N','S','Fe','P'))
cdata_metab$compound_order=factor(cdata_metab$compound, levels=c('tyramine','N-Methyltyramine','hordenine','candicine'))
ggplot(cdata_metab, aes(x=treatment, y=mean))+
  geom_bar(position=position_dodge2(), stat="identity")+
  scale_y_continuous(limits = c(0,60000000))+
  geom_text(aes(label=t_test, vjust=-3.5),size=5, show.legend = FALSE)+
  geom_errorbar(aes(x=treatment_order, ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(0.9))+
  facet_wrap(~compound_order, nrow = 2)+
  #facet_grid(rows = 1)+
  ggtitle ("Time course Targeted Metabolomica (Ultimo campionamento TC)")+
  ylab("Mean (Area)")+
  #ylim(0,15)+
  xlab("Deficiency")+
  theme_bw()+
  theme(strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=11))

#
##
###
#### TC FOLD CHANGE METABOLOMICA

TC_FC_metabolomica<-read.csv('Fold_change_metabolomica_carenze.csv',header = TRUE, sep=";")

sapply(TC_FC_metabolomica, class)

sum_FC_meta <- ddply(TC_FC_metabolomica, c("category","treatment"), summarise,
                    N = sum(!is.na(Log_FC)),
                    sum_FC = sum(Log_FC, na.rm=TRUE))

sum_FC_meta

ggplot(sum_FC_meta, aes(category, sum_FC, fill=treatment))+
  geom_bar(position=position_dodge2(), stat="identity")+
  ggtitle("Biosynthesis datatable - Sum FC ultimo campionamento TC")+
  xlab("Category")+
  ylab("Sum Fold Change")+
  theme_bw()+
  scale_fill_discrete(name = "   ")+
  coord_flip()

#
##
###
#### TABLE EXUDATES

tab_exu<-read.csv('TC_exudates.csv',header = TRUE, sep=";")

names(tab_exu) <- c("Days", "Tyramine (nmol g^-1^)", "N-methyltyramine (nmol g^-1^)","Hordenine (nmol g^-1^)")

kable(tab_exu, escape = FALSE) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 15) %>%
  row_spec(1:2, align = "c", color ="black")
