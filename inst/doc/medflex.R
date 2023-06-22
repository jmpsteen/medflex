### R code from vignette source 'medflex.Rnw'

###################################################
### code chunk number 1: medflex.Rnw:208-213
###################################################
options(prompt = "R> ", continue = "+  ", width = 76, digits = 3, useFancyQuotes = FALSE)
  listing <- function(x, options) {
    paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n",
      x, "\\end{lstlisting}\n", sep = "")
  }


###################################################
### code chunk number 2: medflex.Rnw:216-219
###################################################
library("medflex")
data("UPBdata")
head(UPBdata)


###################################################
### code chunk number 3: medflex.Rnw:292-294
###################################################
medFit <- glm(negaff ~ factor(attbin) + gender + educ + age, 
  family = gaussian, data = UPBdata)


###################################################
### code chunk number 4: medflex.Rnw:297-298
###################################################
expData <- neWeight(medFit)


###################################################
### code chunk number 5: medflex.Rnw:301-302
###################################################
head(expData, 4)


###################################################
### code chunk number 6: medflex.Rnw:306-308
###################################################
expData <- neWeight(negaff ~ factor(attbin) + gender + educ + age,
  data = UPBdata)


###################################################
### code chunk number 7: medflex.Rnw:312-314
###################################################
w <- weights(expData)
head(w, 10)


###################################################
### code chunk number 8: medflex.Rnw:352-355
###################################################
neMod1 <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age,
  family = binomial("logit"), expData = expData, se = "robust")
summary(neMod1)


###################################################
### code chunk number 9: medflex.Rnw:368-369
###################################################
exp(confint(neMod1)[c("attbin01", "attbin11"), ])


###################################################
### code chunk number 10: medflex.Rnw:402-404
###################################################
impFit <- glm(UPB ~ factor(attbin) + negaff + gender + educ + age,
  family = binomial("logit"), data = UPBdata)


###################################################
### code chunk number 11: medflex.Rnw:408-409
###################################################
expData <- neImpute(impFit)


###################################################
### code chunk number 12: medflex.Rnw:412-414
###################################################
expData <- neImpute(UPB ~ factor(attbin) + negaff + gender + educ + age,
  family = binomial("logit"), data = UPBdata)


###################################################
### code chunk number 13: medflex.Rnw:417-418
###################################################
head(expData, 4)


###################################################
### code chunk number 14: medflex.Rnw:423-425
###################################################
neMod1 <- neModel(UPB ~ attbin0 + attbin1 + gender + educ + age,
  family = binomial("logit"), expData = expData, se = "robust")


###################################################
### code chunk number 15: medflex.Rnw:428-429
###################################################
summary(neMod1)


###################################################
### code chunk number 16: medflex.Rnw:461-464
###################################################
expData <- neImpute(UPB ~ attcat + negaff + gender + educ + age,
  family = binomial("logit"), data = UPBdata)
head(expData)


###################################################
### code chunk number 17: medflex.Rnw:467-470
###################################################
neMod <- neModel(UPB ~ attcat0 + attcat1 + gender + educ + age,
  family = binomial("logit"), expData = expData, se = "robust")
summary(neMod)


###################################################
### code chunk number 18: medflex.Rnw:473-475
###################################################
library("car")
Anova(neMod)


###################################################
### code chunk number 19: medflex.Rnw:483-486
###################################################
expData <- neImpute(UPB ~ att + negaff + gender + educ + age, 
  family = binomial("logit"), data = UPBdata, nRep = 3)
head(expData)


###################################################
### code chunk number 20: medflex.Rnw:489-492
###################################################
neMod1 <- neModel(UPB ~ att0 + att1 + gender + educ + age,
  family = binomial("logit"), expData = expData, se = "robust")
summary(neMod1)


###################################################
### code chunk number 21: medflex.Rnw:535-540
###################################################
expData <- neImpute(UPB ~ att * negaff + gender + educ + age,
  family = binomial("logit"), data = UPBdata)
neMod2 <- neModel(UPB ~ att0 * att1 + gender + educ + age,
  family = binomial("logit"), expData = expData, se = "robust")
summary(neMod2)


###################################################
### code chunk number 22: medflex.Rnw:550-555
###################################################
impData <- neImpute(UPB ~ (att + negaff) * gender + educ + age,
  family = binomial("logit"), data = UPBdata)
neMod3 <- neModel(UPB ~ att0 + att1 * gender + educ + age,
  family = binomial("logit"), expData = impData, se = "robust")
summary(neMod3)


###################################################
### code chunk number 23: medflex.Rnw:563-567
###################################################
impData <- neImpute(UPB ~ (att + negaff) * educ + gender + age,
  family = binomial("logit"), data = UPBdata)
neMod4 <- neModel(UPB ~ (att0 + att1) * educ + gender + age,
  family = binomial("logit"), expData = impData, se = "robust")


###################################################
### code chunk number 24: medflex.Rnw:578-580
###################################################
lht <- neLht(neMod2, linfct = c("att0 + att0:att1 = 0",
  "att1 + att0:att1 = 0", "att0 + att1 + att0:att1 = 0"))


###################################################
### code chunk number 25: medflex.Rnw:583-584
###################################################
exp(cbind(coef(lht), confint(lht)))


###################################################
### code chunk number 26: medflex.Rnw:587-588
###################################################
summary(lht)


###################################################
### code chunk number 27: medflex.Rnw:594-596
###################################################
effdecomp <- neEffdecomp(neMod2)
summary(effdecomp)


###################################################
### code chunk number 28: medflex.Rnw:622-624
###################################################
neEffdecomp(neMod3)
neEffdecomp(neMod3, covLev = c(gender = "M"))


###################################################
### code chunk number 29: medflex.Rnw:629-631
###################################################
modmed <- neLht(neMod4, linfct = c("att1:educM = 0", "att1:educH = 0"))
summary(modmed, test = Chisqtest())


###################################################
### code chunk number 30: plot2
###################################################
par(mfrow = c(1, 2))
plot(neMod2, xlab = "log odds ratio")
plot(neMod2, xlab = "odds ratio", transf = exp)


###################################################
### code chunk number 31: fig2
###################################################
par(mfrow = c(1, 2))
plot(neMod2, xlab = "log odds ratio")
plot(neMod2, xlab = "odds ratio", transf = exp)


###################################################
### code chunk number 32: medflex.Rnw:651-654 (eval = FALSE)
###################################################
## par(mfrow = c(1, 2))
## plot(neMod2, xlab = "log odds ratio")
## plot(neMod2, xlab = "odds ratio", transf = exp)


###################################################
### code chunk number 33: medflex.Rnw:661-662
###################################################
expFit <- glm(att ~ gender + educ + age, data = UPBdata)


###################################################
### code chunk number 34: medflex.Rnw:665-667
###################################################
impData <- neImpute(UPB ~ att + negaff + gender + educ + age,
  family = binomial("logit"), data = UPBdata)


###################################################
### code chunk number 35: medflex.Rnw:674-677
###################################################
neMod5 <- neModel(UPB ~ att0 + att1, family = binomial("logit"),
  expData = impData, xFit = expFit, se = "robust")
summary(neMod5)


###################################################
### code chunk number 36: medflex.Rnw:680-682
###################################################
de <- neLht(neMod5, c("att0 = 0"))
ie <- neLht(neMod5, c("att1 = 0"))


###################################################
### code chunk number 37: medflex.Rnw:742-744
###################################################
impData <- neImpute(UPB ~ att + initiator * negaff + gender + educ + age,
  family = binomial("logit"), nMed = 2, data = UPBdata)


###################################################
### code chunk number 38: medflex.Rnw:747-750
###################################################
neMod6 <- neModel(UPB ~ att0 + att1 + gender + educ + age,
  family = binomial("logit"), expData = impData, se = "robust")
summary(neMod6)


