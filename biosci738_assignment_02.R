data <- readxl::read_xlsx("heather_all.xlsx", skip = 2, sheet = 3)
data

library(tidyverse)
data <- data %>%
  select(c("Block!", "Treatment!","Year!", "Total natives")) %>%
  rename("block" = "Block!", "treatment" = "Treatment!", "year" = "Year!",
         "natives" = "Total natives") %>%
  mutate(treatment = as_factor(treatment))


## making the violin plots 
library(ggplot2)
library(ggpubr)


# Basic violin plot per year 
#2008
violplot_eight <- ggplot(subset(data, year %in% c("2008")))+ 
                           geom_violin(aes(x=as.factor("2008"), y=natives, fill=treatment))+
  labs(x="", y="total natives", title="1. Data Spread for 2008")
violplot_eight

#2009
violplot_nine <- ggplot(subset(data, year %in% c("2009")))+ 
  geom_violin(aes(x=as.factor("2009"), y=natives, fill=treatment))+
  labs(x="", y="total natives", title="2. Data Spread for 2009")+
  coord_cartesian(ylim = c(4, 15))
violplot_nine

#2010
violplot_ten <- ggplot(subset(data, year %in% c("2010")))+ 
  geom_violin(aes(x=as.factor("2010"), y=natives, fill=treatment))+
  labs(x="", y="total natives", title="3. Data Spread for 2010")+
  coord_cartesian(ylim = c(4, 15))
violplot_ten

#2011
violplot_twelve <- ggplot(subset(data, year %in% c("2012")))+ 
  geom_violin(aes(x=as.factor("2012"), y=natives, fill=treatment))+
  labs(x="", y="total natives", title="4. Data Spread for 2012")
violplot_twelve

#merge the violin plots
ggarrange(violplot_eight, violplot_nine, violplot_ten, violplot_twelve,
          legend = c("right"),
          legend.grob = get_legend(violplot_eight, position = c("right")))


## making my linear model - this is a one-way model
oneway_model <- lm(natives ~ treatment, data = data)  
summary(oneway_model$coefficients, digits = 2)
gglm::gglm(oneway_model)

with(data, tapply(natives, treatment, mean))

## pairwise comparison of means 
library("predictmeans")
pm <- predictmeans::predictmeans(oneway_model, modelterm = "treatment", pairwise = TRUE, plot = FALSE)

my_pm_comparisons_function <- function(pm){
  ## treatment names
  nms <- colnames(pm$`Pairwise LSDs`)
  ## comparison names
  c_names <- outer(nms, nms, function(x, y) paste(x, "-", y, sep = ""))
  ## unique comparison names
  names <- c_names[upper.tri(c_names)]
  ## extract bit of predictmeant outpur
  LSD_matrix <- pm[[5]]
  pwise_p <- pm$`Pairwise p-value`
  LSD <- LSD_matrix[lower.tri(LSD_matrix)]
  diffs <- LSD_matrix[upper.tri(t(LSD_matrix))]
  tstats <- pwise_p[upper.tri(pwise_p)]
  pvalues <- pwise_p[lower.tri(pwise_p)]
  SED <- pm$`Standard Error of Differences`
  L95    <- diffs-LSD
  U95    <- diffs+LSD
  ## create dataframe with needed elements
  results <- data.frame(differences = diffs, 
                        LSD = LSD, 
                        tstatistic = tstats, 
                        SED = SED, 
                        Lower95CI = L95,
                        Upper95CI = U95,
                        pvalue = pvalues)
  ## append comparison names as rownames
  rownames(results) <- names
  ## return data frame
  return(results)
}
my_pm_comparisons_function(pm)

pm_plot <- predictmeans::predictmeans(oneway_model, modelterm = "treatment", pairwise = TRUE)
pm_plot$predictmeansPlot #produces a plot of predicted means and LSD
pm$`Predicted Means`

## fitting two interaction models where year is treated as a factor 
#1. treatment*year 
twoway_model <- lm(natives ~ treatment*year, data = data)
anova(twoway_model)

#2. year*treatment
twoway_year <- lm(natives ~ year*treatment, data = data)
anova(twoway_year)

#type II anova on the two way model
car::Anova(twoway_model, type = 2)

#creating four models for multi level marketing i mean linear mixed-effect models
lmer_one <- lme4::lmer(natives ~ treatment + (1|block), data = data, REML = FALSE)
summary(lmer_one)
plot(lmer_one)

lmer_two <- lme4::lmer(natives ~ treatment + year + (1|block), data = data)
summary(lmer_two)

lmer_three <- lme4::lmer(natives ~ treatment:year + (1|block), data = data)
summary(lmer_three)

lmer_four <- lme4::lmer(natives ~ treatment*year + (1|block), data = data)
summary(lmer_four)

#anova-ing the above models 
anova(lmer_one, lmer_two, lmer_three, lmer_four)


##use predictmeans to do pairwise on model i. 
pm1 <- predictmeans::predictmeans(lmer_one, modelterm = "treatment", pairwise = TRUE, plot = FALSE)
my_pm_comparisons_function(pm1)
print(pm1)

#subset the herbicide v herbcide and biocontorl
herbcide_hb_data <- data[!(data$treatment =="Control" | data$treatment=="Biocontrol"),]%>%
  mutate(year = as_factor(year))
herbcide_hb_data

#mdoel i for subset data 
lmer_herbicide <- lme4::lmer(natives ~ treatment + (1|block), data = herbcide_hb_data)


#do a bonferroni's adjustment 
library("predictmeans")
alpha.adj <- 0.05/choose(2,2)
bonferroni <-  predictmeans::predictmeans(lmer_herbicide, 
                                          modelterm = "treatment", 
                                          adj = "bonferroni",
                                          level = alpha.adj,
                                          pairwise = TRUE, 
                                          plot = FALSE)


my_pm_comparisons_function(bonferroni)





#tukeys HSD
aov_mod <- aov(natives ~ treatment + (1|block), data = data)
summary(aov_mod)

qtukey(p = 1- 0.05, nmeans = 4, df = 87)
x
sqrt(2 * 6.214 / 6)

HSD <- (qtukey(p = 1- 0.05, nmeans = 4, df = 87)/(sqrt(2)))*1.439213

HSD

tukey_one <- TukeyHSD(aov_mod)
tukey_two <- HSD.test(aov_mod, trt='treatment')
summary(tukey_two)

#fitting the tukeys HSD into the pairwise function 
tukey <- predictmeans::predictmeans(lmer_one,
                                   modelterm = "treatment", 
                                   adj = "tukey",
                                   level = alpha.adj,
                                   pairwise = TRUE, 
                                   plot = FALSE)

all_diffs <- outer(tukey$`Predicted Means`, tukey$`Predicted Means`, "-")
comparison_table <- data.frame(differences = all_diffs[lower.tri(all_diffs)]) %>%
  mutate(upper = differences + HSD,
         lower = differences - HSD,
         HSD = HSD)
comparison_table

tukeyhsd <- TukeyHSD(aov(natives~treatment+(1|block), data = data))

print(tukey$`Pairwise p-value`)

em <- emmeans::emmeans(lmer_one, specs =  "treatment")
pairs(em, adjust = "tukey")
