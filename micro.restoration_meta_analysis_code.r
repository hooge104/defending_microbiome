library('tidyverse')
library('metafor')
library("ggExtra")

#read in our dataset
dataset <- read_csv('Supplementary_Data_File_2._Studies_and_Effect_Sizes.csv')

#assign vi of zero the minimum value we found
dataset$vi[which(dataset$vi == 0)] <- min(filter(dataset, vi != 0)$vi)

#run the random effects model with metafor
modelOutput <- rma.mv(yi, vi, data = dataset, 
                      random = ~ 1 | DOI)

#results?
summary(modelOutput)

#we exponentiate the model results which are in log space
plotData <- data.frame('value' = (exp(modelOutput$beta)*100)-100,
                       'confLow' = (exp(modelOutput$ci.lb)*100)-100,
                       'confHigh' = (exp(modelOutput$ci.ub)*100)-100)

#turn the yi into % biomass increase
dataset$outcome <- (exp(dataset$yi) * 100)-100

#make a plot
metaPlot <- ggplot(dataset) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_jitter(aes(x = 1, y = outcome,
                  size = 1/vi, color = as.factor(outcome > 0),
                  fill = as.factor(outcome > 0),
                  shape = as.factor(outcome > 0)), 
              width = 0.1,
              height = 0,
              alpha = 0.3) +
  geom_point(aes(x = 1, y = value),
             size = 4, data = plotData,
             shape = 16, color = '#50725A') +
  geom_errorbar(aes(x = 1, ymin = 
                      confLow, ymax = confHigh),
                width = 0, data = plotData,
                color = '#50725A') +
  xlim(0.75,1.25) +
  coord_flip() +
  theme_minimal() +
  scale_color_manual(values = c('black', '#50725A')) +
  scale_fill_manual(values = c('white', '#50725A')) +
  scale_shape_manual(values = c(21, 16)) +
  ylab("Difference in biomass (%)") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")

finalPlot <- 
  ggMarginal(metaPlot, data=dataset, type = 'histogram',
             margins = 'x', size = 2, fill = '#50725A', xparams = list(binwidth = 20))

