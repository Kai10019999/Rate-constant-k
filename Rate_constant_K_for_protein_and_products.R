library(ggpubr)
library(ggplot2)
#install.packages("extrafont")
library(extrafont)
library(nls2)
library(lme4)
library(minpack.lm)


#For parent protein
Func_decrease <- function(data_df, title, interce) {
  func = nlsLM(y ~ a * exp(b * x),data=subset(data_df, y > -1),algorithm = "LM",start = list(a = 66.88, b = -0.026273))  
  d <- summary(func)  #initial value a = 66.88, b = -0.026273
  
  a_func <- d[[10]][[1]]
  b_func <- d[[10]][[2]]
  slope0 <- a_func*b_func*exp(0)
  print(d) 
  print(slope0) #calculate initial digestion rate
  summary(func)["r.squared"] 
  #rate constant k=-b
  
  #plot the curve fitting
  data_df.plot <- ggplot(data=data_df,
                         aes(x,y)) +
    geom_point() +
    xlab("Time (min)") +
    ylab("Relative band intensity (%)") +
    ggtitle(title) +
    theme_pubr(base_size = 15, base_family = "Times New Roman",  border = T) +
   
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(color = "grey20", size = 15), 
          axis.text.y = element_text(color = "grey20", size = 15))
  
  data_df.plot + stat_smooth(data = subset(data_df, y > -1), fullrange = T, method = "nls",
                             method.args = list(formula= y ~ a * exp(b * x),
                                                start = list(a = 66.88, b = -0.026273)), #initial value a = 66.88, b = -0.026273
                             se= F, size=1,
                             color="red", show.legend = T) +
    coord_cartesian(xlim=c(0,120),ylim=c(0, 120))
  
}

#For digestion products
Func_increase <- function(data_df, title){
  func = nlsLM(y ~ a- a * exp(b * x),data = subset(data_df),algorithm = "LM",start = list(a = 126, b = -0.1231)) #a = 5.621, b = -0.021
  d <- summary(func)
  
  a_func <- d[[10]][[1]]
  b_func <- d[[10]][[2]]
  slope0 <- -a_func*b_func*exp(0)
  intercept0 <- 0
  print(d)
  print(slope0) #calculate initial production rate
  summary(func)["r.squared"] 
  #rate constant k=-b
  
  #plot the curve fitting
  data_df.plot <- ggplot(data=data_df,
                         aes(x,y))+
    geom_point()+ 
    xlab("Timepoints (min)")+
    ylab("Relative abundance (%)")+
    geom_point(data = data_df, color = "black") + 
    ggtitle(title) +
    theme_pubr(base_size = 14, base_family = "Times New Roman",  border = T) +
    
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(color = "grey20", size = 12), 
          axis.text.y = element_text(color = "grey20", size = 12))
  
  data_df.plot+stat_smooth(data = subset(data_df), fullrange = T, method = "nls",
                           method.args = list(formula= y ~ a-a * exp(b * x),
                                              start = list(a = 126, b = -0.1231)), #a = 41.25, b = -11.54
                           se= F, size=1,
                           color="red", show.legend = T)+
    coord_cartesian(xlim=c(0,120),ylim=c(0, 120))+
    geom_abline(size = 1, intercept=intercept0, slope= slope0, col = "blue")
  
}

