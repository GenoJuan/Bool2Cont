source("Bool2Cont.r")
library(BoolNet)
library(ggplot2)
library(gridExtra)
library(reshape)

#########################################################################################################################
#########################################################################################################################
# Developed by: Juan A. Arias Del Angel (PhD Candidate) - Instituto de Ecología, Universidad Nacional Autónoma de México
# Advisor: Mariana Benítez (PhD) - Instituto de Ecología, Universidad Nacional Autónoma de México
# Last modificacion: May 12th 2018

# Description: It is shown a comparison between a full Boolean network model and a simplified version of itself.
# The Boolean network model (14 nodes) corresponds to the cell cycle in Arabidopsis thaliana.
# The model is taken from Ortíz-Gutierrez et al., 2016 (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004486).
# The state space converges to a single oscillatory attractor (11 states) under synchronous update scheme.
# The simplified version (8 nodes) preserves the oscillatory attractor (7 states) under synchronous upodate scheme as well as in a Bool-to-continuous mapping strategy.
# In the simplified version, the system converges to a complex attractor under asynchronous update scheme.
# This result is easily obtained in the simplified version in contrast to the full version.
# Note that the simplified version is obtained through an heuristic search so it is not guaranted that this is the optimal simplified version.
# The strategy was taken from Naldi et al., 2011 (https://www.sciencedirect.com/science/article/pii/S0304397510005839)
# The Bool2Cont was performed using the strategy presented in Villarreal et al., 2012 (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.118102)
#########################################################################################################################
#########################################################################################################################

dec2bin = function(dec,len){
    bin = rep(0,len)
    
    bin = .C("dec2bin",as.integer(bin),as.integer(dec),as.integer(len),NAOK=TRUE)[[1]]
    return(bin)
}

Attr2df = function(x){
    attrs = data.frame(Attr = c(), Variable = c(), Time = c(), State = c())
    genes = x$stateInfo$genes
    for(attr in 1:length(x$attractors)){
        for(s in 1:length(x$attractors[[attr]]$involvedStates)){
            attr_info = data.frame(Attr = rep(attr, length(genes)), Variable = genes, Time = rep(s, length(genes)), State = dec2bin(x$attractors[[attr]]$involvedStates[s], length(genes)))
            attrs = rbind.data.frame(attrs, attr_info)
        }        
    }
    return(attrs)
}

#################################################################################################################
#################################################################################################################
net = loadNetwork("Arabidopsis_thaliana_CC.net")
attr = getAttractors(net)
attr = Attr2df(attr)

Boolean_whole = ggplot(attr, aes(x = factor(Time), y = Variable, fill = factor(State)), colour = "black") + 
    geom_tile(width=0.9, height=0.9, size=3.0) +
    scale_fill_manual(values = c("#8DD35F", "#FF5555")) +
    #facet_grid(.~Attr, scales="free", space = "free") +  # Activate in case of multiple attractors
    theme_classic() +
    xlab("") +
    theme(axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    labs(fill = "State") +
    theme(legend.title = element_text(face = "bold", size = 9.5)) +
    ggtitle("a) Boolean dynamic for whole A. thaliana cell cycle model") +
    theme(plot.title = element_text(size = 9, face = "bold"))

net = loadNetwork("Arabidopsis_thaliana_CC_simplified.net")
attr = getAttractors(net)
attr = Attr2df(attr)

Boolean_simplified = ggplot(attr, aes(x = factor(Time), y = Variable, fill = factor(State)), colour = "black") + 
    geom_tile(width=0.9, height=0.9, size=3.0) +
    scale_fill_manual(values = c("#8DD35F", "#FF5555")) +
    # facet_grid(.~Attr, scales="free", space = "free") +  # Activate in case of multiple attractors
    theme_classic() +
    xlab("") +
    theme(axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    labs(fill = "State") +
    theme(legend.title = element_text(face = "bold", size = 9.5)) +
    ggtitle("c) Boolean dynamic for simplified A. thaliana cell cycle model") +
    theme(plot.title = element_text(size = 9, face = "bold"))


odeNet = BoolNet2ToContinuous("Arabidopsis_thaliana_CC.net")
times = seq(0, 35, 0.05)
o = ode(y = state, func = odeNet, parms = parameters, times = times)

o = as.data.frame(o)
o = melt(o, id = c("time"))
cont_whole = ggplot(subset(o, variable %in% c("KRP1", "CYCD31", "CYCB11")), aes(x = time, y = value, colour = variable)) + 
  geom_line() + 
  theme_linedraw() +
  ylab("Relative concentration") +
  xlab("Time (a.u.)") +
  theme(legend.title = element_text(face = "bold", size = 9.5)) +
  labs(colour = "Node", linetype = "Node") +
  ylim(0,1) + 
  scale_colour_brewer(palette="Set1") +
  ggtitle("b) Continuous dynamic for A. thaliana cell cycle model") +
  theme(plot.title = element_text(size = 9, face = "bold")) +
  theme(axis.title = element_text(size=7.5)) #+
  #facet_grid(variable~.)

odeNet = BoolNet2ToContinuous("Arabidopsis_thaliana_CC_simplified.net")
times = seq(0, 35, 0.05)
o = ode(y = state, func = odeNet, parms = parameters, times = times)


o = as.data.frame(o)
o = melt(o, id = c("time"))
cont_simplified = ggplot(subset(o, variable %in% c("CYCD31", "CYCB11", "KRP1")), aes(x = time, y = value, colour = variable)) + 
  geom_line() + 
  theme_linedraw() +
  ylab("Relative concentration") +
  xlab("Time (a.u.)") +
  theme(legend.title = element_text(face = "bold", size = 9.5)) +
  labs(colour = "Node", linetype = "Node") +
  ylim(0,1) + scale_colour_brewer(palette="Set1") +
  ggtitle("d) Continuous dynamic for simplified A. thaliana cell cycle model") +
  theme(plot.title = element_text(size = 9, face = "bold")) +
  theme(axis.title = element_text(size=7.5)) 
  
  
grid.arrange(arrangeGrob(Boolean_whole + theme(legend.position="none"),
                         Boolean_simplified + theme(legend.position="bottom"),
                         nrow = 2, heights = c(1.75,1.25)), 
             arrangeGrob(cont_whole + theme(legend.position="none"),
                         cont_simplified + theme(legend.position="bottom"),
                         nrow = 2, heights = c(1, 1.25)),
                         ncol = 2)
